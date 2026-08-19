#!/usr/bin/env python3
"""
pp_visualizer.py -- timing-diagram generator for X-Pulse pulse programmes.

Runs a real sequences/**/*.py file's run() through an INSTRUMENTED mock
firebird package (built on the same pattern as the mock hardware harness
in the xpulse-pulse-programmes skill / pulse-programme-parameters skill),
recording every hardware call as a timed event, then renders the result
as a self-contained HTML timing diagram with one lane per channel (H/F,
X, and one lane per gradient axis actually used).

Callable API (for interactive use -- REPL, notebook, another script):

    from pp_visualizer import visualize
    result = visualize("sequences/diffusion/LWG_PGSTE_H.py",
                        overrides={"DELTA": 20000}, open_browser=True)
    # result = {'svg':..., 'html':..., 'events':[...], 'total_time_us':...,
    #           'sequence_name':..., 'output_path':...}

CLI:
    python3 tools/pp_visualizer.py sequences/diffusion/LWG_PGSTE_H.py
    python3 tools/pp_visualizer.py sequences/diffusion/LWG_PGSTE_H.py -o out.html
    python3 tools/pp_visualizer.py sequences/relaxation/LWG_CPMG_H.py --set EchoNumber=4 --open

Only ONE scan's worth of timing is traced (NumScans/DS are forced to
1/0 before run() is called) -- the per-scan structure is what's
interesting to visualise; the outer scan-repeat loop is not.

Design notes on the timing model:
  - Hardware instruction "duration" arguments count as literal elapsed
    wall-clock time (per this repo's own timing-model convention -- see
    the xpulse-pulse-programmes skill sec.4) -- Delay() and every
    duration-arg call advance a single logical clock.
  - `with parallel:` blocks contain multiple `with sequential:`
    sub-branches that Python executes one after another in REAL source
    order (there is no actual concurrency in this trace -- "parallel" is
    a hardware compilation directive, not a Python runtime property).
    The instrumented sequential/parallel context managers below
    reproduce the INTENDED hardware timing anyway: entering a
    `sequential` block that is the immediate child of a `parallel` block
    rewinds the clock back to the parallel block's start time; exiting a
    parallel block advances the clock by the MAX of its branches'
    elapsed time (not the sum) -- matching how this repo's own
    time_calculation() functions already reason about parallel blocks.
  - Long idle gaps (chiefly the relaxation delay, RD, which is typically
    500,000-2,000,000us vs. low-hundreds-of-us pulses/gradients) would
    make a literal linear time axis useless -- actual RF/gradient/ACQU
    events are rendered TRUE TO SCALE relative to each other, but idle
    gaps between them are compressed on a log scale (still monotonic --
    a longer gap always gets modestly more width -- and always
    exactly-labelled with its real duration) so the interesting
    structure and the RD delay can both be seen on one diagram.
  - GRADIENT PULSES ARE TRAPEZOIDS. A real gradient event is
    ramp-up -- flat plateau (held for the "delta"/GradientOnTime
    duration) -- ramp-down. In the pp source this looks like
    Gradient3(RampTime, G1); Delay(delta); Gradient3(RampTime, 0) -- the
    plateau is just a generic Delay(), invisible to a naive tracer. The
    mock below tracks each gradient axis's held value explicitly and
    synthesises an intermediate 'plateau' event spanning that Delay(),
    then the renderer merges contiguous ramp/plateau events on one axis
    into a single filled trapezoid polygon (with a delta<->delta arrow
    labelling the hold time), instead of showing two disconnected
    diagonal lines with a gap in between.

Known limitations (v1):
  - Shaped pulses (shaped_pulse()'s `with parallel:` sub-branches of
    per-microsecond Transmit1SetScale()/Channel1SetBasePhase() calls)
    are rendered as a single flat pulse block, not with their true
    amplitude/phase envelope -- accurate for every HARD-pulse sequence
    in this repo (PGSTE/PGSE/CPMG/Hahn-echo/etc.), a simplification for
    the shaped/selective-pulse family (LWG_Selective-Echo_H.py etc.).
  - Only traces one code path through run() -- e.g. an `if int(P.WetOn)`
    branch is traced with whatever WetOn value the Parameters carry
    (default, or your --set/overrides value); it does not show both
    branches.
"""

import argparse
import math
import os
import runpy
import sys
import types
import webbrowser


# =============================================================================
# Clock: tracks elapsed time, handling with sequential:/with parallel:
# exactly as this repo's own time_calculation() functions reason about
# them (parallel block's own elapsed time = MAX of its branches').
# =============================================================================

class Clock(object):
    def __init__(self):
        self.t = 0.0
        self.stack = []  # dicts: {kind: 'sequential'|'parallel', start, max_elapsed}

    def advance(self, dt):
        if dt:
            self.t += float(dt)
        return self.t

    def enter_sequential(self):
        if self.stack and self.stack[-1]['kind'] == 'parallel':
            self.t = self.stack[-1]['start']
        self.stack.append({'kind': 'sequential', 'start': self.t})

    def exit_sequential(self):
        ctx = self.stack.pop()
        elapsed = self.t - ctx['start']
        if self.stack and self.stack[-1]['kind'] == 'parallel':
            self.stack[-1]['max_elapsed'] = max(self.stack[-1]['max_elapsed'], elapsed)

    def enter_parallel(self):
        self.stack.append({'kind': 'parallel', 'start': self.t, 'max_elapsed': 0.0})

    def exit_parallel(self):
        ctx = self.stack.pop()
        self.t = ctx['start'] + ctx['max_elapsed']


# =============================================================================
# Instrumented mock firebird package -- installed directly into
# sys.modules (no files written to disk) before runpy.run_path()'ing the
# target pp file.
# =============================================================================

def install_mock_firebird(clock, events):
    """Build and register sys.modules['firebird']/['firebird.applications'],
    wired to append timed events to `events` (a list of dicts) and advance
    `clock`. Returns the mock's mutable state dict (mostly useful for
    debugging); side-effects sys.modules."""

    state = {
        'phase': {1: 0.0, 2: 0.0},         # Channel1/2SetBasePhase
        'scale': {1: 1.0, 2: 1.0},         # Transmit1/2SetScale
        'rx_phase': {1: 0.0, 2: 0.0},      # Receiver1/2Phase
        'grad': {1: 0.0, 2: 0.0, 3: 0.0},        # current Gradient1/2/3 value
        'grad_since': {1: 0.0, 2: 0.0, 3: 0.0},  # time current value began
    }

    def emit(lane, kind, t_start, duration, **meta):
        events.append({
            'lane': lane, 'kind': kind,
            't_start': t_start, 't_end': t_start + duration,
            'duration': duration, 'meta': meta,
        })

    fb = types.ModuleType('firebird')
    fbapp = types.ModuleType('firebird.applications')

    # ---- Delay ----
    def Delay(value):
        clock.advance(value)
    fb.Delay = Delay

    # ---- Channel (RF frequency/phase) ----
    def _make_channel_freq(ch):
        def fn(duration, freq):
            t0 = clock.advance(duration) - duration
            emit('HF' if ch == 1 else 'X', 'config', t0, duration, label='SetFrequency', freq=freq)
        return fn

    def _make_channel_restart(ch):
        def fn(duration):
            clock.advance(duration)
        return fn

    def _make_channel_phase(ch):
        def fn(duration, phase):
            state['phase'][ch] = phase
            clock.advance(duration)
        return fn

    fb.Channel1SetFrequency = _make_channel_freq(1)
    fb.Channel2SetFrequency = _make_channel_freq(2)
    fb.Channel1RestartSynth = _make_channel_restart(1)
    fb.Channel2RestartSynth = _make_channel_restart(2)
    fb.Channel1SetBasePhase = _make_channel_phase(1)
    fb.Channel2SetBasePhase = _make_channel_phase(2)

    # ---- Transmit (RF pulses) ----
    def _make_transmit(ch):
        def fn(width):
            t0 = clock.advance(width) - width
            emit('HF' if ch == 1 else 'X', 'pulse', t0, width,
                 phase=state['phase'][ch], scale=state['scale'][ch])
        return fn

    def _make_set_scale(ch):
        def fn(duration, scale):
            state['scale'][ch] = scale
            clock.advance(duration)
        return fn

    def _make_blanking(ch, on):
        def fn(duration):
            clock.advance(duration)
        return fn

    def _make_select_port(ch):
        def fn(duration, port):
            clock.advance(duration)
        return fn

    def _make_lp_enable(ch):
        def fn(duration, enable):
            clock.advance(duration)
        return fn

    def _make_shaped_pulse(ch):
        def fn(width, phase, shapename):
            t0 = clock.advance(width) - width
            emit('HF' if ch == 1 else 'X', 'pulse', t0, width,
                 phase=phase, scale=state['scale'][ch], label='shaped:' + str(shapename))
        return fn

    fb.Transmit1 = _make_transmit(1)
    fb.Transmit2 = _make_transmit(2)
    fb.Transmit1SetScale = _make_set_scale(1)
    fb.Transmit2SetScale = _make_set_scale(2)
    fb.Transmit1BlankingOn = _make_blanking(1, True)
    fb.Transmit1BlankingOff = _make_blanking(1, False)
    fb.Transmit2BlankingOn = _make_blanking(2, True)
    fb.Transmit2BlankingOff = _make_blanking(2, False)
    fb.Transmit1SelectPort = _make_select_port(1)
    fb.Transmit2SelectPort = _make_select_port(2)
    fb.Transmit1LPEnable = _make_lp_enable(1)
    fb.Transmit2LPEnable = _make_lp_enable(2)
    fb.Transmit1ShapedPulse = _make_shaped_pulse(1)
    fb.Transmit2ShapedPulse = _make_shaped_pulse(2)

    # ---- Receiver ----
    def _make_rx_config(ch):
        def fn(duration, arg=None):
            clock.advance(duration)
        return fn

    def _make_rx_phase(ch):
        def fn(duration, phase):
            state['rx_phase'][ch] = phase
            clock.advance(duration)
        return fn

    def _make_receiver(ch):
        def fn(duration, points):
            t0 = clock.advance(duration) - duration
            emit('HF' if ch == 1 else 'X', 'acqu', t0, duration,
                 phase=state['rx_phase'][ch], points=points)
        return fn

    fb.Receiver1Preamp = _make_rx_config(1)
    fb.Receiver2Preamp = _make_rx_config(2)
    fb.Receiver1Filter = _make_rx_config(1)
    fb.Receiver2Filter = _make_rx_config(2)
    fb.Receiver1FilterFlush = _make_rx_config(1)
    fb.Receiver2FilterFlush = _make_rx_config(2)
    fb.Receiver1Phase = _make_rx_phase(1)
    fb.Receiver2Phase = _make_rx_phase(2)
    fb.Receiver1 = _make_receiver(1)
    fb.Receiver2 = _make_receiver(2)

    # ---- Gradients: track held value per axis, synthesise a 'plateau'
    # event spanning whatever Delay() happens while a nonzero value is
    # held, so ramp+plateau+ramp can be rendered as one true trapezoid. ----
    def _make_gradient(axis):
        lane = 'GRAD{0}'.format(axis)

        def fn(ramp_time, value):
            t_ramp_start = clock.t
            v0 = state['grad'][axis]
            if value == v0:
                # Degenerate/safety call (e.g. "Gradient3(1,0)" when
                # already at 0) -- just advance the clock, no event, and
                # don't disturb the held-since bookkeeping.
                clock.advance(ramp_time)
                return
            since = state['grad_since'][axis]
            if v0 != 0 and t_ramp_start > since:
                emit(lane, 'plateau', since, t_ramp_start - since, value=v0)
            clock.advance(ramp_time)
            t_ramp_end = clock.t
            emit(lane, 'ramp', t_ramp_start, ramp_time, value_from=v0, value_to=value)
            state['grad'][axis] = value
            state['grad_since'][axis] = t_ramp_end
        return fn

    def _make_gradient_slew(axis):
        def fn(duration, rate):
            clock.advance(duration)
        return fn

    fb.Gradient1 = _make_gradient(1)
    fb.Gradient2 = _make_gradient(2)
    fb.Gradient3 = _make_gradient(3)
    fb.Gradient1SlewRate = _make_gradient_slew(1)
    fb.Gradient2SlewRate = _make_gradient_slew(2)
    fb.Gradient3SlewRate = _make_gradient_slew(3)

    def GradientMatrix(duration, matrix):
        clock.advance(duration)
    fb.GradientMatrix = GradientMatrix

    # ---- Triggers ----
    def _make_trigger(n):
        def fn():
            emit('TRIGGER', 'marker', clock.t, 0.0, label='ExternalTrigger{0}'.format(n))
        return fn
    fb.ExternalTrigger1 = _make_trigger(1)
    fb.ExternalTrigger2 = _make_trigger(2)
    fb.ExternalTrigger3 = _make_trigger(3)

    # ---- End-of-scan housekeeping (no timing significance for one scan) ----
    fb.WriteToHardware = lambda *a, **kw: None
    fb.start = lambda *a, **kw: None
    fb.get_single_scan_execution_time = lambda: 0
    fb.Get_Compilation_Time = lambda: 0

    class _TXStub(object):
        def setup_receive(self, *a, **kw):
            pass

        def wait_for_data(self, scan, total, callback):
            import numpy as np
            npoints = 16
            fake = (np.zeros(npoints, dtype=np.int32), type('C', (), {'value': 0})())
            callback(fake, scan)
    fb.TX0 = _TXStub()
    fb.TX1 = _TXStub()

    # ---- with sequential:/with parallel: ----
    class _Sequential(object):
        def __enter__(self):
            clock.enter_sequential()
            return self

        def __exit__(self, *exc):
            clock.exit_sequential()
            return False

    class _Parallel(object):
        def __enter__(self):
            clock.enter_parallel()
            return self

        def __exit__(self, *exc):
            clock.exit_parallel()
            return False

    fb.sequential = _Sequential()
    fb.parallel = _Parallel()

    def sequential_main(scan, total):
        return _Sequential()
    fb.sequential_main = sequential_main

    def repeat(n):
        return _Sequential()
    fb.repeat = repeat

    # ---- Phases ----
    class PhaseListContainer(object):
        def __init__(self, phase_str):
            self.phases = [float(x) for x in str(phase_str).split(",") if x.strip() != ""] or [0.0]
            self.idx = 0

        def Reset(self):
            self.idx = 0

        def Inc(self):
            v = self.phases[self.idx % len(self.phases)]
            self.idx += 1
            return v
    fb.PhaseListContainer = PhaseListContainer

    class PhasesManager(object):
        def __init__(self, P):
            self.P = P
            cls = P if isinstance(P, type) else type(P)
            self._containers = {}
            for name in dir(cls):
                if name.startswith("_"):
                    continue
                if not (name.startswith("PH") or name.endswith("Phase")):
                    continue
                value = getattr(cls, name)
                if not isinstance(value, str):
                    raise AttributeError(
                        "{0!r} object has no attribute 'split' -- Parameter "
                        "attribute '{1}' looks like a phase-cycle Parameter "
                        "but its value {2!r} is not a string.".format(
                            type(value).__name__, name, value))
                self._containers[name] = PhaseListContainer(value)

        def Reset(self):
            for c in self._containers.values():
                c.Reset()

        def Incd(self):
            return dict((name, c.Inc()) for name, c in self._containers.items())
    fb.PhasesManager = PhasesManager

    # ---- Filter / ParameterTypes / Parameter / ParameterBlock ----
    class Filter(object):
        def __init__(self, spec):
            self.dwell = 1.0
            self.dead_time = 10.0
            self.group_delay = 5.0
            self.gain = 1.0
    fb.Filter = Filter
    fb.ChooseFilter = lambda spec: spec

    class ParameterTypes(object):
        Double = float
        Int32 = int
        String = str
    fb.ParameterTypes = ParameterTypes

    class Parameter(object):
        def __init__(self, code, default, ptype, desc, validator=None, **kw):
            self.code, self.default = code, default

        def __get__(self, obj, objtype=None):
            return self.default

        def __set__(self, obj, value):
            self.default = value
    fb.Parameter = Parameter

    def ParameterBlock(cls):
        return cls
    fb.ParameterBlock = ParameterBlock

    sys.modules['firebird'] = fb
    sys.modules['firebird.applications'] = fbapp
    fb.applications = fbapp

    return state


class FakeComms(object):
    def log(self, msg):
        pass

    def send_data(self, *a, **kw):
        pass


# =============================================================================
# Tracing
# =============================================================================

def _apply_overrides(P, overrides):
    if not overrides:
        return
    for key, value in overrides.items():
        if not hasattr(P, key):
            raise ValueError("Unknown Parameter '{0}' (not found on this file's "
                              "Parameters class).".format(key))
        current = getattr(P, key)
        caster = type(current) if not isinstance(current, bool) else str
        try:
            setattr(P, key, caster(value))
        except (TypeError, ValueError):
            setattr(P, key, value)


def _load_parameters(pp_path, overrides=None):
    """runpy-load pp_path fresh (mock must already be installed in
    sys.modules) and return (module_dict, Parameters class), with
    NumScans/DS forced to 1/0 and overrides applied."""
    mod = runpy.run_path(pp_path)
    P = mod['Parameters']
    P.NumScans = 1
    P.DS = 0
    _apply_overrides(P, overrides)
    return mod, P


def dump_params(pp_path, overrides=None):
    """Return a name->value dict of every Parameter default on pp_path,
    with overrides applied -- does NOT need the mock installed first
    (safe to call standalone)."""
    clock = Clock()
    install_mock_firebird(clock, [])
    mod, P = _load_parameters(pp_path, overrides)
    cls = type(P)
    out = {}
    for name in sorted(dir(cls)):
        if name.startswith('_'):
            continue
        val = getattr(P, name)
        if callable(val):
            continue
        out[name] = val
    return out


def trace_pp_file(pp_path, overrides=None):
    """Run one scan of pp_path's run() through the instrumented mock and
    return (events, total_time_us, sequence_name)."""
    clock = Clock()
    events = []
    install_mock_firebird(clock, events)

    mod, P = _load_parameters(pp_path, overrides)
    mod['run'](FakeComms())

    seq_name = getattr(P, 'Sequence', os.path.basename(pp_path))
    return events, clock.t, seq_name


def merge_gradient_trapezoids(events):
    """Merge contiguous ramp/plateau events on each GRAD lane into single
    filled-trapezoid composite events (kind='trapezoid'), tracing the
    actual v(t) path of the gradient coil -- ramp up, flat top for the
    held ("delta") duration, ramp down -- instead of two disconnected
    diagonal 'ramp' lines with an invisible gap where the plateau was."""
    by_lane = {}
    passthrough = []
    for e in events:
        if e['lane'].startswith('GRAD') and e['kind'] in ('ramp', 'plateau'):
            by_lane.setdefault(e['lane'], []).append(e)
        else:
            passthrough.append(e)

    trapezoids = []
    for lane, evs in by_lane.items():
        evs = sorted(evs, key=lambda e: e['t_start'])
        run = []
        for e in evs:
            if run and abs(e['t_start'] - run[-1]['t_end']) < 1e-6:
                run.append(e)
            else:
                if run:
                    trapezoids.append(_trapezoid_from_run(lane, run))
                run = [e]
        if run:
            trapezoids.append(_trapezoid_from_run(lane, run))
    return passthrough + trapezoids


def _trapezoid_from_run(lane, run):
    verts = []
    for e in run:
        if e['kind'] == 'ramp':
            if not verts:
                verts.append((e['t_start'], e['meta']['value_from']))
            verts.append((e['t_end'], e['meta']['value_to']))
        else:  # plateau
            v = e['meta']['value']
            if not verts:
                verts.append((e['t_start'], v))
            verts.append((e['t_end'], v))
    t_start, t_end = run[0]['t_start'], run[-1]['t_end']
    peak = max((abs(v) for _, v in verts), default=0.0)
    plateaus = [(e['t_start'], e['t_end'], e['meta']['value'])
                for e in run if e['kind'] == 'plateau' and e['duration'] > 0]
    return {
        'lane': lane, 'kind': 'trapezoid',
        't_start': t_start, 't_end': t_end, 'duration': t_end - t_start,
        'meta': {'vertices': verts, 'peak': peak, 'plateaus': plateaus},
    }


# =============================================================================
# Rendering (HTML/SVG)
# =============================================================================

LANE_ORDER_HINT = ['HF', 'X', 'GRAD1', 'GRAD2', 'GRAD3', 'TRIGGER']
LANE_LABELS = {
    'HF': 'H/F',
    'X': 'X',
    'GRAD1': 'Grad 1 (x)',
    'GRAD2': 'Grad 2 (y)',
    'GRAD3': 'Grad 3 (z)',
    'TRIGGER': 'Trigger',
}

EVENT_PX_PER_US = 0.10          # true-to-scale width for actual events
GAP_BASE_PX = 16.0               # minimum width for any idle gap
GAP_LOG_PX = 14.0                # log-scale coefficient for idle gap width
GAP_COMPRESS_THRESHOLD = 60.0    # gaps shorter than this render at true scale
LANE_HEIGHT = 52
LANE_GAP = 18
LEFT_MARGIN = 110
TOP_MARGIN = 78
RIGHT_MARGIN = 40
BOTTOM_MARGIN = 30

# A calm, colorblind-considerate categorical palette (dark-theme-first;
# same hues carry into the light-theme CSS variables below).
COLORS = {
    'pulse': '#6C8EF5',
    'pulse_text': '#0b1220',
    'acqu': '#F0955C',
    'grad': '#3FC08A',
    'config': '#8A93A6',
    'marker': '#E0B84B',
}


def _fmt_us(v):
    v = float(v)
    if v >= 1.0e6:
        return "{0:.3f} s".format(v / 1.0e6)
    if v >= 1000.0:
        return "{0:.3f} ms".format(v / 1000.0)
    return "{0:.1f} µs".format(v)


def _esc(s):
    return (str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            .replace('"', '&quot;'))


class TimeAxis(object):
    """Piecewise time->pixel mapping: actual event durations render true
    to scale; idle gaps beyond GAP_COMPRESS_THRESHOLD render log-
    compressed (still monotonic, always exactly labelled)."""

    def __init__(self, events, total_time):
        occupied = sorted((e['t_start'], e['t_end']) for e in events if e['duration'] > 0)
        merged = []
        for s, en in occupied:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], en))
            else:
                merged.append((s, en))

        segments = []
        cursor = 0.0
        for s, en in merged:
            if s > cursor:
                segments.append((cursor, s, False))
            segments.append((s, en, True))
            cursor = en
        if cursor < total_time:
            segments.append((cursor, total_time, False))

        self.px_map = []
        px_cursor = 0.0
        for t0, t1, is_event in segments:
            dur = t1 - t0
            if dur <= 0:
                continue
            if is_event or dur <= GAP_COMPRESS_THRESHOLD:
                width = dur * EVENT_PX_PER_US
            else:
                width = GAP_BASE_PX + GAP_LOG_PX * math.log10(1.0 + dur)
            self.px_map.append((t0, t1, px_cursor, px_cursor + width))
            px_cursor += width
        self.total_px = px_cursor

    def px(self, t):
        for (a, b, pa, pb) in self.px_map:
            if a <= t <= b + 1e-9:
                return pa if b <= a else pa + (pb - pa) * (t - a) / (b - a)
        return self.total_px

    def span(self, t0, t1):
        x0, x1 = self.px(t0), self.px(t1)
        return x0, max(x1 - x0, 0.6)

    def gaps(self):
        return [(a, b, pa, pb) for (a, b, pa, pb) in self.px_map if (b - a) > GAP_COMPRESS_THRESHOLD]


def _lane_y_for_value(v, lane_top, lane_bot):
    v = max(-1.0, min(1.0, float(v)))
    return lane_bot - (v + 1.0) / 2.0 * (lane_bot - lane_top)


def render_svg(events, total_time, seq_name, src_name):
    events = merge_gradient_trapezoids(events)

    active_lanes = []
    seen = set()
    for hint in LANE_ORDER_HINT:
        if any(e['lane'] == hint for e in events) and hint not in seen:
            active_lanes.append(hint)
            seen.add(hint)
    for e in events:
        if e['lane'] not in seen:
            active_lanes.append(e['lane'])
            seen.add(e['lane'])
    if not active_lanes:
        active_lanes = ['HF']

    axis = TimeAxis(events, total_time)
    width = LEFT_MARGIN + axis.total_px + RIGHT_MARGIN
    height = TOP_MARGIN + len(active_lanes) * (LANE_HEIGHT + LANE_GAP) + BOTTOM_MARGIN

    svg = []
    svg.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {0:.0f} {1:.0f}" '
                'width="100%" style="min-width:{0:.0f}px" '
                'font-family="-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif">'
                .format(width, height))
    svg.append('<rect x="0" y="0" width="{0:.0f}" height="{1:.0f}" fill="var(--pv-bg)"/>'
                .format(width, height))

    svg.append('<text x="{0}" y="26" font-size="16" font-weight="600" fill="var(--pv-fg)">{1}</text>'
                .format(LEFT_MARGIN, _esc(seq_name)))
    svg.append('<text x="{0}" y="44" font-size="11" fill="var(--pv-muted)">{1} — total {2} '
                '(one scan; idle gaps &gt;{3:.0f}µs log-compressed, hover any element for exact values)</text>'
                .format(LEFT_MARGIN, _esc(src_name), _fmt_us(total_time), GAP_COMPRESS_THRESHOLD))

    # Legend
    legend = [('pulse', 'RF pulse'), ('acqu', 'Acquisition'), ('grad', 'Gradient'), ('marker', 'Trigger')]
    lx = LEFT_MARGIN
    for kind, label in legend:
        svg.append('<rect x="{0}" y="54" width="10" height="10" rx="2" fill="{1}"/>'
                    .format(lx, COLORS[kind]))
        svg.append('<text x="{0}" y="63" font-size="10" fill="var(--pv-muted)">{1}</text>'
                    .format(lx + 14, _esc(label)))
        lx += 20 + 7 * len(label)

    # Gap break markers + duration labels along the top
    for (a, b, pa, pb) in axis.gaps():
        cx = LEFT_MARGIN + (pa + pb) / 2.0
        svg.append('<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="var(--pv-grid)" '
                    'stroke-width="1" stroke-dasharray="2,3"/>'
                    .format(cx, TOP_MARGIN - 14, height - BOTTOM_MARGIN))
        svg.append('<text x="{0:.1f}" y="{1}" font-size="10" fill="var(--pv-muted)" '
                    'text-anchor="middle">{2}</text>'
                    .format(cx, TOP_MARGIN - 18, _fmt_us(b - a)))

    for i, lane in enumerate(active_lanes):
        y0 = TOP_MARGIN + i * (LANE_HEIGHT + LANE_GAP)
        y_mid = y0 + LANE_HEIGHT / 2.0
        is_grad = lane.startswith('GRAD')
        label = LANE_LABELS.get(lane, lane)

        svg.append('<rect x="0" y="{0:.1f}" width="{1:.0f}" height="{2}" fill="var(--pv-lane-bg)"/>'
                    .format(y0 - 4, width, LANE_HEIGHT + 8))
        svg.append('<text x="{0}" y="{1:.1f}" font-size="12" font-weight="600" fill="var(--pv-fg)" '
                    'text-anchor="end" dominant-baseline="middle">{2}</text>'
                    .format(LEFT_MARGIN - 14, y_mid, _esc(label)))

        if is_grad:
            lane_top, lane_bot = y0 + 6, y0 + LANE_HEIGHT - 6
            for tick_v, tick_label in ((1.0, '+1'), (0.0, '0'), (-1.0, '−1')):
                ty = _lane_y_for_value(tick_v, lane_top, lane_bot)
                svg.append('<line x1="{0}" y1="{1:.1f}" x2="{2:.0f}" y2="{1:.1f}" '
                            'stroke="var(--pv-grid)" stroke-width="1" '
                            'stroke-dasharray="{3}"/>'
                            .format(LEFT_MARGIN, ty, width - RIGHT_MARGIN,
                                    "1,0" if tick_v == 0 else "2,4"))
                svg.append('<text x="{0}" y="{1:.1f}" font-size="9" fill="var(--pv-muted)" '
                            'text-anchor="end" dominant-baseline="middle">{2}</text>'
                            .format(LEFT_MARGIN - 3, ty, tick_label))
        else:
            svg.append('<line x1="{0}" y1="{1:.1f}" x2="{2:.0f}" y2="{1:.1f}" '
                        'stroke="var(--pv-grid)" stroke-width="1"/>'
                        .format(LEFT_MARGIN, y_mid, width - RIGHT_MARGIN))

        for e in events:
            if e['lane'] != lane:
                continue
            _render_event(svg, e, axis, y0, y_mid)

    svg.append('</svg>')
    return "\n".join(svg)


def _render_event(svg, e, axis, y0, y_mid):
    x, w = axis.span(e['t_start'], e['t_end'])
    x += LEFT_MARGIN
    kind = e['kind']
    title_bits = ["t={0}".format(_fmt_us(e['t_start'])), "dur={0}".format(_fmt_us(e['duration']))]
    for k, v in e['meta'].items():
        if k in ('vertices', 'plateaus'):
            continue
        title_bits.append("{0}={1}".format(k, v))
    title = _esc(" ".join(title_bits))

    if kind == 'trapezoid':
        lane_top, lane_bot = y0 + 6, y0 + LANE_HEIGHT - 6
        baseline_y = _lane_y_for_value(0, lane_top, lane_bot)
        verts = e['meta']['vertices']
        pts = [(axis.px(verts[0][0]) + LEFT_MARGIN, baseline_y)]
        for (t, v) in verts:
            pts.append((axis.px(t) + LEFT_MARGIN, _lane_y_for_value(v, lane_top, lane_bot)))
        pts.append((axis.px(verts[-1][0]) + LEFT_MARGIN, baseline_y))
        pts_str = " ".join("{0:.1f},{1:.1f}".format(px, py) for px, py in pts)
        svg.append('<polygon points="{0}" fill="{1}" fill-opacity="0.55" '
                    'stroke="{1}" stroke-width="1.5" stroke-linejoin="round">'
                    '<title>{2}</title></polygon>'.format(pts_str, COLORS['grad'], title))
        for (pt0, pt1, pv) in e['meta']['plateaus']:
            px0, pw = axis.span(pt0, pt1)
            px0 += LEFT_MARGIN
            if pw < 14:
                continue
            arrow_y = _lane_y_for_value(pv, lane_top, lane_bot) + (-10 if pv >= 0 else 10)
            svg.append('<line x1="{0:.1f}" y1="{1:.1f}" x2="{2:.1f}" y2="{1:.1f}" '
                        'stroke="var(--pv-fg)" stroke-width="1" marker-start="url(#pv-arrow)" '
                        'marker-end="url(#pv-arrow)"/>'.format(px0, arrow_y, px0 + pw))
            if pw > 22:
                svg.append('<text x="{0:.1f}" y="{1:.1f}" font-size="9" fill="var(--pv-fg)" '
                            'text-anchor="middle">δ</text>'
                            .format(px0 + pw / 2.0, arrow_y - 4))
        return

    if kind == 'marker':
        svg.append('<circle cx="{0:.1f}" cy="{1:.1f}" r="4" fill="{2}"><title>{3}</title></circle>'
                    .format(x, y_mid, COLORS['marker'], title))
        return

    if kind == 'config':
        svg.append('<rect x="{0:.1f}" y="{1:.1f}" width="{2:.1f}" height="4" fill="{3}" '
                    'opacity="0.6"><title>{4}</title></rect>'
                    .format(x, y0 + LANE_HEIGHT - 4, max(w, 1.5), COLORS['config'], title))
        return

    # pulse / acqu
    ry, rh = y0 + 6, LANE_HEIGHT - 12
    color = COLORS['pulse'] if kind == 'pulse' else COLORS['acqu']
    style = 'fill="{0}"'.format(color) if kind == 'pulse' else \
        'fill="{0}" fill-opacity="0.35" stroke="{0}" stroke-width="1.5"'.format(color)
    svg.append('<rect x="{0:.1f}" y="{1:.1f}" width="{2:.1f}" height="{3:.1f}" rx="4" {4}>'
                '<title>{5}</title></rect>'
                .format(x, ry, max(w, 1.5), rh, style, title))
    if w > 26:
        label_txt = e['meta'].get('label')
        if not label_txt:
            if kind == 'pulse':
                label_txt = "φ={0:g}°".format(e['meta'].get('phase', 0))
            else:
                label_txt = "ACQU"
        text_color = COLORS['pulse_text'] if kind == 'pulse' else 'var(--pv-fg)'
        svg.append('<text x="{0:.1f}" y="{1:.1f}" font-size="9.5" '
                    'fill="{2}" text-anchor="middle" dominant-baseline="middle" '
                    'pointer-events="none">{3}</text>'
                    .format(x + w / 2.0, y_mid, text_color, _esc(label_txt)))


def render_html(svg, seq_name, pp_path, param_summary):
    return """<!doctype html>
<title>{title}</title>
<style>
  :root {{
    --pv-bg:#0f1117; --pv-fg:#e7e9ee; --pv-muted:#8b93a5;
    --pv-grid:#262b38; --pv-lane-bg:#161a24;
  }}
  @media (prefers-color-scheme: light) {{
    :root:not([data-theme="dark"]) {{
      --pv-bg:#ffffff; --pv-fg:#15181f; --pv-muted:#5b6472;
      --pv-grid:#e6e9f0; --pv-lane-bg:#f6f7fb;
    }}
  }}
  :root[data-theme="light"] {{
    --pv-bg:#ffffff; --pv-fg:#15181f; --pv-muted:#5b6472;
    --pv-grid:#e6e9f0; --pv-lane-bg:#f6f7fb;
  }}
  body {{ margin:0; background:var(--pv-bg); color:var(--pv-fg);
          font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif; }}
  .wrap {{ padding:20px; overflow-x:auto; }}
  .params {{ padding:0 20px 24px; font-family:ui-monospace,Menlo,monospace;
             font-size:12px; color:var(--pv-muted); white-space:pre-wrap; }}
  h2 {{ font-size:13px; color:var(--pv-muted); font-weight:600; margin:0 0 8px; }}
  svg text {{ user-select:none; }}
</style>
<svg width="0" height="0" style="position:absolute">
  <defs>
    <marker id="pv-arrow" viewBox="0 0 8 8" refX="4" refY="4" markerWidth="5" markerHeight="5" orient="auto-start-reverse">
      <path d="M0,0 L8,4 L0,8" fill="none" stroke="var(--pv-fg)" stroke-width="1.5"/>
    </marker>
  </defs>
</svg>
<div class="wrap">
{svg}
</div>
<div class="params">
<h2>Parameters used for this trace ({pp_path})</h2>{params}
</div>
""".format(title=_esc(seq_name), svg=svg, pp_path=_esc(pp_path), params=_esc(param_summary))


# =============================================================================
# Public callable API
# =============================================================================

def visualize(pp_path, overrides=None, output=None, open_browser=False):
    """Trace pp_path and render an HTML timing diagram.

    overrides: dict of {ParameterName: value} to apply before tracing
        (e.g. {"DELTA": 20000, "EchoNumber": 4}).
    output: path to write the HTML to. Defaults to
        '<pp_path without .py>.timing.html'. Pass output=False to skip
        writing to disk (useful in a notebook: just use the returned
        'html'/'svg' strings, e.g. IPython.display.HTML(result['html'])).
    open_browser: if True, opens the written file in the default browser.

    Returns a dict: {svg, html, events, total_time_us, sequence_name,
    output_path}. 'events' are the RAW traced events (before the
    gradient-trapezoid merge the renderer applies internally) -- useful
    if you want to build your own rendering (e.g. the TikZ generator in
    pp_tikz.py reuses this).
    """
    events, total_time, seq_name = trace_pp_file(pp_path, overrides)
    svg = render_svg(events, total_time, seq_name, os.path.basename(pp_path))
    params = dump_params(pp_path, overrides)
    param_summary = "\n".join("{0} = {1}".format(k, v) for k, v in params.items())
    html = render_html(svg, seq_name, pp_path, param_summary)

    out_path = None
    if output is not False:
        out_path = output or (os.path.splitext(pp_path)[0] + ".timing.html")
        with open(out_path, "w") as f:
            f.write(html)
        if open_browser:
            webbrowser.open("file://" + os.path.abspath(out_path))

    return {
        'svg': svg, 'html': html, 'events': events,
        'total_time_us': total_time, 'sequence_name': seq_name,
        'output_path': out_path,
    }


# =============================================================================
# CLI
# =============================================================================

def _parse_overrides(pairs):
    out = {}
    for p in pairs or []:
        if '=' not in p:
            raise SystemExit("--set expects KEY=VALUE, got: {0}".format(p))
        k, v = p.split('=', 1)
        out[k.strip()] = v.strip()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pp_file', help="Path to a sequences/**/*.py pulse programme")
    ap.add_argument('-o', '--output', help="Output HTML path (default: alongside input, .timing.html)")
    ap.add_argument('--set', action='append', dest='overrides', metavar='KEY=VALUE',
                     help="Override a Parameter before tracing (repeatable), e.g. --set EchoNumber=4")
    ap.add_argument('--open', action='store_true', help="Open the result in your default browser")
    args = ap.parse_args()

    if not os.path.isfile(args.pp_file):
        raise SystemExit("No such file: {0}".format(args.pp_file))

    overrides = _parse_overrides(args.overrides)
    result = visualize(args.pp_file, overrides=overrides, output=args.output, open_browser=args.open)
    print("Wrote {0} ({1} events, total {2})".format(
        result['output_path'], len(result['events']), _fmt_us(result['total_time_us'])))


if __name__ == '__main__':
    main()
