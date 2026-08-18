#!/usr/bin/env python3
"""
pp_visualizer.py -- timing-diagram generator for X-Pulse pulse programmes.

Runs a real sequences/**/*.py file's run() through an INSTRUMENTED mock
firebird package (built on the same pattern as the mock hardware harness
in the xpulse-pulse-programmes skill / pulse-programme-parameters skill),
recording every hardware call as a timed event, then renders the result
as a self-contained HTML timing diagram with one lane per channel (H/F,
X, and one lane per gradient axis actually used).

Usage:
    python3 tools/pp_visualizer.py sequences/diffusion/LWG_PGSTE_H.py
    python3 tools/pp_visualizer.py sequences/diffusion/LWG_PGSTE_H.py -o out.html
    python3 tools/pp_visualizer.py sequences/relaxation/LWG_CPMG_H.py --set EchoNumber=4 --set Tau=1000

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

Known limitations (v1):
  - Shaped pulses (shaped_pulse()'s `with parallel:` sub-branches of
    per-microsecond Transmit1SetScale()/Channel1SetBasePhase() calls)
    are rendered as a single flat pulse block, not with their true
    amplitude/phase envelope -- accurate for every HARD-pulse sequence
    in this repo (PGSTE/PGSE/CPMG/Hahn-echo/etc.), a simplification for
    the shaped/selective-pulse family (LWG_Selective-Echo_H.py etc.).
  - Only traces one code path through run() -- e.g. an `if int(P.WetOn)`
    branch is traced with whatever WetOn value the Parameters carry
    (default, or your --set override); it does not show both branches.
"""

import argparse
import json
import math
import os
import runpy
import sys
import types


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
    `clock`. Returns nothing; side-effects sys.modules."""

    state = {
        'phase': {1: 0.0, 2: 0.0},       # Channel1/2SetBasePhase
        'scale': {1: 1.0, 2: 1.0},       # Transmit1/2SetScale
        'rx_phase': {1: 0.0, 2: 0.0},    # Receiver1/2Phase
        'grad': {1: 0.0, 2: 0.0, 3: 0.0},  # current Gradient1/2/3 value
    }

    def emit(lane, kind, t_start, duration, **meta):
        events.append({
            'lane': lane, 'kind': kind,
            't_start': t_start, 't_end': t_start + duration,
            'duration': duration, 'meta': meta,
        })

    def _noop_advance(dt=0):
        def fn(*a, **kw):
            clock.advance(dt)
        return fn

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
                 phase=phase, scale=state['scale'][ch], label='shaped:'+str(shapename))
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
    def _make_rx_config(ch, label):
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

    fb.Receiver1Preamp = _make_rx_config(1, 'Preamp')
    fb.Receiver2Preamp = _make_rx_config(2, 'Preamp')
    fb.Receiver1Filter = _make_rx_config(1, 'Filter')
    fb.Receiver2Filter = _make_rx_config(2, 'Filter')
    fb.Receiver1FilterFlush = _make_rx_config(1, 'FilterFlush')
    fb.Receiver2FilterFlush = _make_rx_config(2, 'FilterFlush')
    fb.Receiver1Phase = _make_rx_phase(1)
    fb.Receiver2Phase = _make_rx_phase(2)
    fb.Receiver1 = _make_receiver(1)
    fb.Receiver2 = _make_receiver(2)

    # ---- Gradients ----
    def _make_gradient(axis):
        def fn(ramp_time, value):
            t0 = clock.advance(ramp_time) - ramp_time
            v0 = state['grad'][axis]
            state['grad'][axis] = value
            emit('GRAD{0}'.format(axis), 'ramp', t0, ramp_time, value_from=v0, value_to=value)
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

def trace_pp_file(pp_path, overrides=None):
    """Run one scan of pp_path's run() through the instrumented mock and
    return (events, total_time_us, sequence_name)."""
    clock = Clock()
    events = []
    install_mock_firebird(clock, events)

    # Force a clean re-import of the target module each call (avoid stale
    # sys.modules entries if this is invoked more than once in-process).
    mod_name = None
    for name in list(sys.modules):
        if name == 'Parameters':
            del sys.modules[name]

    mod = runpy.run_path(pp_path)
    P = mod['Parameters']

    # Only one scan's worth of timing is interesting.
    P.NumScans = 1
    P.DS = 0

    if overrides:
        for key, value in overrides.items():
            if not hasattr(P, key):
                raise SystemExit("Unknown Parameter '{0}' for --set (not found on "
                                  "this file's Parameters class).".format(key))
            current = getattr(P, key)
            caster = type(current) if not isinstance(current, bool) else str
            try:
                setattr(P, key, caster(value))
            except (TypeError, ValueError):
                setattr(P, key, value)

    mod['run'](FakeComms())

    seq_name = getattr(P, 'Sequence', os.path.basename(pp_path))
    return events, clock.t, seq_name


# =============================================================================
# Rendering
# =============================================================================

LANE_ORDER_HINT = ['HF', 'X', 'GRAD1', 'GRAD2', 'GRAD3', 'TRIGGER']
LANE_LABELS = {
    'HF': 'H/F channel',
    'X': 'X channel',
    'GRAD1': 'Gradient 1 (typ. X)',
    'GRAD2': 'Gradient 2 (typ. Y)',
    'GRAD3': 'Gradient 3 (typ. Z)',
    'TRIGGER': 'Trigger',
}

EVENT_PX_PER_US = 0.09      # true-to-scale width for actual events
GAP_BASE_PX = 14.0          # minimum width for any idle gap
GAP_LOG_PX = 13.0           # log-scale coefficient for idle gap width
GAP_COMPRESS_THRESHOLD = 60.0  # gaps shorter than this render at true scale
LANE_HEIGHT = 46
LANE_GAP = 14
LEFT_MARGIN = 170
TOP_MARGIN = 70
RIGHT_MARGIN = 40
BOTTOM_MARGIN = 60

LANE_COLORS = {
    'pulse': '#4C78E8',
    'acqu': '#E8674C',
    'ramp': '#3FAE6B',
    'config': '#B0B7C3',
    'marker': '#C9A227',
}


def _fmt_us(v):
    v = float(v)
    if v >= 1.0e6:
        return "{0:.3f} s".format(v / 1.0e6)
    if v >= 1000.0:
        return "{0:.3f} ms".format(v / 1000.0)
    return "{0:.1f} us".format(v)


def build_time_transform(events, total_time):
    """Return (px_for(t_start, duration) -> (x, width), total_px), using
    true-to-scale widths for event durations and log-compressed widths
    for idle gaps between them."""
    boundaries = set([0.0, total_time])
    for e in events:
        boundaries.add(e['t_start'])
        boundaries.add(e['t_end'])
    boundaries = sorted(boundaries)

    occupied = []
    for e in events:
        if e['duration'] > 0:
            occupied.append((e['t_start'], e['t_end']))
    occupied.sort()
    merged = []
    for s, en in occupied:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], en))
        else:
            merged.append((s, en))

    segments = []  # (t0, t1, is_event)
    cursor = 0.0
    for s, en in merged:
        if s > cursor:
            segments.append((cursor, s, False))
        segments.append((s, en, True))
        cursor = en
    if cursor < total_time:
        segments.append((cursor, total_time, False))

    px_map = []  # (t0, t1, px0, px1)
    px_cursor = 0.0
    for t0, t1, is_event in segments:
        dur = t1 - t0
        if dur <= 0:
            continue
        if is_event:
            width = dur * EVENT_PX_PER_US
        elif dur <= GAP_COMPRESS_THRESHOLD:
            width = dur * EVENT_PX_PER_US
        else:
            width = GAP_BASE_PX + GAP_LOG_PX * math.log10(1.0 + dur)
        px_map.append((t0, t1, px_cursor, px_cursor + width))
        px_cursor += width

    def transform(t0, t1):
        for (a, b, pa, pb) in px_map:
            if a <= t0 <= b + 1e-9:
                if b > a:
                    x0 = pa + (pb - pa) * (t0 - a) / (b - a)
                else:
                    x0 = pa
                x1 = pa + (pb - pa) * (min(t1, b) - a) / (b - a) if b > a else pb
                if t1 > b:
                    for (a2, b2, pa2, pb2) in px_map:
                        if a2 <= t1 <= b2 + 1e-9:
                            x1 = pa2 + (pb2 - pa2) * (t1 - a2) / (b2 - a2) if b2 > a2 else pb2
                            break
                return x0, max(x1 - x0, 0.5)
        return px_cursor, 0.5

    return transform, px_cursor, px_map, segments


def render_svg(events, total_time, seq_name, src_name):
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

    transform, total_px, px_map, segments = build_time_transform(events, total_time)
    width = LEFT_MARGIN + total_px + RIGHT_MARGIN
    height = TOP_MARGIN + len(active_lanes) * (LANE_HEIGHT + LANE_GAP) + BOTTOM_MARGIN

    svg = []
    svg.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {0:.0f} {1:.0f}" '
                'width="100%" style="min-width:{0:.0f}px" font-family="ui-monospace,Menlo,monospace">'
                .format(width, height))
    svg.append('<rect x="0" y="0" width="{0:.0f}" height="{1:.0f}" fill="var(--bg,#0e1117)"/>'
                .format(width, height))

    svg.append('<text x="{0}" y="28" font-size="16" font-weight="600" fill="var(--fg,#e6e6e6)">{1}</text>'
                .format(LEFT_MARGIN, seq_name))
    svg.append('<text x="{0}" y="46" font-size="11" fill="var(--muted,#9aa4b2)">{1} -- total {2} '
                '(one scan; idle gaps &gt;{3:.0f}us log-compressed, hover for exact values)</text>'
                .format(LEFT_MARGIN, src_name, _fmt_us(total_time), GAP_COMPRESS_THRESHOLD))

    # gap break markers + axis ticks along the top
    for (a, b, pa, pb) in px_map:
        dur = b - a
        if dur > GAP_COMPRESS_THRESHOLD:
            cx = LEFT_MARGIN + (pa + pb) / 2.0
            svg.append('<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="var(--grid,#2a2f3a)" '
                        'stroke-width="1" stroke-dasharray="3,3"/>'
                        .format(cx, TOP_MARGIN - 14, height - BOTTOM_MARGIN))
            svg.append('<text x="{0:.1f}" y="{1}" font-size="10" fill="var(--muted,#9aa4b2)" '
                        'text-anchor="middle">{2}</text>'
                        .format(cx, TOP_MARGIN - 18, _fmt_us(dur)))

    for i, lane in enumerate(active_lanes):
        y0 = TOP_MARGIN + i * (LANE_HEIGHT + LANE_GAP)
        y_mid = y0 + LANE_HEIGHT / 2.0
        label = LANE_LABELS.get(lane, lane)
        svg.append('<text x="{0}" y="{1:.1f}" font-size="12" fill="var(--fg,#e6e6e6)" '
                    'text-anchor="end" dominant-baseline="middle">{2}</text>'
                    .format(LEFT_MARGIN - 12, y_mid, label))
        svg.append('<line x1="{0}" y1="{1:.1f}" x2="{2:.0f}" y2="{1:.1f}" stroke="var(--grid,#2a2f3a)" '
                    'stroke-width="1"/>'.format(LEFT_MARGIN, y_mid, width - RIGHT_MARGIN))

        for e in events:
            if e['lane'] != lane:
                continue
            x, w = transform(e['t_start'], e['t_end'])
            x += LEFT_MARGIN
            color = LANE_COLORS.get(e['kind'], '#888')
            title_bits = ["t={0}".format(_fmt_us(e['t_start'])), "dur={0}".format(_fmt_us(e['duration']))]
            for k, v in e['meta'].items():
                title_bits.append("{0}={1}".format(k, v))
            title = " ".join(title_bits)

            if e['kind'] == 'ramp':
                v0 = e['meta'].get('value_from', 0.0)
                v1 = e['meta'].get('value_to', 0.0)
                lane_top = y0 + 6
                lane_bot = y0 + LANE_HEIGHT - 6
                def y_for(v):
                    v = max(-1.0, min(1.0, float(v)))
                    return lane_bot - (v - (-1.0)) / 2.0 * (lane_bot - lane_top)
                y_a, y_b = y_for(v0), y_for(v1)
                svg.append('<polyline points="{0:.1f},{1:.1f} {2:.1f},{3:.1f}" stroke="{4}" '
                            'stroke-width="3" fill="none"><title>{5}</title></polyline>'
                            .format(x, y_a, x + w, y_b, color, title))
            elif e['kind'] == 'marker':
                svg.append('<circle cx="{0:.1f}" cy="{1:.1f}" r="4" fill="{2}"><title>{3}</title></circle>'
                            .format(x, y_mid, color, title))
            elif e['kind'] == 'config':
                svg.append('<rect x="{0:.1f}" y="{1:.1f}" width="{2:.1f}" height="4" fill="{3}" '
                            'opacity="0.6"><title>{4}</title></rect>'
                            .format(x, y0 + LANE_HEIGHT - 4, max(w, 1.5), color, title))
            else:
                ry = y0 + 6 if e['kind'] == 'pulse' else y0 + 6
                rh = LANE_HEIGHT - 12
                rx_style = 'fill="{0}"'.format(color) if e['kind'] == 'pulse' else \
                           'fill="none" stroke="{0}" stroke-width="2" stroke-dasharray="4,2"'.format(color)
                svg.append('<rect x="{0:.1f}" y="{1:.1f}" width="{2:.1f}" height="{3:.1f}" rx="3" {4}>'
                            '<title>{5}</title></rect>'
                            .format(x, ry, max(w, 1.5), rh, rx_style, title))
                if w > 24:
                    label_txt = e['meta'].get('label')
                    if not label_txt:
                        if e['kind'] == 'pulse':
                            label_txt = "ph={0:g}".format(e['meta'].get('phase', 0))
                        elif e['kind'] == 'acqu':
                            label_txt = "ACQU"
                    if label_txt:
                        svg.append('<text x="{0:.1f}" y="{1:.1f}" font-size="9" '
                                    'fill="{2}" text-anchor="middle" dominant-baseline="middle" '
                                    'pointer-events="none">{3}</text>'
                                    .format(x + w / 2.0, y_mid, '#0e1117' if e['kind'] == 'pulse' else color, label_txt))

    svg.append('</svg>')
    return "\n".join(svg)


def render_html(svg, seq_name, pp_path, param_summary):
    return """<!doctype html>
<title>{title}</title>
<style>
  :root {{ --bg:#0e1117; --fg:#e6e6e6; --muted:#9aa4b2; --grid:#2a2f3a; }}
  @media (prefers-color-scheme: light) {{
    :root:not([data-theme="dark"]) {{ --bg:#ffffff; --fg:#1a1a1a; --muted:#6b7280; --grid:#e2e5eb; }}
  }}
  :root[data-theme="light"] {{ --bg:#ffffff; --fg:#1a1a1a; --muted:#6b7280; --grid:#e2e5eb; }}
  body {{ margin:0; background:var(--bg); color:var(--fg); font-family: -apple-system, sans-serif; }}
  .wrap {{ padding: 20px; overflow-x: auto; }}
  .params {{ padding: 0 20px 20px; font-family: ui-monospace, Menlo, monospace; font-size: 12px; color: var(--muted); white-space: pre-wrap; }}
  h2 {{ font-size: 13px; color: var(--muted); font-weight: 500; margin: 0 0 6px; }}
</style>
<div class="wrap">
{svg}
</div>
<div class="params">
<h2>Parameters used for this trace ({pp_path})</h2>{params}
</div>
""".format(title=seq_name, svg=svg, pp_path=pp_path, params=param_summary)


def summarize_params(pp_path):
    mod = runpy.run_path(pp_path, run_name="__paramdump__")
    P = mod['Parameters']
    lines = []
    for name in sorted(dir(P)):
        if name.startswith('_'):
            continue
        val = getattr(P, name)
        if callable(val):
            continue
        lines.append("{0} = {1}".format(name, val))
    return "\n".join(lines)


def parse_overrides(pairs):
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
    args = ap.parse_args()

    if not os.path.isfile(args.pp_file):
        raise SystemExit("No such file: {0}".format(args.pp_file))

    overrides = parse_overrides(args.overrides)
    events, total_time, seq_name = trace_pp_file(args.pp_file, overrides)
    svg = render_svg(events, total_time, seq_name, os.path.basename(args.pp_file))

    # Re-run once more (fresh mock, undisturbed) just to dump the final
    # Parameter values actually used, for the caption under the diagram.
    clock2 = Clock()
    install_mock_firebird(clock2, [])
    mod = runpy.run_path(args.pp_file)
    P = mod['Parameters']
    P.NumScans, P.DS = 1, 0
    if overrides:
        for key, value in overrides.items():
            current = getattr(P, key)
            caster = type(current) if not isinstance(current, bool) else str
            try:
                setattr(P, key, caster(value))
            except (TypeError, ValueError):
                setattr(P, key, value)
    lines = []
    cls = type(P)
    for name in sorted(dir(cls)):
        if name.startswith('_'):
            continue
        val = getattr(P, name)
        if callable(val):
            continue
        lines.append("{0} = {1}".format(name, val))
    param_summary = "\n".join(lines)

    html = render_html(svg, seq_name, args.pp_file, param_summary)

    out_path = args.output or (os.path.splitext(args.pp_file)[0] + ".timing.html")
    with open(out_path, "w") as f:
        f.write(html)
    print("Wrote {0} ({1} events, total {2})".format(out_path, len(events), _fmt_us(total_time)))


if __name__ == '__main__':
    main()
