#!/usr/bin/env python3
"""
pp_tikz.py -- publication-quality TikZ pulse-sequence diagram generator.

Reuses pp_visualizer.py's tracer (same instrumented mock, same
with-parallel:/with-sequential: clock semantics, same gradient-trapezoid
merging) to emit a standalone, compilable TikZ figure, styled to match
the user-supplied reference (Tex_figures/tex_files/ste_PFG.tex): thick
black rectangles for RF pulses, filled cyan trapezoids for gradient
pulses with a delta<->delta arrow over the flat top, a <-> arrow over
the most notable inter-pulse gap, and a decaying-oscillation
"Detection" squiggle after the final acquisition.

Callable API:
    from pp_tikz import visualize_tikz
    result = visualize_tikz("sequences/diffusion/LWG_PGSTE_H.py",
                             overrides={"DELTA": 20000})
    # result = {'tex': str, 'output_path': ...}

CLI:
    python3 tools/pp_tikz.py sequences/diffusion/LWG_PGSTE_H.py
    python3 tools/pp_tikz.py sequences/diffusion/LWG_PGSTE_H.py -o out.tex

Output is a standalone document (\\documentclass[margin=4pt]{standalone})
using the same tikz libraries as the reference file -- compile directly
with pdflatex/lualatex, no extra setup needed beyond a standard TeX Live
with the tikz/circuitikz/tikz-3dplot/chemformula packages the reference
already depends on (chemformula/tikz-3dplot/circuitikz aren't actually
USED by anything this generator draws, so they're OMITTED here -- only
tikz + its listed libraries, which is everything this file's own output
needs).

Known v1 simplifications (schematic, not literally to scale -- matching
the reference document's own convention, where pulse width is a fixed
schematic \\pi2, not a literal microsecond-accurate draw):
  - RF pulses are drawn with a MINIMUM visible width regardless of their
    true (often <20us) duration -- true-to-scale would make most hard
    pulses invisible hairlines at any reasonable page width. Gradient
    pulses/plateaus ARE drawn proportionally to their real duration
    (relative to each other and to the labelled gap) -- that IS the
    quantity being illustrated.
  - Pulses are labelled with a flip angle ONLY when their traced width
    matches a Parameter whose name contains '90' or '180' (e.g. P90,
    P180) to within 0.05us -- otherwise labelled generically (P1, P2,
    ...). This is a best-effort match, not guaranteed for every file.
  - Only ONE inter-event gap gets an explicit <-> arrow (the largest
    gap that sits between two gradient trapezoids if there are any,
    otherwise the largest gap between two RF pulses) -- labelled with
    its numeric duration, not the originating Parameter's name (that
    mapping isn't recoverable from the trace alone).
"""

import argparse
import math
import os

from pp_visualizer import (
    Clock, install_mock_firebird, trace_pp_file, merge_gradient_trapezoids,
    dump_params, GAP_COMPRESS_THRESHOLD,
)


# =============================================================================
# Layout constants (TikZ units, ~= cm at scale=1)
# =============================================================================

EVENT_UNITS_PER_US = 0.0006     # true-to-scale width for gradient trapezoids
GAP_BASE_UNITS = 0.3
GAP_LOG_UNITS = 0.35
MIN_PULSE_WIDTH = 0.4            # schematic minimum, regardless of true duration
MAX_EVENT_WIDTH = 2.0            # schematic cap for pulses/ACQU (e.g. a 30ms
                                  # acquisition window isn't meant to be shown
                                  # true-to-scale any more than a 10us hard
                                  # pulse is -- both are schematic in a
                                  # pulse-sequence figure; gradients are the
                                  # exception, see 'grad' branch below)
MIN_GRAD_RAMP = 0.12
LANE_SPACING = 3.0
PULSE_HEIGHT = 1.6
GRAD_HEIGHT = 1.0
RF_LABEL_HINTS = ['90', '180', '270']


def _fmt_us(v):
    v = float(v)
    if v >= 1.0e6:
        return "{0:.3g} s".format(v / 1.0e6)
    if v >= 1000.0:
        return "{0:.3g} ms".format(v / 1000.0)
    return "{0:.3g} \\textmu s".format(v)


class TikzAxis(object):
    """Same piecewise true-scale/log-compressed time mapping as
    pp_visualizer.TimeAxis, but calibrated to TikZ units and using a
    MINIMUM width floor for non-gradient events (pulses/acqu), since a
    true-to-scale ~10us hard pulse would be an invisible hairline in a
    printed figure -- gradient trapezoids keep their real proportional
    width, which is the point of drawing them at all."""

    def __init__(self, events, total_time):
        occupied = sorted((e['t_start'], e['t_end']) for e in events if e['duration'] > 0)
        merged = []
        for s, en in occupied:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], en))
            else:
                merged.append((s, en))

        # Which merged occupied spans are "gradient-only" (true-scale, no
        # min-width floor) vs "RF/other" (schematic min-width floor
        # applied) -- a merged span counts as gradient-only if it
        # exactly matches one gradient trapezoid's own [t_start,t_end),
        # i.e. nothing on another lane overlaps it in time. Any span
        # that doesn't exactly match (e.g. a 'with parallel:' block
        # where RF/ACQU timing genuinely overlaps a gradient) falls back
        # to the 'event' (min-width-floor) treatment -- safe, just
        # slightly less proportionally-accurate for that rarer case.
        grad_spans = set()
        for e in events:
            if e['lane'].startswith('GRAD') and e['duration'] > 0:
                grad_spans.add((e['t_start'], e['t_end']))

        segments = []
        cursor = 0.0
        for s, en in merged:
            if s > cursor:
                segments.append((cursor, s, 'gap'))
            segments.append((s, en, 'grad' if (s, en) in grad_spans else 'event'))
            cursor = en
        if cursor < total_time:
            segments.append((cursor, total_time, 'gap'))

        self.px_map = []
        cursor_u = 0.0
        for t0, t1, kind in segments:
            dur = t1 - t0
            if dur <= 0:
                continue
            if kind == 'grad':
                width = max(dur * EVENT_UNITS_PER_US, MIN_GRAD_RAMP)
            elif kind == 'event':
                width = min(max(dur * EVENT_UNITS_PER_US, MIN_PULSE_WIDTH), MAX_EVENT_WIDTH)
            elif dur <= GAP_COMPRESS_THRESHOLD:
                width = dur * EVENT_UNITS_PER_US
            else:
                width = GAP_BASE_UNITS + GAP_LOG_UNITS * math.log10(1.0 + dur)
            self.px_map.append((t0, t1, cursor_u, cursor_u + width))
            cursor_u += width
        self.total_units = cursor_u

    def px(self, t):
        for (a, b, pa, pb) in self.px_map:
            if a <= t <= b + 1e-9:
                return pa if b <= a else pa + (pb - pa) * (t - a) / (b - a)
        return self.total_units

    def span(self, t0, t1):
        x0, x1 = self.px(t0), self.px(t1)
        return x0, max(x1 - x0, 0.05)

    def gaps(self):
        return [(a, b, pa, pb) for (a, b, pa, pb) in self.px_map if (b - a) > GAP_COMPRESS_THRESHOLD]

    def headline_gaps(self):
        """Only the gap(s) worth an explicit break-marker in a printed
        figure -- i.e. within 50% of the single largest gap's duration.
        In practice this is just RD (500us-2s), which dominates every
        other structural delay in this repo by 1-3 orders of magnitude;
        marking every >60us gap (as the interactive HTML view does) would
        bury a static figure in clutter."""
        all_gaps = self.gaps()
        if not all_gaps:
            return []
        max_dur = max(b - a for (a, b, pa, pb) in all_gaps)
        return [g for g in all_gaps if (g[1] - g[0]) >= max_dur * 0.5]


def _guess_pulse_label(width, params, seq_counters, lane):
    # Tolerance 0.55us, not 0.05: shaped_pulse() rounds its Parameter
    # duration to whole microseconds internally (n_steps = int(round(
    # duration))) before synthesising the envelope, so a traced width can
    # be up to 0.5us off its originating Parameter's exact value.
    for name, val in params.items():
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            continue
        if abs(float(val) - width) < 0.55:
            low = name.lower()
            for hint in RF_LABEL_HINTS:
                if hint in low:
                    return "{0}\\textdegree".format(hint)
    seq_counters[lane] = seq_counters.get(lane, 0) + 1
    return "P{0}".format(seq_counters[lane])


def _guess_shape_name(label, params, used_shape_keys):
    """Best-effort match of a shaped pulse's flip-angle label (e.g.
    '90\\textdegree', from _guess_pulse_label()) to the *Shape Parameter
    that actually generated it -- this repo's selective/imaging-selective
    family consistently names these 'ExcitationShape' (90) and
    'RefocusShape' (180/270), see e.g. LWG_Selective-Echo_H.py. Returns
    None (no label suffix) rather than guessing wrong when the hint
    doesn't match anything -- e.g. WetShape, which isn't tied to a flip
    angle this way."""
    candidates = {k: v for k, v in params.items()
                  if 'shape' in k.lower() and k not in used_shape_keys
                  and isinstance(v, str) and v}
    low_label = label.lower()
    for key in sorted(candidates):
        klow = key.lower()
        if ('90' in low_label and 'excit' in klow) or \
           (('180' in low_label or '270' in low_label) and 'refocus' in klow):
            used_shape_keys.add(key)
            return candidates[key]
    return None


def _downsample_profile(profile, max_points=48):
    n = len(profile)
    if n <= max_points:
        return profile
    step = (n - 1) / float(max_points - 1)
    idx = sorted(set(int(round(i * step)) for i in range(max_points)))
    return [profile[i] for i in idx]


def _lane_v_to_y(v, y_base):
    """Gradient value (-1..1) -> y-offset from the lane baseline, using
    GRAD_HEIGHT as the full-scale rise for |v|=1 -- matches the
    reference figure's own trapezoid convention exactly ((0,b)--(.,b+1))."""
    return y_base + max(-1.0, min(1.0, float(v))) * GRAD_HEIGHT


LANE_ORDER = ['HF', 'X', 'GRAD1', 'GRAD2', 'GRAD3']
LANE_LABEL = {'HF': 'H/F', 'X': 'X', 'GRAD1': 'G$_x$', 'GRAD2': 'G$_y$', 'GRAD3': 'G$_z$'}


def render_tikz(events, total_time, seq_name, params, show_values=False):
    events = merge_gradient_trapezoids(events)

    active_lanes = [lane for lane in LANE_ORDER if any(e['lane'] == lane for e in events)]
    for e in events:
        if e['lane'] not in active_lanes and e['lane'] != 'TRIGGER':
            active_lanes.append(e['lane'])
    if not active_lanes:
        active_lanes = ['HF']

    axis = TikzAxis(events, total_time)

    # y-baseline for each lane, RF lanes on top, gradient lanes below,
    # top-to-bottom in LANE_ORDER (matching NMR pulse-sequence convention).
    y_of = {}
    for i, lane in enumerate(active_lanes):
        y_of[lane] = -i * LANE_SPACING

    lines = []
    lines.append("\\documentclass[margin=4pt]{standalone}")
    lines.append("\\usepackage[euler]{textgreek}")
    lines.append("\\usepackage{sansmath}")
    lines.append("\\sansmath")
    lines.append("\\usepackage{tikz}")
    lines.append("\\tikzstyle{every picture}+=[font=\\sffamily\\bfseries]")
    lines.append("\\usetikzlibrary{arrows,calc,positioning}")
    lines.append("\\tikzset{>=latex}")
    lines.append("")
    lines.append("% Auto-generated by tools/pp_tikz.py from {0}".format(_esc(seq_name)))
    lines.append("% -- schematic, not literally to scale (see module docstring)")
    lines.append("")
    lines.append("\\begin{document}")
    lines.append("\\begin{tikzpicture}[scale=0.7, every node/.style={scale=0.9}]")
    lines.append("")

    total_x = axis.total_units
    seq_counters = {}
    used_shape_keys = set()

    # ---- Lane baselines + labels ----
    for lane in active_lanes:
        y = y_of[lane]
        lines.append("\\draw[ultra thick] (0,{0:.3f}) node[above]{{\\Large {1}}} "
                      "(0,{0:.3f})--({2:.3f},{0:.3f});".format(y, LANE_LABEL.get(lane, lane), total_x + 0.6))
    lines.append("")

    # ---- Idle-gap break markers -- RD only, see headline_gaps() ----
    # Labelled with the actual duration only when show_values=True; by
    # default this shows the generic "Recycle Delay" name instead, since
    # this gap is (per headline_gaps()'s own docstring) always RD in
    # practice, and a schematic figure usually wants the delay NAMED, not
    # its specific numeric value.
    for (a, b, pa, pb) in axis.headline_gaps():
        cx = (pa + pb) / 2.0
        top_y = max(y_of.values()) + PULSE_HEIGHT + 1.0
        bot_y = min(y_of.values()) - GRAD_HEIGHT - 0.6
        gap_label = _fmt_us(b - a) if show_values else "Recycle Delay"
        lines.append("\\draw[densely dotted] ({0:.3f},{1:.3f})--({0:.3f},{2:.3f});"
                      .format(cx, top_y, bot_y))
        lines.append("\\node[above] at ({0:.3f},{1:.3f}) {{\\small {2}}};"
                      .format(cx, top_y, gap_label))

    # ---- Per-lane events ----
    last_acqu = None  # (lane, x_end)
    for lane in active_lanes:
        y = y_of[lane]
        for e in events:
            if e['lane'] != lane:
                continue
            if e['kind'] == 'pulse':
                x0, w = axis.span(e['t_start'], e['t_end'])
                label = _guess_pulse_label(e['duration'], params, seq_counters, lane)
                profile = e['meta'].get('shape_profile')
                if profile and len(profile) >= 4:
                    shape_name = _guess_shape_name(label, params, used_shape_keys)
                    full_label = "{0} {1}".format(label, _esc(shape_name)) if shape_name else label
                    pts = _downsample_profile(profile)
                    peak = max(abs(v) for (_, v) in pts) or 1.0
                    curve_pts = ["({0:.3f},{1:.3f})".format(axis.px(e['t_start'] + t),
                                                              y + PULSE_HEIGHT * (v / peak))
                                 for (t, v) in pts]
                    path = ("({0:.3f},{1:.3f})--".format(x0, y) + "--".join(curve_pts) +
                            "--({0:.3f},{1:.3f})--cycle".format(x0 + w, y))
                    lines.append("\\draw[thick, fill=black!30] {0};".format(path))
                    lines.append("\\node[above] at ({0:.3f},{1:.3f}) {{{2}}};"
                                  .format(x0 + w / 2.0, y + PULSE_HEIGHT + 0.15, full_label))
                else:
                    lines.append("\\draw[fill=black] ({0:.3f},{1:.3f}) node[above]{{{2}}} "
                                  "({3:.3f},{4:.3f})rectangle({5:.3f},{6:.3f});"
                                  .format(x0 + w / 2.0, y + PULSE_HEIGHT + 0.15, label,
                                          x0, y, x0 + w, y + PULSE_HEIGHT))
            elif e['kind'] == 'acqu':
                x0, w = axis.span(e['t_start'], e['t_end'])
                last_acqu = (lane, x0 + w)
            elif e['kind'] == 'trapezoid':
                verts = e['meta']['vertices']
                pts = ["({0:.3f},{1:.3f})".format(axis.px(t), _lane_v_to_y(v, y)) for (t, v) in verts]
                path = "--".join(pts)
                lines.append("\\draw[ultra thick, fill=cyan] {0};".format(path))
                for (pt0, pt1, pv) in e['meta']['plateaus']:
                    px0, pw = axis.span(pt0, pt1)
                    if pw < 0.15:
                        continue
                    ay = _lane_v_to_y(pv, y) + (0.28 if pv >= 0 else -0.28)
                    label_side = "above" if pv >= 0 else "below"
                    lines.append("\\draw[<->] ({0:.3f},{1:.3f})--({2:.3f},{1:.3f}) "
                                  "node[{3}=0.03,midway]{{\\large \\textdelta}};"
                                  .format(px0, ay, px0 + pw, label_side))
        lines.append("")

    # ---- One headline <-> arrow for the most notable inter-event gap ----
    big_gap = _pick_headline_gap(events, axis)
    if big_gap:
        (t0, t1, lane_hint) = big_gap
        x0, x1 = axis.px(t0), axis.px(t1)
        y = min(y_of.values()) - GRAD_HEIGHT - 1.1 if lane_hint == 'below' else \
            max(y_of.values()) + PULSE_HEIGHT + 1.6
        delta_label = "\\textDelta\\ ({0})".format(_fmt_us(t1 - t0)) if show_values else "\\textDelta"
        lines.append("\\draw[<->, ultra thick] ({0:.3f},{1:.3f})--({2:.3f},{1:.3f}) "
                      "node[below=0.05,midway]{{\\large {3}}};"
                      .format(x0, y, x1, delta_label))
        lines.append("")

    # ---- Decaying "Detection" squiggle after the final acquisition ----
    # (same functional form as the reference figure's own FID squiggle:
    # a decaying two-frequency cosine/sine sum, e^(-x/4) envelope).
    if last_acqu:
        lane, x_start = last_acqu
        y = y_of[lane]
        x_expr = "(\\x/pi)+{0:.3f}".format(x_start + 0.3)
        y_expr = ("(1*cos(15+8*\\x r)+0.5*sin(75+10*\\x r))*e^(-\\x/4)+{0:.3f}".format(y))
        lines.append("% Detection")
        lines.append("\\draw[ultra thick] plot[domain=0:3.6*pi, samples=800] "
                      "({{{0}}},{{{1}}});".format(x_expr, y_expr))
        lines.append("")

    lines.append("\\end{tikzpicture}")
    lines.append("")
    lines.append("\\end{document}")
    return "\n".join(lines)


def _pick_headline_gap(events, axis):
    """Return (t0, t1, 'below'|'above') for the single gap worth an
    explicit <-> arrow: the largest gap between two gradient trapezoids
    if any exist (that's the diffusion-time-style interval a reference
    figure like ste_PFG.tex highlights), else the largest gap between
    two RF pulses."""
    grad_events = sorted([e for e in events if e['kind'] == 'trapezoid'], key=lambda e: e['t_start'])
    if len(grad_events) >= 2:
        best = None
        for a, b in zip(grad_events, grad_events[1:]):
            gap = (a['t_end'], b['t_start'])
            if gap[1] > gap[0] and (best is None or gap[1] - gap[0] > best[1] - best[0]):
                best = gap
        if best:
            return (best[0], best[1], 'below')

    pulses = sorted([e for e in events if e['kind'] == 'pulse'], key=lambda e: e['t_start'])
    best = None
    for a, b in zip(pulses, pulses[1:]):
        gap = (a['t_end'], b['t_start'])
        if best is None or gap[1] - gap[0] > best[1] - best[0]:
            best = gap
    if best and best[1] > best[0]:
        return (best[0], best[1], 'above')
    return None


def _esc(s):
    return str(s).replace('_', '\\_').replace('&', '\\&').replace('%', '\\%')


def visualize_tikz(pp_path, overrides=None, output=None, show_values=False):
    """Trace pp_path and render a standalone TikZ (.tex) pulse-sequence
    diagram, styled to match Tex_figures/tex_files/ste_PFG.tex.

    overrides: dict of {ParameterName: value} to apply before tracing.
    output: path to write the .tex to. Defaults to
        '<pp_path without .py>.pp.tex'. Pass output=False to skip
        writing to disk.
    show_values: if True, label the recycle-delay break marker with its
        actual duration (e.g. "1.92 s") instead of the default generic
        "Recycle Delay" -- off by default (numbers are opt-in, not the
        default, for a schematic figure). Gradient delta/Delta labels
        always show their actual duration regardless of this flag.

    Returns {'tex': str, 'output_path': str or None, 'sequence_name': str}.
    """
    events, total_time, seq_name = trace_pp_file(pp_path, overrides)
    params = dump_params(pp_path, overrides)
    tex = render_tikz(events, total_time, seq_name, params, show_values=show_values)

    out_path = None
    if output is not False:
        out_path = output or (os.path.splitext(pp_path)[0] + ".pp.tex")
        with open(out_path, "w") as f:
            f.write(tex)

    return {'tex': tex, 'output_path': out_path, 'sequence_name': seq_name}


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
    ap.add_argument('-o', '--output', help="Output .tex path (default: alongside input, .pp.tex)")
    ap.add_argument('--set', action='append', dest='overrides', metavar='KEY=VALUE',
                     help="Override a Parameter before tracing (repeatable)")
    ap.add_argument('--show-values', action='store_true',
                     help="Label the recycle-delay break marker with its actual "
                          "duration instead of the default generic 'Recycle Delay'")
    args = ap.parse_args()

    if not os.path.isfile(args.pp_file):
        raise SystemExit("No such file: {0}".format(args.pp_file))

    overrides = _parse_overrides(args.overrides)
    result = visualize_tikz(args.pp_file, overrides=overrides, output=args.output,
                             show_values=args.show_values)
    print("Wrote {0}".format(result['output_path']))


if __name__ == '__main__':
    main()
