#!/usr/bin/env python3
"""
spinflow_bvalue_diff.py -- write a b-value table (.diff file) into a
SpinFlow experiment's own parameter folder.

WHY THIS EXISTS (and why the pulse programme can't do it itself)
---------------------------------------------------------------
A pulse programme's run() does NOT execute on the Windows PC. SpinFlow
copies the sequence to the spectrometer console and runs it there --
visible in the instrument's own error log as

    File "/tmp/tmpt65cdn.py", line 600, in run
    File "/home/fire/PythonMonitor/firebird/parallel/utils.py", ...

while the parameter files live on the host at

    C:\\Users\\Public\\Documents\\SpinFlow\\Parameters\\<serial>\\User\\<nuc>\\<Sequence>.par

A file written from run() would land on the console's filesystem, not next
to your parameters. So the sequences log their b-values and stamp them into
the saved JCAMP metadata, and THIS script -- run on the SpinFlow PC --
writes the .diff file where you actually want it.

SINGLE SOURCE OF TRUTH
----------------------
This script does NOT re-implement the b-value maths. It imports
bvalue_segments(), integrate_bvalue(), MAXGRAD_TABLE and
NUCLEUS_GAMMA_MHZ_PER_T straight out of the pulse programme .py itself
(under a stub `firebird`, so no instrument libraries are needed), then
overlays the .par's saved values on top of the sequence's own Parameter
defaults. Whatever the sequence would have computed is what you get here --
there is nothing to keep in sync.

USAGE
-----
    # point it at a saved .par; finds the matching sequence automatically
    python spinflow_bvalue_diff.py "C:\\Users\\Public\\Documents\\SpinFlow\\Parameters\\82863\\User\\1H\\LWG_PFGSTE_H.par"

    # or search the whole Parameters tree for a sequence by name
    python spinflow_bvalue_diff.py --find LWG_DPFGSTE_H

    # override the gradient list without editing the .par
    python spinflow_bvalue_diff.py LWG_PFGSTE_H.par --grad-list 0.05,0.3,0.6,0.95

    # no .par yet -- just use the sequence's own defaults
    python spinflow_bvalue_diff.py --sequence ../sequences/diffusion/LWG_PFGSTE_H.py --defaults

Writes <Sequence>.diff beside the .par (or beside the sequence, with
--defaults). Add --difflist to also write <Sequence>.difflist, a bare
one-value-per-line list of the relative gradients.
"""

import argparse
import datetime
import importlib.util
import os
import sys
import types


# ---------------------------------------------------------------------------
# Stub firebird, so a pulse programme can be imported off-instrument
# ---------------------------------------------------------------------------

_FIREBIRD_NAMES = [
    "Delay", "Filter", "ChooseFilter", "PhasesManager", "PhaseListContainer",
    "GradientMatrix", "Gradient1", "Gradient2", "Gradient3",
    "Gradient1SlewRate", "Gradient2SlewRate", "Gradient3SlewRate",
    "TX0", "TX1", "sequential", "sequential_main", "parallel", "start",
    "WriteToHardware", "get_single_scan_execution_time",
    "Get_Compilation_Time", "ExternalTrigger1", "ExternalTrigger2",
    "ExternalTrigger3",
]
for _ch in ("1", "2"):
    for _pre, _sufs in (
            ("Channel", ["SetFrequency", "RestartSynth", "SetBasePhase"]),
            ("Transmit", ["SelectPort", "LPEnable", "SetScale", "BlankingOn",
                          "BlankingOff", ""]),
            ("Receiver", ["Preamp", "Filter", "FilterFlush", "Phase", ""])):
        for _s in _sufs:
            _FIREBIRD_NAMES.append(_pre + _ch + _s)


class _Stub(object):
    def __call__(self, *a, **k): return None
    def __enter__(self): return self
    def __exit__(self, *a): return False
    def __getattr__(self, n): return _Stub()


class ParameterTypes(object):
    String = "String"
    Double = "Double"
    Int32 = "Int32"


class _Param(object):
    __slots__ = ("key", "value", "ptype", "desc")

    def __init__(self, key, value, ptype, desc):
        self.key = key
        self.value = value
        self.ptype = ptype
        self.desc = desc


def _Parameter(key, default, ptype, desc, validator=None, **kw):
    return _Param(key, default, ptype, desc)


def _ParameterBlock(cls):
    """Collapse each Parameter attribute to its default value, keeping a
    _params map (attribute name -> _Param) so .par keys can be matched
    against SpinFlow short codes as well as attribute names."""
    params = {}
    for name in list(vars(cls)):
        value = getattr(cls, name)
        if isinstance(value, _Param):
            params[name] = value
            setattr(cls, name, value.value)
    cls._params = params
    return cls


def load_sequence(path):
    """Import a pulse programme's module-level code under the stub."""
    fb = types.ModuleType("firebird")
    fb.Parameter = _Parameter
    fb.ParameterBlock = _ParameterBlock
    fb.ParameterTypes = ParameterTypes
    for name in _FIREBIRD_NAMES:
        setattr(fb, name, _Stub())
    fb.__all__ = ["Parameter", "ParameterBlock", "ParameterTypes"] + _FIREBIRD_NAMES
    apps = types.ModuleType("firebird.applications")
    apps.__all__ = []
    fb.applications = apps
    sys.modules["firebird"] = fb
    sys.modules["firebird.applications"] = apps

    spec = importlib.util.spec_from_file_location(
        "pp_" + os.path.splitext(os.path.basename(path))[0], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for needed in ("Parameters", "bvalue_segments", "integrate_bvalue",
                   "MAXGRAD_TABLE", "NUCLEUS_GAMMA_MHZ_PER_T"):
        if not hasattr(mod, needed):
            raise SystemExit(
                "{0} has no {1}() -- this script only works with the house-style "
                "diffusion sequences (LWG_PFGSTE_*/LWG_DPFGSTE_*), not the vendor "
                "originals.".format(os.path.basename(path), needed))
    return mod


# ---------------------------------------------------------------------------
# .par parsing -- format auto-detected (XML first, then key/value lines)
# ---------------------------------------------------------------------------

def parse_par(path):
    """Return {key: string_value} from a SpinFlow .par file.

    The exact .par serialisation has not been confirmed here, so both plausible
    shapes are tried and whichever yields more keys wins. Run with --show-par
    to print what was actually recovered before trusting the numbers.
    """
    raw = open(path, "rb").read()
    text = raw.decode("utf-8-sig", errors="replace")

    xml_values = {}
    try:
        import xml.etree.ElementTree as ET
        root = ET.fromstring(text)
        for el in root.iter():
            name = el.get("name") or el.get("Name") or el.get("key") or el.get("Key")
            value = el.get("value")
            if value is None:
                value = el.get("Value")
            if name is None:
                child_name = el.find("Name")
                child_value = el.find("Value")
                if child_name is not None and child_value is not None:
                    name = (child_name.text or "").strip()
                    value = (child_value.text or "").strip()
            if name and value is not None:
                xml_values[str(name).strip()] = str(value).strip()
            elif name and el.text and el.text.strip():
                xml_values[str(name).strip()] = el.text.strip()
    except Exception:
        pass

    flat_values = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", ";", "<", "--")):
            continue
        for sep in ("=", "\t", ":"):
            if sep in line:
                key, _, value = line.partition(sep)
                key = key.strip()
                value = value.strip().strip('"')
                if key and not key.startswith("<"):
                    flat_values[key] = value
                break

    return xml_values if len(xml_values) >= len(flat_values) else flat_values


def overlay(Parameters, par_values, verbose=False):
    """Apply .par values onto a copy of the sequence's Parameters defaults,
    matching by Python attribute name OR by SpinFlow short code, and casting
    to the type of the existing default."""
    class P(object):
        pass
    p = P()
    for name in dir(Parameters):
        if not name.startswith("_"):
            setattr(p, name, getattr(Parameters, name))

    # Match CASE-SENSITIVELY first. This matters: the house rename gives these
    # sequences both a "DELTA" (diffusion time) and a "delta" (gradient
    # plateau) short code, which differ ONLY in case -- a case-insensitive
    # lookup silently writes one over the other. The case-insensitive pass is
    # a fallback, and only for names that are unambiguous when folded.
    declared = getattr(Parameters, "_params", {})
    exact = {}
    for attr, meta in declared.items():
        exact[str(meta.key)] = attr
        exact.setdefault(attr, attr)
    folded_count = {}
    for name in exact:
        folded_count[name.lower()] = folded_count.get(name.lower(), 0) + 1
    folded = dict((name.lower(), attr) for name, attr in exact.items()
                  if folded_count[name.lower()] == 1)

    applied, unknown, ambiguous = [], [], []
    for key, value in par_values.items():
        key = str(key)
        attr = exact.get(key)
        if attr is None:
            attr = folded.get(key.lower())
            if attr is None and key.lower() in folded_count:
                ambiguous.append(key)
                continue
        if attr is None:
            unknown.append(key)
            continue
        current = getattr(p, attr)
        try:
            if isinstance(current, bool):
                cast = value.strip().lower() in ("1", "true", "yes", "on")
            elif isinstance(current, int):
                cast = int(float(value))
            elif isinstance(current, float):
                cast = float(value)
            else:
                cast = value
        except (TypeError, ValueError):
            unknown.append(key + " (uncastable: " + repr(value) + ")")
            continue
        setattr(p, attr, cast)
        applied.append("{0} ({1}) = {2}".format(attr, key, cast))

    if ambiguous:
        print("  WARNING: .par key(s) {0} match a Parameter only when case is "
              "ignored, and more than one Parameter folds to that name (e.g. "
              "DELTA vs delta) -- NOT applied. Check the .par's exact "
              "spelling.".format(", ".join(sorted(ambiguous))))
    if verbose:
        print("  applied from .par ({0}):".format(len(applied)))
        for line in applied:
            print("     ", line)
        if unknown:
            print("  ignored keys ({0}): {1}".format(
                len(unknown), ", ".join(sorted(unknown)[:20])))
    return p


# ---------------------------------------------------------------------------
# b-value table
# ---------------------------------------------------------------------------

def base_frequency(P):
    """SF for an _H sequence, SFX for an _X one -- the vendor's _X files use
    a different short code and attribute name for the observe channel, and
    the house ports keep that (see the _X files' own design notes)."""
    for name in ("FrequencyBase", "FrequencyBaseX", "FrequencyX"):
        if hasattr(P, name):
            return float(getattr(P, name))
    raise SystemExit("No base-frequency Parameter (SF/SFX) found on this sequence.")


def normalise_key(value):
    return str(value).upper().replace("/", "").replace("-", "").replace(" ", "").replace("_", "")


def build_table(mod, P, grad_list=None):
    """Return (rows, meta). rows = [(index, g_rel, G_T_per_m, q, b_s_per_m2)]."""
    axis = str(getattr(P, "Axis", "z")).strip().lower()
    probe = normalise_key(getattr(P, "Probe", ""))

    if axis == "none":
        raise SystemExit("GradAxis is 'none' in this parameter set -- no diffusion "
                         "gradient is applied, so every b-value is 0.")
    if getattr(P, "GradMaxOverride", 0.0) > 0.0:
        maxgrad = float(P.GradMaxOverride)
        maxgrad_source = "GradMax Parameter"
    elif probe not in mod.MAXGRAD_TABLE:
        raise SystemExit("Probe '{0}' is not in MAXGRAD_TABLE (known: {1}). Set the "
                         "Probe Parameter in SpinFlow, or pass --probe or --gmax."
                         .format(getattr(P, "Probe", ""), sorted(mod.MAXGRAD_TABLE)))
    else:
        maxgrad = mod.MAXGRAD_TABLE[probe].get(axis)
        maxgrad_source = "MAXGRAD_TABLE[{0}][{1}]".format(probe, axis)
        if maxgrad is None:
            raise SystemExit("Probe '{0}' has no calibration for the {1} axis. "
                             "Pass --gmax to supply one."
                             .format(getattr(P, "Probe", ""), axis))

    if getattr(P, "GammaOverride", 0.0) > 0.0:
        gamma_over_2pi = float(P.GammaOverride)
        gamma_source = "GammaOverride"
    else:
        nuc = normalise_key(getattr(P, "Nucleus", ""))
        if nuc not in mod.NUCLEUS_GAMMA_MHZ_PER_T:
            raise SystemExit("Nucleus '{0}' is not in the gamma table (known: {1}). "
                             "Set Nucleus, or GammaOverride, or pass --nucleus."
                             .format(getattr(P, "Nucleus", ""),
                                     sorted(mod.NUCLEUS_GAMMA_MHZ_PER_T)))
        gamma_over_2pi = mod.NUCLEUS_GAMMA_MHZ_PER_T[nuc]
        gamma_source = "Nucleus=" + str(getattr(P, "Nucleus", ""))
    gamma = 2.0e6 * 3.141592653589793 * gamma_over_2pi

    fp = {"x": P.XGradNorm, "y": P.YGradNorm, "z": P.ZGradNorm}[axis]

    if grad_list is None:
        grad_list = [float(x) for x in str(P.GradList).split(",") if x.strip()]

    rows = []
    for i, rel in enumerate(grad_list):
        G = abs(rel) * fp * maxgrad * mod.GAUSS_PER_CM_TO_T_PER_M
        b = mod.integrate_bvalue(mod.bvalue_segments(P, G), gamma)
        segs = mod.bvalue_segments(P, G)
        q = gamma * G * ((P.delta + P.RampTime) * 1.0e-6)
        rows.append((i + 1, rel, G, q, b))

    implied_B0 = base_frequency(P) / gamma_over_2pi
    # b = gamma^2 * G_abs^2 * S0 -- see bvalue_shape_factor() in the sequence.
    S0 = mod.bvalue_shape_factor(P)
    b_per_G2 = gamma * gamma * S0
    b_per_Grel2 = b_per_G2 * (fp * maxgrad * mod.GAUSS_PER_CM_TO_T_PER_M) ** 2
    meta = {
        "axis": axis, "probe": getattr(P, "Probe", ""), "maxgrad_G_per_cm": maxgrad,
        "maxgrad_source": maxgrad_source,
        "fp_scaler": fp, "gamma_over_2pi_MHz_per_T": gamma_over_2pi,
        "gamma_source": gamma_source, "gamma_rad_per_s_per_T": gamma,
        "implied_B0_T": implied_B0,
        "shape_factor_s3": S0, "b_per_G2": b_per_G2, "b_per_Grel2": b_per_Grel2,
        "total_gradient_time_us": len(segs) and sum(s[0] for s in segs),
    }
    return rows, meta


def write_diff(path, sequence, rows, meta, P, source_note):
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = []
    lines.append("# b-value table for {0}".format(sequence))
    lines.append("# generated {0} by tools/spinflow_bvalue_diff.py".format(stamp))
    lines.append("# {0}".format(source_note))
    lines.append("#")
    lines.append("# Sequence timing")
    lines.append("#   DELTA (diffusion time)      = {0} us".format(P.DELTA))
    lines.append("#   delta (gradient plateau)    = {0} us".format(P.delta))
    lines.append("#   RampTime                    = {0} us".format(P.RampTime))
    lines.append("#   delta_eff (plateau + ramp)  = {0} us".format(P.delta + P.RampTime))
    if hasattr(P, "EddyCurrentDelay"):
        lines.append("#   EddyCurrentDelay (LED)      = {0} us".format(P.EddyCurrentDelay))
    lines.append("#   GradSpoil / GradSettle / PreGrad = {0} / {1} / {2} us".format(
        P.GradSpoil, P.GradSettle, P.PreGrad))
    lines.append("#")
    lines.append("# Gradient calibration")
    lines.append("#   Probe                       = {0}".format(meta["probe"]))
    lines.append("#   GradAxis                    = {0}".format(meta["axis"]))
    lines.append("#   max gradient (this axis)    = {0} G/cm  [{1}]".format(
        meta["maxgrad_G_per_cm"], meta["maxgrad_source"]))
    lines.append("#   FP scaler (this axis)       = {0}".format(meta["fp_scaler"]))
    lines.append("#")
    lines.append("# Nucleus")
    lines.append("#   gamma/2pi                   = {0} MHz/T ({1})".format(
        meta["gamma_over_2pi_MHz_per_T"], meta["gamma_source"]))
    lines.append("#   SF                          = {0} MHz".format(base_frequency(P)))
    lines.append("#   implied B0                  = {0:.4f} T".format(meta["implied_B0_T"]))
    lines.append("#")
    lines.append("# Factorisation: b = gamma^2 * G_abs^2 * S0, G_abs = G_rel * FP * GradMax")
    lines.append("#   S0 (timing only)            = {0:.9e} s^3".format(meta["shape_factor_s3"]))
    lines.append("#   b / G_abs^2                 = {0:.9e} s/m^2 per (T/m)^2".format(meta["b_per_G2"]))
    lines.append("#   b / G_rel^2                 = {0:.9e} s/m^2".format(meta["b_per_Grel2"]))
    lines.append("# All three are CONSTANT over this whole ramp. To rebuild every b")
    lines.append("# from the relative-gradient list alone:  b = (b/G_rel^2) * G_rel^2.")
    lines.append("# To rebuild it under a DIFFERENT gradient calibration:")
    lines.append("#   b = (b/G_abs^2) * (G_rel * FP * GradMax_in_T_per_m)^2")
    lines.append("# Note the SQUARE -- b/G_rel (first power) is not constant.")
    lines.append("#")
    lines.append("# b is the full integral of q(t)^2 dt over this sequence's own")
    lines.append("# effective gradient waveform -- trapezoidal ramps, RF pulse widths")
    lines.append("# and frozen-q storage periods all included. Fit ln(S/S0) vs -b.")
    lines.append("#")
    lines.append("# {0:>3s} {1:>10s} {2:>12s} {3:>14s} {4:>16s} {5:>14s}".format(
        "idx", "G_rel", "G_T_per_m", "q_per_m", "b_s_per_m2", "b_s_per_mm2"))
    for idx, rel, G, q, b in rows:
        lines.append("  {0:>3d} {1:>10.5f} {2:>12.6f} {3:>14.6e} {4:>16.8e} {5:>14.6f}".format(
            idx, rel, G, q, b, b * 1.0e-6))
    text = "\n".join(lines) + "\n"
    with open(path, "w", newline="\r\n") as handle:
        handle.write(text)
    return text


# ---------------------------------------------------------------------------

def solve_gmax(meta, slope, known_d):
    """Invert a measured diffusion decay to get the gradient calibration.

    For a ramp in G_rel at fixed timing, ln(S/S0) = -b*D = -(b/G_rel^2)*D*G_rel^2,
    so a straight-line fit of ln(S/S0) against G_rel^2 has

        slope = -(b/G_rel^2) * D

    and since b/G_rel^2 = gamma^2 * S0 * (FP * GradMax_T_per_m)^2, the measured
    slope pins GradMax down directly -- no gradient calibration assumed anywhere.
    This is the honest way to settle the vendor-vs-MAXGRAD_TABLE disagreement:
    run a ramp on a sample whose D you know, fit, and read the answer off here.
    """
    if slope > 0:
        print("NOTE: --solve-gmax expects the slope of ln(S/S0) vs G_rel^2, which "
              "is NEGATIVE for a decay. Using |slope|.")
    measured_b_per_Grel2 = abs(slope) / known_d
    fp = meta["fp_scaler"]
    denom = meta["b_per_G2"] * fp * fp
    gmax_T_per_m = (measured_b_per_Grel2 / denom) ** 0.5
    gmax_G_per_cm = gmax_T_per_m / 0.01
    assumed = meta["maxgrad_G_per_cm"]
    print()
    print("--- gradient calibration implied by your measurement -----------------")
    print("  fitted slope of ln(S/S0) vs G_rel^2 : {0:.6g}".format(slope))
    print("  reference D                         : {0:.6g} m^2/s".format(known_d))
    print("  => measured b/G_rel^2               : {0:.6g} s/m^2"
          .format(measured_b_per_Grel2))
    print("  => implied max gradient             : {0:.4f} G/cm = {1:.5f} T/m"
          .format(gmax_G_per_cm, gmax_T_per_m))
    print("  currently assumed ({0})".format(meta["maxgrad_source"]))
    print("                                      : {0:.4f} G/cm".format(assumed))
    print("  ratio implied/assumed               : {0:.4f}  (b would scale by {1:.4f})"
          .format(gmax_G_per_cm / assumed, (gmax_G_per_cm / assumed) ** 2))
    print("----------------------------------------------------------------------")
    print()


def find_sequence_file(sequence_name, hint_dir=None):
    roots = []
    if hint_dir:
        roots.append(hint_dir)
    here = os.path.dirname(os.path.abspath(__file__))
    roots.append(os.path.join(here, "..", "sequences"))
    roots.append(os.path.join(here, "..", "sequences", "diffusion"))
    for root in roots:
        root = os.path.abspath(root)
        if not os.path.isdir(root):
            continue
        for dirpath, _dirnames, filenames in os.walk(root):
            for name in filenames:
                if name == sequence_name + ".py":
                    return os.path.join(dirpath, name)
    return None


DEFAULT_PAR_ROOT = r"C:\Users\Public\Documents\SpinFlow\Parameters"


def find_par(sequence_name, root=None):
    root = root or DEFAULT_PAR_ROOT
    if not os.path.isdir(root):
        raise SystemExit("Parameters root not found: {0}\nPass --par-root.".format(root))
    hits = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            if name.lower() == (sequence_name + ".par").lower():
                hits.append(os.path.join(dirpath, name))
    if not hits:
        raise SystemExit("No {0}.par found under {1} -- run the sequence in SpinFlow "
                         "once so it saves a parameter file, or use --defaults."
                         .format(sequence_name, root))
    if len(hits) > 1:
        print("NOTE: {0} matches, using the most recently modified:".format(len(hits)))
        for h in hits:
            print("   ", h)
    return max(hits, key=os.path.getmtime)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("par", nargs="?", help="path to a saved SpinFlow .par file")
    ap.add_argument("--find", metavar="SEQUENCE",
                    help="search the Parameters tree for <SEQUENCE>.par")
    ap.add_argument("--par-root", default=DEFAULT_PAR_ROOT,
                    help="root of the SpinFlow Parameters tree (default: %(default)s)")
    ap.add_argument("--sequence", metavar="PATH",
                    help="pulse programme .py to take the maths and defaults from "
                         "(default: found by name next to this script)")
    ap.add_argument("--defaults", action="store_true",
                    help="ignore any .par and use the sequence's own defaults")
    ap.add_argument("--grad-list", help="comma-separated relative gradients, "
                                        "overriding the GradList Parameter")
    ap.add_argument("--probe", help="override the Probe Parameter")
    ap.add_argument("--nucleus", help="override the Nucleus Parameter")
    ap.add_argument("--gmax", type=float, metavar="G_PER_CM",
                    help="override the max gradient for the selected axis [G/cm at "
                         "|G1|=1.0], instead of the Probe table or the GradMax "
                         "Parameter -- use to re-cost an experiment under a "
                         "different gradient calibration")
    ap.add_argument("--solve-gmax", type=float, metavar="SLOPE",
                    help="work BACKWARDS from a measurement: given the fitted slope "
                         "of ln(S/S0) vs G_rel^2 (negative) and --known-d, print the "
                         "gradient calibration that measurement implies")
    ap.add_argument("--known-d", type=float, metavar="D",
                    help="diffusion coefficient of the reference sample [m^2/s], "
                         "e.g. 2.299e-9 for pure H2O at 25 C -- used with --solve-gmax")
    ap.add_argument("--out", help="output path (default: <Sequence>.diff beside the .par)")
    ap.add_argument("--difflist", action="store_true",
                    help="also write <Sequence>.difflist (bare relative gradients)")
    ap.add_argument("--show-par", action="store_true",
                    help="print every key recovered from the .par and what it mapped to")
    args = ap.parse_args(argv)

    par_path = None
    if args.find:
        sequence_name = args.find
        if not args.defaults:
            par_path = find_par(sequence_name, args.par_root)
    elif args.par:
        par_path = os.path.abspath(args.par)
        if not os.path.isfile(par_path):
            raise SystemExit("No such file: " + par_path)
        sequence_name = os.path.splitext(os.path.basename(par_path))[0]
    elif args.sequence:
        sequence_name = os.path.splitext(os.path.basename(args.sequence))[0]
    else:
        ap.error("give a .par path, or --find <SEQUENCE>, or --sequence <PATH> --defaults")

    seq_path = args.sequence or find_sequence_file(sequence_name)
    if not seq_path or not os.path.isfile(seq_path):
        raise SystemExit("Could not find the pulse programme for '{0}'. Pass "
                         "--sequence <path to {0}.py>.".format(sequence_name))
    seq_path = os.path.abspath(seq_path)
    print("sequence : {0}".format(seq_path))
    mod = load_sequence(seq_path)

    if par_path and not args.defaults:
        print("par file : {0}".format(par_path))
        par_values = parse_par(par_path)
        if not par_values:
            print("  WARNING: no key/value pairs recovered from the .par -- falling "
                  "back to the sequence's own defaults. Re-run with --show-par.")
        P = overlay(mod.Parameters, par_values, verbose=args.show_par)
        source_note = "parameters from " + par_path
        out_dir = os.path.dirname(par_path)
    else:
        print("par file : (none -- using the sequence's own defaults)")
        P = overlay(mod.Parameters, {}, verbose=False)
        source_note = "sequence defaults (no .par applied)"
        out_dir = os.path.dirname(seq_path)

    if args.probe:
        P.Probe = args.probe
    if args.nucleus:
        P.Nucleus = args.nucleus
    if args.gmax is not None:
        P.GradMaxOverride = args.gmax
    grad_list = None
    if args.grad_list:
        grad_list = [float(x) for x in args.grad_list.split(",") if x.strip()]

    rows, meta = build_table(mod, P, grad_list)

    if args.solve_gmax is not None:
        if not args.known_d:
            raise SystemExit("--solve-gmax also needs --known-d (the reference "
                             "sample's diffusion coefficient, m^2/s).")
        solve_gmax(meta, args.solve_gmax, args.known_d)

    out_path = args.out or os.path.join(out_dir, sequence_name + ".diff")
    text = write_diff(out_path, sequence_name, rows, meta, P, source_note)
    print("wrote    : {0}".format(out_path))
    print()
    print(text)

    if args.difflist:
        list_path = os.path.splitext(out_path)[0] + ".difflist"
        with open(list_path, "w", newline="\r\n") as handle:
            handle.write("\n".join("{0}".format(r[1]) for r in rows) + "\n")
        print("wrote    : {0}".format(list_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
