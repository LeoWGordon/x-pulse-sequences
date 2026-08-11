#-------------------------------------------------------------------------------
# Name:        LWG_CPMG-Image-Echo_H.py
# Purpose:     Multi-echo (CPMG) 1D MRI / image-echo sequence for the Oxford
#              Instruments X-Pulse Broadband Benchtop NMR Spectrometer
#              (1H/19F channel). Applies NECH refocusing pulses per scan and
#              images ONLY the last one -- i.e. a full T2(x) decay curve is
#              built by running this sequence as a "quasi-2D" experiment,
#              sweeping NECH across a list of separate scans, exactly like
#              the vendor's own (non-imaging) CPMG_H.py. Built to answer:
#              "does the effective-T2 gradient near an interface come from
#              T2 or from diffusion?" (see chat discussion, 2026-08-11) -- a
#              genuine CPMG T2 series with short, FIXED echo spacing
#              (2*TAU), reached by increasing NECH rather than TAU, which
#              suppresses diffusion attenuation relative to the existing
#              diffprof.lwg-style single-echo TAU array (where increasing
#              TAU directly, rather than the echo count, is what increases
#              TE -- see LWG_1D-Image-Echo-GradAfterRefocus_H.py's own
#              "IMPLICATIONS" section for the diffusion-weighting mechanism
#              this is trying to suppress).
#
# Author:      Claude for L. Gordon (DTU)
#
# Created:     11/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     2.0 (draft) -- REVISED from v1.0 after L. Gordon supplied the
#              vendor's own CPMG_H.py as a reference (see "WHAT CHANGED IN
#              v2.0" below). v1.0's core mechanism (multiple Receiver1()
#              calls within one scan, to acquire every echo from a single
#              excitation) is NOT how CPMG_H.py -- the one real, presumably-
#              running multi-echo sequence available as a reference on this
#              instrument -- does it, and v1.0 is not carried forward.
# Status:      draft -- NOT YET RUN ON HARDWARE. See "BEFORE YOU RUN THIS"
#              below and docs/validation-notes/LWG_CPMG-Image-Echo_H.md.
#              Per this repo's README: no AI-assisted sequence is used for
#              real data collection until a named person has reviewed and
#              validated it on the instrument.
#
# ===============================================================================
# WHAT CHANGED IN v2.0 (vs the original v1.0 draft)
# ===============================================================================
#
# v1.0 called Receiver1() NECH times within a single scan, betting that the
# framework would deliver one NECH*ReceiverPoints array per scan (echo-
# major). There was no precedent for this anywhere in this repo. The
# vendor's own CPMG_H.py (supplied by L. Gordon, 2026-08-11) resolves this
# directly by demonstrating the ACTUAL supported pattern: apply NECH
# refocusing pulses in a plain loop, call Receiver1() exactly ONCE, right
# after the LAST one -- only the final echo is ever acquired in a given
# scan. Per CPMG_H.py's own header comment, the full decay curve is
# obtained by running the sequence as a "quasi 2D experiment (FFT data
# [ppm] vs. EchoNumber)" -- i.e. a VC-list/array sweep of NECH across
# separate scans, structurally the same "sweep one parameter across
# separate acquisitions" idea already used by diffprof.lwg's TAU array.
#
# This version follows that proven pattern instead: NECH-1 plain "blind"
# refocusing pulses (no gradient, no acquisition -- single sequential
# timeline, exactly like CPMG_H.py's own loop, including its exact
# Delay() formulas and its 8-step Meiboom-Gill-style phase cycle for P1/
# P180/receiver, all reused verbatim) followed by ONE final refocusing
# pulse that -- unlike the blind ones -- gets the gradient dephase/readout
# treatment and the acquisition, structured exactly like
# LWG_1D-Image-Echo-GradAfterRefocus_H.py's own single echo (diffusion-
# compensated gradient geometry, same post-refocus wait formula, which is
# NOT the same formula CPMG_H.py uses for its own pre-acquisition wait --
# see "WHY THE FINAL ECHO DOESN'T REUSE CPMG_H.py's OWN PRE-ACQUISITION
# FORMULA" below for why these two proven formulas were kept separate
# rather than blended). This removes v1.0's entire multi-Receiver1-per-scan
# assumption AND its hand-derived inter-echo branch-padding math -- there
# is now only ONE parallel (gradient vs RF) block per scan, for the final
# echo only, exactly as validated in the imaging family already.
#
# WHY THE FINAL ECHO DOESN'T REUSE CPMG_H.py's OWN PRE-ACQUISITION FORMULA:
# CPMG_H.py's post-refocus wait (Tau - P2/2 - dead_time + group_delay - 3)
# is tuned for a plain spectroscopy FID acquisition (compensating the
# receiver filter's group delay so the FID itself starts cleanly) -- it
# has no reason to centre a GRADIENT echo (k=0) within the acquisition
# window, because there is no gradient in that sequence. Blending that
# formula with this sequence's gradient readout lobe would not obviously
# still centre the gradient echo correctly (worked through by hand while
# writing this -- the two formulas' implicit alignment targets don't
# obviously coincide for arbitrary Dead1/group_delay/RampTime values).
# Instead, the final echo's post-refocus RF wait reuses
# LWG_1D-Image-Echo-GradAfterRefocus_H.py's OWN formula (Tau - P180/2 +
# PreGrad - dead_time), which IS specifically derived to align with this
# same gradient lobe geometry -- already validated in that file. Net
# effect: the NECH-1 blind pulses are timed exactly like CPMG_H.py; the
# final, imaged echo is timed exactly like the existing single-echo
# imaging sequences. Nothing here blends the two references' formulas.
#
# ===============================================================================
# BEFORE YOU RUN THIS -- read in full
# ===============================================================================
#
# (1) Run as a quasi-2D / VC-list sweep of NECH (e.g. 1, 2, 4, 8, 16, ...),
#     exactly per CPMG_H.py's own header note -- this sequence images only
#     the LAST of NECH refocusing pulses each scan, on purpose (see "WHAT
#     CHANGED IN v2.0"). A single run at fixed NECH gives you one profile
#     at one TE = (see comms.log at run start for the exact value, same
#     "echo time isn't just 2*TAU*NECH" caveat as v1.0 -- real hardware
#     overhead adds to it).
#
# (2) UNVERIFIED: the NECH-1 blind refocusing pulses' timing formulas are
#     copied verbatim from CPMG_H.py (a real vendor sequence), so they
#     should be reliable AS TIMING FORMULAS -- but this file has still
#     never been compiled or run, and the transition point (going from the
#     last blind pulse into the final, gradient-carrying echo) is new,
#     not present in CPMG_H.py. Check gradient/RF synchronisation on the
#     final echo specifically (e.g. by comparing the profile shape/
#     position against a same-Parameters single-echo run of
#     LWG_1D-Image-Echo-GradAfterRefocus_H.py at NECH=1, where the two
#     sequences' final-echo timing should be nearly identical) before
#     trusting NECH>1 data.
#
# (3) DUTY CYCLE: RF duty cycle now scales with NECH (P90 once + P180 x
#     NECH per scan); gradient duty cycle does NOT scale with NECH (only
#     the final echo ever drives the gradient coil) -- see
#     estimate_duty_cycles() below, which reflects this asymmetry
#     (different from v1.0, which scaled both with NumEchoes).
#
# Basis:       NECH-1 blind-pulse loop and its phase-cycle strings (PH1/
#              PH2/PHRX) copied verbatim from the vendor's CPMG_H.py
#              (Created 10/09/2013, Author AS). Final echo's gradient
#              geometry and RF timing copied verbatim from
#              LWG_1D-Image-Echo-GradAfterRefocus_H.py (itself derived from
#              LWG_1D-Image-Echo_H.py v3.1). Non-timing helpers
#              (apply_gradient, pulse, safe_delay, mains_lock_trigger,
#              estimate_duty_cycles, report_probe_gradient, MAXGRAD_TABLE,
#              get_permutation) carried over from the imaging family
#              unchanged. No composite (90-270) refocusing pulse is used
#              anywhere in this file -- CPMG_H.py uses a plain single P180
#              for every refocusing pulse, including the last, so this
#              file does too (a deliberate divergence from the rest of the
#              imaging family, which uses a composite pulse for its single
#              refocusing event).
#-------------------------------------------------------------------------------

from firebird import *
from firebird.applications import *
import numpy as np
import time
global BLP

BLP = 0

def get_permutation(ds, ns, nb):
    first_scans = []
    scans_to_save = []
    for n in range(ds+ns):
        if (n-ds)%(ns/nb) == 0 and n >=ds-1:
            first_scans.append(n)
        if (n-ds+1)%(ns/nb) == 0 and n >ds-1:
            scans_to_save.append(n)
    return first_scans, scans_to_save

class CallBack1D(object):
#   Callback to handle one-dimensional sequences. Each scan returns exactly
#   ReceiverPoints samples (the FINAL echo only) -- same shape/semantics as
#   every other sequence in this repo, NOT NECH*ReceiverPoints (see "WHAT
#   CHANGED IN v2.0" -- unlike v1.0, this version images only one echo per
#   scan, so no multi-echo reshaping is needed here).
    def __init__(self, Params, comms, ReceiverFilter, jcamp_meta):

        self.comms = comms                                      # Comms to the clients
        self.P = Params
        self.NS = self.P.NumScans
        self.NBlock = self.P.NumScans                           # default NBLOCK to 1
        self.jcamp_meta = jcamp_meta
        self.start_time = time.time()
        self.end_time = 0
        if getattr(Params, "NBlock", 'nofound') != "nofound":
            self.NBlock = getattr(Params, "NBlock")             # if defined, set to the value defined
        self.DS = 0                                             # default DS to 0
        if getattr(Params, "DS", 'nofound') != "nofound":
            self.DS = getattr(Params, "DS")                     # if defined, set to the value defined
        self.rcv_filter = ReceiverFilter                        # receiver filter
        if getattr(Params, "TD2", 'nofound') != "nofound":
            Points = getattr(Params, "TD2")
        elif getattr(Params, "ReceiverPoints", 'nofound') != "nofound":
            Points = getattr(Params, "ReceiverPoints")
        else:
            comms.log("{0}".format(getattr(Params, "ReceiverPoints", 'nofound')))

        self.times = np.arange(0.0, Points*self.rcv_filter.dwell, self.rcv_filter.dwell) / 1.0e6

        self.dummys = range(self.DS)                            # Dummy scan
        self.first_scans, self.scan_to_save = get_permutation(self.DS, self.NS, self.NBlock)
        self.final_scan = (self.DS + self.NS) -1                # Last scan in sequence

    def process_data(self, scan, data):
        scan_data, clipped = data
        scaled_data = np.int32(np.round(scan_data / self.rcv_filter.gain))

        if (clipped.value == 1):    #clipped
            self.comms.log(">>>>> !! Data is Clipped!! <<<<<")
        elif (clipped.value == 100):
            self.comms.log(">>>>> !! IPC failed!! <<<<<")

        self.acc_data = np.array(scaled_data)

        if scan in self.dummys:                                                   # sort out dummy scan
            self.comms.send_data(self.acc_data, scan+1, self.times, is_clipped=clipped,
                is_last_update=False, metadata={"DummyScan":"true"})
        elif (scan not in self.scan_to_save):
            self.comms.send_data(self.acc_data, scan + 1, self.times, is_clipped=clipped,
                                 metadata={"JCAMP":self.jcamp_meta, "nosave":"true","Preview":2})
            if scan==self.final_scan:
                self.comms.send_data(self.acc_data, scan + 1, self.times, is_clipped=clipped, is_last_update=(scan==self.final_scan),
                                metadata={"JCAMP":self.jcamp_meta})
        else: #scans_per_block'th scan, This data will be saved accumulator will be reset
            self.comms.send_data(self.acc_data, scan + 1, self.times, is_clipped=clipped, is_last_update=(scan==self.final_scan),
                            metadata={"JCAMP":self.jcamp_meta, "Preview":2})
            self.acc_data = None
        if scan == (self.DS+self.NS-1):
            self.end_time = time.time()
            self.comms.log("seqTime LWG_CPMG-Image-Echo_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


def check_range(pvalue, **kwargs):
    params = {}
    for key, value in kwargs.items():
        params[key] = value
    if params and "min" in params and "max" in params:
        if pvalue < params['min'] or pvalue > params['max']:
            return False, "Out of range."
    return True, "OK"

RD = check_range
RFA = check_range


def _lobe_elapsed(P, delta):
    """Elapsed time of a single apply_gradient() lobe call with plateau
    duration 'delta' -- see LWG_1D-Image-Echo-GradAfterRefocus_H.py."""
    return P.PreGrad + 2.0*P.RampTime + delta + P.GradSettle


def _first_gap(P):
    """Delay between P1 (excitation) and the first refocusing pulse.
    Verbatim from CPMG_H.py."""
    return P.Tau - P.TXEnableTime - ((P.P90+P.P180)/2.0) - 12.0


def _subsequent_gap(P):
    """Delay between one refocusing pulse and the next. Verbatim from
    CPMG_H.py."""
    return 2.0*P.Tau - P.TXEnableTime - P.P180 - 12.0


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)
    N = int(P.EchoNumber)

    # NECH-1 blind pulses, each (gap + pulse). Gap formula is _first_gap
    # only for the very first pulse (N==1 case: that pulse IS the final,
    # imaged one, handled below instead); every other blind pulse uses
    # _subsequent_gap.
    blind_time = 0.0
    if N > 1:
        blind_time += _first_gap(P) + (P.TXEnableTime + P.P180)
        blind_time += (N-2) * (_subsequent_gap(P) + (P.TXEnableTime + P.P180))

    final_gap = _first_gap(P) if N == 1 else _subsequent_gap(P)

    # Final echo body (post-refocus wait + ACQU on the RF branch; dephase +
    # inter-lobe idle + readout on the gradient branch) -- same closed-form
    # simplification as LWG_1D-Image-Echo-GradAfterRefocus_H.py.
    RefocusPulseWidth = P.TXEnableTime + P.P180  # plain pulse(), not composite
    B = P.Tau - P.P180/2.0 + P.PreGrad - ReceiverFilter.dead_time
    T_rf_body = B + 1.0 + 1.0 + ReceiverFilter.dead_time + P.Dead1 + ReceiverTime
    T_grad_body = P.Tau + 4.0*P.RampTime + 2.0*P.GradientOnTime + P.GradSettle
    final_echo_time = final_gap + RefocusPulseWidth + max(T_rf_body, T_grad_body)

    t_pulse_gradient_event = (P.TXEnableTime + P.P90) + blind_time + final_echo_time

    t_scanTime = P.RecycleDelay + 200.0 + t_pulse_gradient_event
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime


def sequence_description():

    seq_desc = ("Multi-echo (CPMG) 1D MRI / image-echo sequence on the "
                "{H/F} channel. Applies NECH refocusing pulses per scan "
                "(blind-pulse loop timed exactly per the vendor's own "
                "CPMG_H.py) and images only the LAST one, with a "
                "diffusion-compensated gradient readout (see "
                "LWG_1D-Image-Echo-GradAfterRefocus_H.py). Run as a "
                "quasi-2D sweep of NECH to build a T2(x) decay curve with "
                "short, fixed echo spacing. DRAFT, not yet run on "
                "hardware -- see file header before use.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,NECH,D70,D73,D71,D75,G1,GradAxis"

    return basic

def get_gradient_functions(axisstr):
    """Map axis number to gradient functions."""
    if axisstr == 'x':
        axis = 1
    elif axisstr == 'y':
        axis = 2
    elif axisstr == 'z':
        axis = 3
    else:
        exit()

    gradient_map = {
        1: (Gradient1SlewRate, Gradient1),
        2: (Gradient2SlewRate, Gradient2),
        3: (Gradient3SlewRate, Gradient3),
    }
    return gradient_map[axis]

def apply_gradient(axis, delta, ramp_time, gradient_value, settle_time, pre_grad):
    """Apply gradient pulses to specified axis. gradient_value is the final
    logical (normalised, -1..1) amplitude sent to hardware -- axis
    calibration (FPX/FPY/FPZ) is applied exactly once, by GradientMatrix()
    in run() below, NOT here."""
    Delay(pre_grad)
    gradient_slew_rate = abs(gradient_value / ramp_time)
    slew_rate_fn, gradient_fn = get_gradient_functions(axis)
    slew_rate_fn(1.0, gradient_slew_rate)
    gradient_fn(ramp_time, gradient_value)
    Delay(delta)
    gradient_fn(ramp_time, 0)
    gradient_fn(1.0, 0)
    Delay(settle_time)

def pulse(length, phase, txenabletime):
    """Apply a pulse of given length and phase."""
    Channel1SetBasePhase(1, phase)
    Transmit1BlankingOn(1)
    Delay(txenabletime)
    Transmit1(length)
    Transmit1BlankingOff(1)


def safe_delay(value, label, comms):
    """Delay() wrapper for every timing-formula-derived (as opposed to
    fixed) delay in the sequence. Raises a clear, specific error (naming
    the exact wait) instead of a silent/confusing failure if the formula
    goes negative."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU and/or PreGrad, or decrease "
               "GradientOnTime/RampTime/GradSettle, then retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger before the first pulse event.
    OFF by default (UseMainsLock=0) -- see LWG_1D-Image-Echo_H.py for full
    reasoning (mains triggering adds a non-deterministic pre-scan delay
    that works against reproducible gradient-echo timing)."""
    if int(P.UseMainsLock) == 0:
        return
    triggers = {1: ExternalTrigger1, 2: ExternalTrigger2, 3: ExternalTrigger3}
    ch = int(P.MainsLockChannel)
    if ch not in triggers:
        comms.log("WARNING: MainsLockChannel={0} is not 1, 2, or 3 -- mains "
                  "lock trigger skipped.".format(ch))
        return
    triggers[ch]()


def estimate_duty_cycles(P, comms):
    """Log estimated RF and gradient duty cycles for this scan and warn if
    they exceed MaxRFDuty/MaxGradDuty.

    UNLIKE v1.0: RF duty cycle scales with NECH (P90 once + P180 x NECH,
    every refocusing pulse transmits regardless of whether it's imaged).
    Gradient duty cycle does NOT scale with NECH -- only the final,
    imaged echo ever drives the gradient coil, exactly one dephase +
    readout lobe pair per scan regardless of how many blind pulses
    preceded it.

    MaxRFDuty/MaxGradDuty are the SAME conservative, NOT vendor-confirmed
    placeholders used throughout this repo -- confirm real limits with
    Oxford Instruments, especially for large NECH.
    """
    TR = float(P.RecycleDelay)
    N = int(P.EchoNumber)
    rf_on_time = P.P90 + N*P.P180
    grad_on_time = P.GradientOnTime + (P.GradientOnTime*2 + P.RampTime*2)

    rf_duty = rf_on_time / TR
    grad_duty = grad_on_time / TR

    comms.log("Estimated duty cycle (NECH={0}) -- RF: {1:.3%} (limit {2:.1%}), "
              "Gradient [{3}-axis, final echo only]: {4:.3%} (limit {5:.1%})".format(
              N, rf_duty, P.MaxRFDuty, P.Axis, grad_duty, P.MaxGradDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} (P90 + {1} x P180) exceeds "
                  "MaxRFDuty ({2:.2%}). Consider a longer RD, fewer echoes, or shorter "
                  "pulses/TXEnable, or confirm with Oxford Instruments that this is "
                  "within the transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, N, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} on the {1}-axis coil "
                  "exceeds MaxGradDuty ({2:.2%}). Consider a longer RD or shorter "
                  "GradientOnTime, or confirm the rated duty cycle for the X-Pulse "
                  "gradient amplifier/coil with Oxford Instruments."
                  .format(grad_duty, P.Axis, P.MaxGradDuty))


# Per-axis max gradient strength (G/cm). Kept in sync with
# LWG_1D-Image-Echo_H.py / leonmr/xpulse_imaging.py's own tables.
MAXGRAD_TABLE = {
    "HFX": {"x": 11.879, "y": 11.978, "z": 57.915},   # G/cm, measured/averaged 06/08/2026
    "LOWGAMMA": {"x": None, "y": None, "z": None},     # NOT YET CALIBRATED
}

def report_probe_gradient(P, comms, axis=None):
    """Log the max gradient strength (G/cm) for the currently-selected Probe
    and gradient axis -- see LWG_1D-Image-Echo_H.py for full rationale."""
    probe_key = str(P.Probe).upper().replace('/', '').replace('-', '').replace(' ', '').replace('_', '')
    axis_key = str(axis if axis is not None else getattr(P, "Axis", "z")).strip().lower()
    if probe_key not in MAXGRAD_TABLE:
        comms.log("WARNING: unrecognised Probe='{0}'. Known probes: {1}. "
                  "Hz->mm conversion will not know this probe's gradient "
                  "calibration.".format(P.Probe, list(MAXGRAD_TABLE)))
        return
    maxgrad = MAXGRAD_TABLE[probe_key].get(axis_key)
    if maxgrad is None:
        comms.log("WARNING: Probe='{0}', axis='{1}' has NO gradient "
                  "calibration yet -- Hz->mm conversion is undefined for "
                  "this data until it is measured.".format(P.Probe, axis_key))
    else:
        comms.log("Probe: {0}. Max gradient strength ({1}-axis) = {2} G/cm."
                  .format(P.Probe, axis_key, maxgrad))


@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_CPMG-Image-Echo_H", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 59.7, ParameterTypes.Double, "H/F Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "H/F Freq Offset [Hz]")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "H/F TX Freq Offset [ppm]")

    ReceiverPoints = Parameter("NP", 1024, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 34, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "100000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans
    NumScans = Parameter("NS", 4, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 0, ParameterTypes.Int32,"Dummy Scans")

    # CPMG train -- run as a quasi-2D sweep of NECH, see file header
    EchoNumber = Parameter("NECH", 8, ParameterTypes.Int32, "Refocusing pulses/scan [last one imaged]")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 100.0, ParameterTypes.Double, "H/F Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 500000, ParameterTypes.Int32, "Relaxation Delay [s]",
                              RD, min=100000, max=2000000000)

    # Hard pulses -- P180 (plain 180, not the composite 90-270 used
    # elsewhere in the imaging family) is used for EVERY refocusing pulse,
    # blind or final, per CPMG_H.py's own convention.
    P90 = Parameter("P90", 9.58, ParameterTypes.Double, "H/F 90&#176; Pulse Width [&#956;s]")
    P180 = Parameter("P180", 19.16, ParameterTypes.Double, "H/F 180&#176; Pulse Width [&#956;s]")
    TXAmplitude = Parameter("RFA0", 0.4, ParameterTypes.Double,
                            "H/F TX Power [0.0&#8230;1.0]", RFA, min=0.0, max=1.0)

    # Gradients -- FPX/FPY/FPZ are the ONLY place per-axis calibration is
    # applied (via GradientMatrix() in run()); G1 is the single logical
    # gradient-strength knob, used only on the final (imaged) echo.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet [HFX(default, calibrated)/LOWGAMMA(not yet calibrated)]")

    # Sequence specific -- same meaning/units as LWG_1D-Image-Echo-
    # GradAfterRefocus_H.py, applied only to the final echo.
    GradientOnTime = Parameter("D71", 1000.0, ParameterTypes.Double, "Gradient Duration [&#956;s] -- dephase lobe plateau (readout = 2*this+2*RampTime), final echo only")
    RampTime = Parameter("D70", 100.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s] -- both lobes incl. polarity reversal")
    Tau = Parameter("TAU", 3000.0, ParameterTypes.Double, "&#964; Delay [&#956;s] -- half the echo spacing, same convention as CPMG_H.py")
    PreGrad = Parameter("D75", 100.0, ParameterTypes.Double, "Pre-Gradient Time [&#956;s]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")

    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rails (conservative, user-adjustable placeholders --
    # see estimate_duty_cycles() for how these are computed for this
    # sequence specifically)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Phases -- verbatim from CPMG_H.py's own 8-step Meiboom-Gill-style
    # cycle. PH2 (P180) is applied identically to EVERY refocusing pulse
    # within a scan (blind and final alike) -- standard CPMG: constant
    # refocus phase relative to excitation across the whole train.
    PH1 = Parameter("PH1", "0,0,180,180,90,90,270,270", ParameterTypes.String, "H/F 90&#176; Pulse Phase")
    PH2 = Parameter("PH2", "90,270,90,270,0,180,0,180", ParameterTypes.String, "H/F 180&#176; Pulse Phase")
    PHRX = Parameter("PHRX", "0,0,180,180,90,90,270,270", ParameterTypes.String, "Acquisition Phase")

def run(comms):

    P = Parameters
    N = int(P.EchoNumber)
    if N < 1:
        raise ValueError("EchoNumber (NECH) must be >= 1.")

    Phases = PhasesManager(P)
    Phases.Reset()

    FPX = P.XGradNorm
    FPY = P.YGradNorm
    FPZ = P.ZGradNorm

    Matrix = np.array([
                        [FPX,  0,  0],
                        [0  ,FPY,  0],
                        [0  ,  0,FPZ]
                        ])

    Frequency = P.FrequencyBase + (P.FrequencyOffset*1.0e-6) + (P.FrequencyBase*P.TxPPM*1.0e-6)
    Channel1SetFrequency(10,Frequency)
    Channel1RestartSynth(10)

    ReceiverFilter = Filter(P.Filter)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)

    final_gap = _first_gap(P) if N == 1 else _subsequent_gap(P)
    RefocusPulseWidth = P.TXEnableTime + P.P180

    T_rf_body = (P.Tau - P.P180/2.0 + P.PreGrad - ReceiverFilter.dead_time) \
                + 1.0 + 1.0 + ReceiverFilter.dead_time + P.Dead1 + ReceiverTime
    T_grad_body = P.Tau + 4.0*P.RampTime + 2.0*P.GradientOnTime + P.GradSettle
    final_echo_TE = final_gap + RefocusPulseWidth + max(T_rf_body, T_grad_body)
    if N > 1:
        prior_TE = _first_gap(P) + (P.TXEnableTime+P.P180) \
                   + (N-2)*(_subsequent_gap(P) + (P.TXEnableTime+P.P180)) if N > 1 else 0.0
    else:
        prior_TE = 0.0
    total_TE = (P.TXEnableTime + P.P90) + prior_TE + final_echo_TE
    comms.log("CPMG-Image-Echo: NECH={0}, imaged echo TE ~= {1:.1f} us ({2:.3f} ms) "
              "relative to P1.".format(N, total_TE, total_TE/1000.0))

    def jcamp_meta():
        global BLP
        if BLP == 1:
            me_mod = 4
            jc_blp = -2*int(round((((ReceiverFilter.dead_time/2.)+P.Dead1)/ReceiverFilter.dwell),0))
        else:
            me_mod = 0
            jc_blp = 0
        jc_data = ["jc_sf={0}".format(Frequency),
                "jc_sf_txppm={0}".format(Frequency+(P.FrequencyBase*P.TxPPM*1.0e-6)),
                "ME_mod={0}".format(me_mod),
                "jc_blp={0}".format(jc_blp),
                "jc_ns=1",
                "EchoNumber={0}".format(N),
                "ImagedEchoTEus={0:.3f}".format(total_TE),
                    ]
        return jc_data

    comms.log("Sequence Time: {0}".format(time_calculation(Parameters)))
    seqAcqu = CallBack1D(Parameters, comms, ReceiverFilter, jcamp_meta())

    def RecvCallback(acqData, scan):
        seqAcqu.process_data(scan, acqData)

    # ---- Read-gradient coverage check (warning only, not a hard stop) ----
    decode_plateau = P.GradientOnTime*2 + P.RampTime*2
    if decode_plateau < ReceiverTime:
        comms.log("WARNING: readout gradient plateau (~{0:.1f} us, from "
                  "GradientOnTime/RampTime) is shorter than the acquisition "
                  "window (ReceiverTime = {1:.1f} us, from NP*dwell). Part "
                  "of the FID will be sampled while the read gradient is "
                  "off or ramping -- increase GradientOnTime (D71), reduce "
                  "NP, or use a wider Filter to avoid truncation/leakage."
                  .format(decode_plateau, ReceiverTime))

    with sequential:

        mains_lock_trigger(P, comms)

        Transmit1SelectPort(1,1)
        Transmit1LPEnable(1,0)
        Delay(10000) # This is required for changing the relay state from tune mode
        Transmit1SetScale(5, P.TXAmplitude)

        Receiver1Preamp(128, P.ReceiverAttenuation)
        Receiver1Filter(200, ReceiverFilter)

        # Gradient Setup -- the ONLY place FPX/FPY/FPZ scale the gradients
        GradientMatrix(200,Matrix.T)

        Phases.Reset()

    report_probe_gradient(P, comms)
    estimate_duty_cycles(P, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver1FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            # ---- P1 excitation ----
            pulse(P.P90, ph["PH1"], P.TXEnableTime)

            # ---- NECH-1 blind refocusing pulses (no gradient, no
            # acquisition) -- single sequential timeline, timed exactly
            # per CPMG_H.py. Skipped entirely when NECH==1 (that single
            # pulse is then the final, imaged one, handled below).
            for echo_idx in range(N-1):
                gap = _first_gap(P) if echo_idx == 0 else _subsequent_gap(P)
                safe_delay(gap, "blind inter-pulse gap, echo {0}".format(echo_idx+1), comms)
                pulse(P.P180, ph["PH2"], P.TXEnableTime)

            # ---- Final refocusing pulse: gradient dephase/readout +
            # acquisition, exactly like LWG_1D-Image-Echo-
            # GradAfterRefocus_H.py's own single echo.
            with parallel:
                with sequential:
                    # Gradient branch: idle through the gap AND the pulse
                    # itself (matching the RF branch's elapsed time up to
                    # the end of the final refocusing pulse), then dephase
                    # (post-refocus, NEGATIVE polarity) + inter-lobe idle +
                    # readout (POSITIVE polarity).
                    safe_delay(final_gap + P.TXEnableTime + P.P180,
                               "gradient pre-final-pulse wait", comms)
                    apply_gradient(P.Axis, P.GradientOnTime, P.RampTime, -P.G1, P.GradSettle, P.PreGrad)
                    safe_delay((P.Tau - P.PreGrad) - _lobe_elapsed(P, P.GradientOnTime),
                               "gradient inter-lobe wait", comms)
                    apply_gradient(P.Axis, P.GradientOnTime*2+P.RampTime*2, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                with sequential:
                    safe_delay(final_gap, "final inter-pulse gap", comms)
                    pulse(P.P180, ph["PH2"], P.TXEnableTime)
                    # Post-refocus wait -- LWG_1D-Image-Echo-GradAfterRefocus_H.py's
                    # own formula (gradient-echo-aligned), NOT CPMG_H.py's
                    # group_delay-based formula -- see "WHY THE FINAL ECHO
                    # DOESN'T REUSE CPMG_H.py's OWN PRE-ACQUISITION FORMULA"
                    # in the file header.
                    safe_delay(P.Tau - P.P180/2.0 + 1.*P.PreGrad - ReceiverFilter.dead_time,
                               "RF post-refocus wait", comms)
                    Channel1SetBasePhase(P.TXEnableTime,0)
                    Receiver1Phase(P.TXEnableTime, ph["PHRX"])
                    Delay(P.Dead1)
                    Delay(ReceiverFilter.dead_time)
                    Receiver1(ReceiverTime, points)

        Delay(9.0e4)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        TX0.setup_receive(points,32,scan=seqScans,circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX0.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))
        comms.log("NECH %s, Encode Time %s" % (N, P.GradientOnTime))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 11/08/26 (v1.0, draft, SUPERSEDED) - Initial version. Called
#      Receiver1() NumEchoes times per scan, acquiring every echo from a
#      single excitation. No precedent for this in the repo; not carried
#      forward once a real reference (CPMG_H.py) became available. Kept in
#      git history, not in this file.
# 2. Claude - 11/08/26 (v2.0, draft, THIS VERSION) - Rebuilt around the
#      vendor's own CPMG_H.py (supplied by L. Gordon), which images only
#      the LAST of NECH refocusing pulses per scan and gets the full decay
#      curve from a quasi-2D sweep of NECH across separate scans. NECH-1
#      blind pulses' timing/phase-cycling copied verbatim from CPMG_H.py;
#      final (imaged) echo's gradient geometry and RF timing copied
#      verbatim from LWG_1D-Image-Echo-GradAfterRefocus_H.py -- see "WHAT
#      CHANGED IN v2.0" and "WHY THE FINAL ECHO DOESN'T REUSE CPMG_H.py's
#      OWN PRE-ACQUISITION FORMULA" at the top of this file for why these
#      two references were kept separate rather than blended. Removes
#      v1.0's multi-Receiver1-per-scan assumption and inter-echo branch-
#      padding math entirely -- there is now only one parallel
#      (gradient vs RF) block per scan. STILL NOT RUN ON HARDWARE -- see
#      "BEFORE YOU RUN THIS" at the top of this file.
#
# -----------------------------------------------------------------------------
