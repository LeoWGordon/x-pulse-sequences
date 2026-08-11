#-------------------------------------------------------------------------------
# Name:        LWG_CPMG-Image-Echo_H.py
# Purpose:     Multi-echo (CPMG) 1D MRI / image-echo sequence for the Oxford
#              Instruments X-Pulse Broadband Benchtop NMR Spectrometer
#              (1H/19F channel). Acquires NumEchoes gradient-refocused spin
#              echoes -- i.e. NumEchoes separate spatial profiles at
#              increasing echo time -- from a SINGLE excitation, instead of
#              LWG_1D-Image-Echo_H.py's one echo per scan. Built specifically
#              to answer: "does the effective-T2 gradient near an interface
#              come from T2 or from diffusion?" (see chat discussion,
#              2026-08-11) -- a genuine multi-echo T2 train, with each echo's
#              own gradient dephase/readout lobes kept diffusion-COMPENSATED
#              (both lobes on the same, post-refocus side, back-to-back --
#              see LWG_1D-Image-Echo-GradAfterRefocus_H.py) rather than
#              diffusion-weighted, minimizes how much diffusion contaminates
#              the resulting T2(x) decay, unlike the diffprof.lwg-style
#              pseudo-2D TAU array (which is deliberately diffusion-weighted
#              as well as T2-weighted). Compare a T2(x) map from THIS
#              sequence against a diffusion-weighted array from the existing
#              sequences to separate the two contributions.
#
# Author:      Claude for L. Gordon (DTU)
#
# Created:     11/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.0 (draft)
# Status:      draft -- NOT YET RUN ON HARDWARE. See "BEFORE YOU RUN THIS"
#              below and docs/validation-notes/LWG_CPMG-Image-Echo_H.md.
#              Per this repo's README: no AI-assisted sequence is used for
#              real data collection until a named person has reviewed and
#              validated it on the instrument.
#
# ===============================================================================
# BEFORE YOU RUN THIS -- read in full, especially item (1)
# ===============================================================================
#
# (1) BIGGEST UNVERIFIED ASSUMPTION, CHECK THIS FIRST: this file calls
#     Receiver1() NumEchoes times within a single scan (once per echo),
#     and correspondingly calls TX0.setup_receive(NumEchoes*ReceiverPoints,
#     ...) instead of the usual TX0.setup_receive(ReceiverPoints, ...).
#     EVERY OTHER sequence in this repo (checked via grep across
#     sequences/**/*.py before writing this file) calls Receiver1() exactly
#     ONCE per scan -- there is no precedent here, and no local copy of the
#     firebird framework's own source to check its behaviour against. It is
#     plausible this works exactly as intended (deliver one NumEchoes*points
#     array per scan, echo-major); it is also plausible the framework
#     expects exactly one Receiver1()/setup_receive() pair per scan and this
#     either fails to compile, silently only captures the first or last
#     echo, or does something else entirely. CallBack1D.process_data() below
#     defensively checks the returned array's length against
#     NumEchoes*ReceiverPoints and comms.log()s a loud, specific error
#     (rather than silently mis-reshaping) if it doesn't match -- but that
#     only catches a LENGTH mismatch, not e.g. echoes silently overwriting
#     each other in a same-length buffer. Recommended first test: set
#     NumEchoes=2, NS=1, DS=0, a short-but-safe RD, and inspect (a) whether
#     it compiles at all, (b) the actual length/shape of what
#     process_data() receives, (c) whether the two echoes' data actually
#     look like two DIFFERENT, correctly-spaced spin echoes (plot both
#     halves of the returned array) before trusting this for N>2 or for
#     real data collection. If this assumption is wrong, the fallback is to
#     acquire each echo as its own scan instead (one Receiver1() per scan,
#     TE varied via the existing NBlock/pseudo-2D mechanism already in
#     CallBack1D) -- i.e. structurally back to a diffprof.lwg-style pseudo-
#     2D array, at the cost of losing the single-shot efficiency that makes
#     this a genuine CPMG train rather than a re-excited series.
#
# (2) INTER-ECHO TIMING: derived (not copied) from LWG_1D-Image-Echo_H.py /
#     LWG_1D-Image-Echo-GradAfterRefocus_H.py's validated per-echo formulas
#     -- see "INTER-ECHO TIMING DERIVATION" comment block in run() below for
#     the full derivation. Each echo after the first pads whichever of the
#     gradient/RF branches finishes first (computed explicitly per-run from
#     the actual Parameter values, not assumed) so both branches reach the
#     next refocusing pulse at exactly the same instant -- this is the
#     mechanism that keeps gradient and RF in sync echo-to-echo, and it has
#     NOT been checked against real hardware timing, only worked through
#     algebraically. Echo spacing (2*TAU + this per-run pad) is reported via
#     comms.log() at compile time -- read it and sanity-check it against
#     what you expect before trusting the TE axis.
#
# (3) DUTY CYCLE: a CPMG train applies NumEchoes refocusing pulses AND
#     NumEchoes dephase/readout gradient lobe pairs per scan instead of one
#     -- both RF and gradient duty cycle scale roughly linearly with
#     NumEchoes for a fixed RD. estimate_duty_cycles() below accounts for
#     this (multiplies per-echo on-time by NumEchoes), but the MaxRFDuty/
#     MaxGradDuty thresholds themselves are still the same conservative,
#     NOT vendor-confirmed placeholders as every other sequence in this
#     repo -- take its warnings as a reason to look closer, not as a
#     guarantee either way, especially for large NumEchoes.
#
# Basis:       Per-echo RF (P1 -- TAU -- refocus -- TAU -- acquire) and
#              gradient (idle -- idle -- dephase -- idle -- readout)
#              structure taken from LWG_1D-Image-Echo-GradAfterRefocus_H.py
#              (itself a variant of LWG_1D-Image-Echo_H.py v3.1) --
#              deliberately the diffusion-COMPENSATED variant, not the
#              original diffusion-weighted one, for the reason in Purpose
#              above. Every non-timing helper (apply_gradient, pulse,
#              ninety270, safe_delay, mains_lock_trigger,
#              estimate_duty_cycles, report_probe_gradient, MAXGRAD_TABLE,
#              get_permutation) is carried over unchanged.
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
#   Callback to handle multi-echo (CPMG) one-dimensional sequences. Each
#   scan is expected to return NumEchoes*ReceiverPoints samples (echo-major:
#   echo 1's ReceiverPoints samples, then echo 2's, ...) -- see "BIGGEST
#   UNVERIFIED ASSUMPTION" at the top of this file. self.times is the
#   WITHIN-ECHO time axis (length ReceiverPoints, same for every echo);
#   the per-echo TE values themselves are computed in run() and logged/
#   included in jcamp_meta, not reconstructed here.
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

        self.points_per_echo = Points
        self.num_echoes = int(getattr(Params, "NumEchoes", 1))
        self.expected_len = self.points_per_echo * self.num_echoes

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

        # Defensive check for the multi-Receiver1()-per-scan assumption --
        # see "BIGGEST UNVERIFIED ASSUMPTION" at the top of this file. This
        # only catches a LENGTH mismatch; it cannot detect e.g. echoes
        # silently overwriting each other in a correctly-SIZED buffer.
        if len(scaled_data) != self.expected_len:
            self.comms.log(
                ">>>>> !! CPMG DATA LENGTH MISMATCH on scan {0}: got {1} samples, "
                "expected NumEchoes({2})*ReceiverPoints({3})={4}. The multi-Receiver1()"
                "-per-scan assumption this sequence relies on may not hold on this "
                "framework -- DO NOT trust this data until this is understood. See "
                "item (1) in this file's header. !! <<<<<".format(
                    scan, len(scaled_data), self.num_echoes, self.points_per_echo,
                    self.expected_len))

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


def _refocus_pulse_width(P):
    """Duration of the composite 90-270 refocusing pulse (ninety270()):
    TXEnableTime + 4*P90. Identical convention to LWG_1D-Image-Echo_H.py."""
    return P.P90*4 + 4 + P.TXEnableTime


def _echo_body_times(P, ReceiverFilter, ReceiverTime):
    """(T_rf_body, T_grad_body) = _echo_body_times(...)

    Elapsed time of each branch's per-echo 'body' -- everything AFTER the
    shared pre-refocus idle wait and the refocusing pulse itself, i.e. from
    the refocusing pulse's end up to (but not including) the next echo's
    shared pre-refocus idle wait. See "INTER-ECHO TIMING DERIVATION" in
    run() for the full derivation of both formulas and why they need to be
    explicitly matched (padded) rather than assumed equal.

    T_rf_body   = post-refocus RF wait (B) + ACQU setup/dead-time/acquire
    T_grad_body = dephase lobe + inter-lobe idle + readout lobe
    """
    RefocusPulseWidth = _refocus_pulse_width(P)
    B = P.Tau - RefocusPulseWidth/2.0 + P.PreGrad - ReceiverFilter.dead_time
    T_rf_body = (B + P.TXEnableTime + ReceiverFilter.dead_time + P.Dead1 + ReceiverTime)

    readout_delta = P.GradientOnTime*2 + P.RampTime*2
    T_grad_body = P.Tau + 4.0*P.RampTime + 2.0*P.GradientOnTime + P.GradSettle
    # (equivalent, expanded form, kept here so the two derivations in the
    # header comment and in this function can be checked against each
    # other term-by-term: T_grad_body =
    #   _lobe_elapsed(P, P.GradientOnTime)                   # dephase lobe
    #   + ((P.Tau - P.PreGrad) - _lobe_elapsed(P, P.GradientOnTime))  # idle
    #   + _lobe_elapsed(P, readout_delta)                    # readout lobe
    # which simplifies to the closed form above.)

    return T_rf_body, T_grad_body


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)
    N = int(P.NumEchoes)

    RefocusPulseWidth = _refocus_pulse_width(P)

    # Echo 1 (P1-adjusted pre-refocus wait, same convention as
    # LWG_1D-Image-Echo_H.py -- see run() for the exact formula).
    T_rf_body, T_grad_body = _echo_body_times(P, ReceiverFilter, ReceiverTime)
    pad = abs(T_rf_body - T_grad_body)
    per_echo_body = max(T_rf_body, T_grad_body)

    A1 = P.Tau - (P.P90+3+P.TXEnableTime)/2.0 - RefocusPulseWidth/2.0
    echo1_time = (P.P90 + P.TXEnableTime) + A1 + RefocusPulseWidth + per_echo_body

    # Echoes 2..N: shared pre-refocus wait (Tau - RefocusPulseWidth/2) +
    # refocus + per-echo body (already padded to match on both branches).
    Ai = P.Tau - RefocusPulseWidth/2.0
    subsequent_echo_period = Ai + RefocusPulseWidth + per_echo_body

    t_pulse_gradient_event = echo1_time + (N-1)*subsequent_echo_period

    t_scanTime = P.RecycleDelay + 200.0 + t_pulse_gradient_event
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime


def sequence_description():

    seq_desc = ("Multi-echo (CPMG) 1D MRI / image-echo sequence on the "
                "{H/F} channel. NumEchoes gradient-refocused spin echoes "
                "per scan (90-[TAU-refocus-TAU-acquire]xNumEchoes), each "
                "echo's dephase+readout gradient lobes back-to-back and "
                "diffusion-compensated (see LWG_1D-Image-Echo-"
                "GradAfterRefocus_H.py) so the resulting T2(x) decay across "
                "echoes is minimally confounded by diffusion. DRAFT, not "
                "yet run on hardware -- see file header before use.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,NumEchoes,D70,D73,D71,D75,G1,GradAxis"

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

def ninety270(length90, phase90, phase270, txenabletime):
    """Composite 90(x)-270(y) refocusing pulse. Total duration is
    TXEnableTime + 4*length90."""
    Channel1SetBasePhase(1, phase90)
    Transmit1BlankingOn(1)
    Delay(txenabletime)
    Transmit1(length90)
    Channel1SetBasePhase(1, phase270)
    Transmit1(length90*3)
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


def estimate_duty_cycles(P, rf_on_time_per_echo, grad_on_time_per_echo, comms):
    """Log estimated RF and gradient duty cycles for this scan and warn if
    they exceed MaxRFDuty/MaxGradDuty. UNLIKE the single-echo sequences,
    on-time here is multiplied by NumEchoes -- a CPMG train applies
    NumEchoes refocusing pulses and gradient lobe pairs per scan, not one.

    MaxRFDuty/MaxGradDuty are the SAME conservative, NOT vendor-confirmed
    placeholders used throughout this repo -- confirm real limits with
    Oxford Instruments, especially for large NumEchoes.
    """
    TR = float(P.RecycleDelay)
    N = int(P.NumEchoes)
    rf_duty = (rf_on_time_per_echo * N) / TR
    grad_duty = (grad_on_time_per_echo * N) / TR

    comms.log("Estimated duty cycle (NumEchoes={0}) -- RF: {1:.3%} (limit {2:.1%}), "
              "Gradient [{3}-axis]: {4:.3%} (limit {5:.1%})".format(
              N, rf_duty, P.MaxRFDuty, P.Axis, grad_duty, P.MaxGradDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} (over {1} echoes) exceeds "
                  "MaxRFDuty ({2:.2%}). Consider a longer RD, fewer echoes, or shorter "
                  "pulses/TXEnable, or confirm with Oxford Instruments that this is "
                  "within the transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, N, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} (over {1} echoes) on "
                  "the {2}-axis coil exceeds MaxGradDuty ({3:.2%}). Consider a longer RD, "
                  "fewer echoes, or shorter GradientOnTime, or confirm the rated duty "
                  "cycle for the X-Pulse gradient amplifier/coil with Oxford Instruments."
                  .format(grad_duty, N, P.Axis, P.MaxGradDuty))


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

    ReceiverPoints = Parameter("NP", 1024, ParameterTypes.Int32, "Acquisition Points [per echo]")
    ReceiverAttenuation = Parameter("RA", 34, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "100000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans
    NumScans = Parameter("NS", 4, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 0, ParameterTypes.Int32,"Dummy Scans")

    # CPMG train
    NumEchoes = Parameter("NumEchoes", 8, ParameterTypes.Int32, "Number of CPMG echoes per scan")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 100.0, ParameterTypes.Double, "H/F Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 500000, ParameterTypes.Int32, "Relaxation Delay [s]",
                              RD, min=100000, max=2000000000)

    # Hard pulses
    P90 = Parameter("P90", 9.58, ParameterTypes.Double, "H/F 90&#176; Pulse Width [&#956;s]")
    P180 = Parameter("P180", 19.16, ParameterTypes.Double, "H/F 180&#176; Pulse Width [&#956;s]")
    TXAmplitude = Parameter("RFA0", 0.4, ParameterTypes.Double,
                            "H/F TX Power [0.0&#8230;1.0]", RFA, min=0.0, max=1.0)

    # Gradients -- FPX/FPY/FPZ are the ONLY place per-axis calibration is
    # applied (via GradientMatrix() in run()); G1 is the single logical
    # gradient-strength knob used for both the dephase and readout lobes,
    # on every echo.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet [HFX(default, calibrated)/LOWGAMMA(not yet calibrated)]")

    # Sequence specific -- same meaning/units as LWG_1D-Image-Echo-
    # GradAfterRefocus_H.py, applied identically to every echo.
    GradientOnTime = Parameter("D71", 1000.0, ParameterTypes.Double, "Gradient Duration [&#956;s] -- dephase lobe plateau (readout = 2*this+2*RampTime), per echo")
    RampTime = Parameter("D70", 100.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s] -- both lobes incl. polarity reversal")
    Tau = Parameter("TAU", 10000, ParameterTypes.Int32, "Echo Half-Spacing [&#956;s] -- echo spacing is approx. 2*TAU, see comms.log at run start for the exact value")
    PreGrad = Parameter("D75", 100.0, ParameterTypes.Double, "Pre-Gradient Time [&#956;s]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")

    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rails (conservative, user-adjustable placeholders --
    # NOTE: for this sequence these are checked against NumEchoes*per-echo
    # on-time, see estimate_duty_cycles())
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Phases -- P1 phase-cycled scan-to-scan as usual (PH1). Every
    # refocusing pulse WITHIN a scan uses the SAME PH2/PH3 (standard CPMG:
    # constant refocus phase relative to excitation across the whole
    # train, not stepped per-echo) and every acquisition within a scan
    # uses the same PHRX.
    PH1 = Parameter("PH1", "0,180,90,270,180,0,270,90", ParameterTypes.String, "H/F 90&#176; Pulse Phase")
    PH2 = Parameter("PH2", "0",ParameterTypes.String, "H/F Refocus Pulse Phase (90 part)")
    PH3 = Parameter("PH3", "180",ParameterTypes.String, "H/F Refocus Pulse Phase (270 part)")
    PHRX = Parameter("PHRX", "0,180,270,90,180,0,90,270", ParameterTypes.String, "Acquisition Phase")

def run(comms):

    P = Parameters
    N = int(P.NumEchoes)
    if N < 1:
        raise ValueError("NumEchoes must be >= 1.")

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

    # Per-echo acquisition-window duration -- same formula as every other
    # sequence in this repo, NOT multiplied by N (each echo's own
    # Receiver1() call uses this).
    ReceiverTime = DW*(points+1)

    RefocusPulseWidth = _refocus_pulse_width(P)

    # -------------------------------------------------------------------
    # INTER-ECHO TIMING DERIVATION
    #
    # Goal: refocusing pulse i+1's centre must land exactly 2*TAU after
    # refocusing pulse i's centre, on BOTH the gradient and RF branches,
    # for every echo after the first (echo 1 is special -- see below).
    #
    # Per echo (i>=2), both branches share an IDENTICAL pre-refocus idle
    # wait of (TAU - RefocusPulseWidth/2) -- this is safe to share exactly
    # because neither branch has done anything asymmetric since the
    # previous refocusing pulse's centre once each branch's own "body"
    # (below) is padded to the same total length.
    #
    # After the shared pre-wait and the refocusing pulse itself, the two
    # branches' remaining work per echo ("body") is NOT the same duration:
    #   T_rf_body   = post-refocus wait (Tau - RefocusPulseWidth/2 +
    #                 PreGrad - dead_time) + TXEnableTime + dead_time +
    #                 Dead1 + ReceiverTime
    #               = Tau - RefocusPulseWidth/2 + PreGrad + TXEnableTime
    #                 + Dead1 + ReceiverTime            [dead_time cancels]
    #   T_grad_body = dephase lobe + inter-lobe idle + readout lobe
    #               = Tau + 4*RampTime + 2*GradientOnTime + GradSettle
    #                                          [see _echo_body_times() for
    #                                           the full expansion/algebra]
    #
    # Whichever is shorter gets an explicit extra Delay() (computed here,
    # from real Parameter values, not assumed) equal to the difference, so
    # BOTH branches take exactly max(T_rf_body, T_grad_body) before the
    # next echo's shared pre-wait begins -- this is what keeps gradient
    # and RF synchronised echo-to-echo. See item (2) in the file header:
    # this has been worked through algebraically but NOT checked against
    # real hardware timing.
    #
    # Echo 1 is different: its pre-refocus wait follows the P1 excitation
    # pulse (not a previous refocusing pulse), so it reuses
    # LWG_1D-Image-Echo_H.py's own P1-adjusted formula
    # (Tau - (P90+3+TXEnableTime)/2 - RefocusPulseWidth/2) instead of the
    # shared (Tau - RefocusPulseWidth/2) used for echoes 2..N. Its body
    # (post-refocus wait + ACQU, and dephase+idle+readout) uses the SAME
    # T_rf_body/T_grad_body formulas and the SAME padding as every other
    # echo, so echo 1's own TE differs from later echoes' spacing only by
    # that small, P90-scale pre-wait adjustment -- reported exactly (not
    # assumed uniform) via comms.log below.
    # -------------------------------------------------------------------
    T_rf_body, T_grad_body = _echo_body_times(P, ReceiverFilter, ReceiverTime)
    pad_rf = max(0.0, T_grad_body - T_rf_body)
    pad_grad = max(0.0, T_rf_body - T_grad_body)
    per_echo_body = max(T_rf_body, T_grad_body)

    A1 = P.Tau - (P.P90+3+P.TXEnableTime)/2.0 - RefocusPulseWidth/2.0
    Ai = P.Tau - RefocusPulseWidth/2.0  # shared pre-wait, echoes 2..N

    echo1_TE = (P.P90 + P.TXEnableTime) + A1 + RefocusPulseWidth + per_echo_body
    subsequent_spacing = Ai + RefocusPulseWidth + per_echo_body

    echo_times_us = [echo1_TE + n*subsequent_spacing for n in range(N)]
    comms.log("CPMG echo times (us) relative to P1: {0}".format(
        ["{0:.1f}".format(t) for t in echo_times_us]))
    comms.log("Per-echo branch pad: RF +{0:.2f} us, Gradient +{1:.2f} us "
              "(whichever branch finishes its per-echo body first gets padded "
              "to match the other -- see INTER-ECHO TIMING DERIVATION above)."
              .format(pad_rf, pad_grad))

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
                "NumEchoes={0}".format(N),
                "EchoTimesUs={0}".format(",".join("{0:.3f}".format(t) for t in echo_times_us)),
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
                  "GradientOnTime/RampTime) is shorter than the per-echo "
                  "acquisition window (ReceiverTime = {1:.1f} us, from NP*dwell). "
                  "Part of each echo's FID will be sampled while the read gradient "
                  "is off or ramping -- increase GradientOnTime (D71), reduce NP, "
                  "or use a wider Filter to avoid truncation/leakage."
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

    # ---- Duty-cycle warnings (computed once, outside the scan loop) ------
    # Per-echo on-time (multiplied by N inside estimate_duty_cycles()).
    rf_on_time_per_echo = (P.P90 + P.TXEnableTime) + RefocusPulseWidth
    grad_on_time_per_echo = (_lobe_elapsed(P, P.GradientOnTime)
                              + _lobe_elapsed(P, P.GradientOnTime*2+P.RampTime*2))
    estimate_duty_cycles(P, rf_on_time_per_echo, grad_on_time_per_echo, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver1FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            with parallel:
                with sequential:
                    # ---- GRADIENT BRANCH ----
                    # 90 delay -- matches RF branch's P1 elapsed time
                    Delay((P.P90+P.TXEnableTime+3))
                    for echo_idx in range(N):
                        # Pre-refocus idle (no gradient) -- echo 1 uses the
                        # SAME formula as every other echo here, since this
                        # branch (unlike the RF branch) never transmits P1
                        # itself, so there's no P1-width asymmetry to
                        # correct for on this branch.
                        safe_delay(P.Tau - RefocusPulseWidth/2.0,
                                   "gradient pre-refocus wait, echo {0}".format(echo_idx+1), comms)
                        # Idle through the refocusing pulse itself (this
                        # branch doesn't transmit it, just waits in step)
                        safe_delay(RefocusPulseWidth,
                                   "gradient refocus-width wait, echo {0}".format(echo_idx+1), comms)
                        # Dephase lobe -- NEGATIVE polarity, immediately
                        # after the refocusing pulse (diffusion-compensated
                        # geometry, see LWG_1D-Image-Echo-
                        # GradAfterRefocus_H.py)
                        apply_gradient(P.Axis, P.GradientOnTime, P.RampTime, -P.G1, P.GradSettle, P.PreGrad)
                        # Remainder of the post-refocus period before readout
                        safe_delay((P.Tau - P.PreGrad) - _lobe_elapsed(P, P.GradientOnTime),
                                   "gradient inter-lobe wait, echo {0}".format(echo_idx+1), comms)
                        # Readout lobe -- POSITIVE polarity, spans this
                        # echo's acquisition window
                        apply_gradient(P.Axis, P.GradientOnTime*2+P.RampTime*2, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                        # Resync pad (see INTER-ECHO TIMING DERIVATION) --
                        # skip after the final echo, nothing follows it.
                        if echo_idx < N-1 and pad_grad > 0:
                            safe_delay(pad_grad, "gradient resync pad, echo {0}".format(echo_idx+1), comms)
                with sequential:
                    # ---- RF BRANCH ----
                    # P1 (excitation)
                    pulse(P.P90, ph["PH1"], P.TXEnableTime)
                    for echo_idx in range(N):
                        if echo_idx == 0:
                            # Echo 1: P1-width-adjusted pre-refocus wait
                            # (same formula as LWG_1D-Image-Echo_H.py)
                            safe_delay(P.Tau - (P.P90+3+P.TXEnableTime)/2.0 - RefocusPulseWidth/2.0,
                                       "RF pre-refocus wait, echo 1", comms)
                        else:
                            # Echoes 2..N: shared pre-refocus wait, same as
                            # the gradient branch's for this echo.
                            safe_delay(P.Tau - RefocusPulseWidth/2.0,
                                       "RF pre-refocus wait, echo {0}".format(echo_idx+1), comms)
                        # Refocusing pulse -- SAME PH2/PH3 every echo
                        # (standard CPMG: constant refocus phase relative
                        # to excitation across the whole train)
                        ninety270(P.P90, ph["PH2"], ph["PH3"], P.TXEnableTime)
                        # Post-refocus wait -- ReceiverFilter.dead_time
                        # reserved out of this wait and paid back explicitly
                        # right before Receiver1(), as in every other
                        # sequence in this repo.
                        safe_delay(P.Tau - RefocusPulseWidth/2.0 + 1.*P.PreGrad - ReceiverFilter.dead_time,
                                   "RF post-refocus wait, echo {0}".format(echo_idx+1), comms)
                        # ACQU for this echo
                        Channel1SetBasePhase(P.TXEnableTime,0)
                        Receiver1Phase(P.TXEnableTime, ph["PHRX"])
                        Delay(P.Dead1)
                        Delay(ReceiverFilter.dead_time)
                        Receiver1(ReceiverTime, points)
                        # Resync pad (see INTER-ECHO TIMING DERIVATION) --
                        # skip after the final echo.
                        if echo_idx < N-1 and pad_rf > 0:
                            safe_delay(pad_rf, "RF resync pad, echo {0}".format(echo_idx+1), comms)

        Delay(9.0e4)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        # ---- NumEchoes*ReceiverPoints per scan -- see item (1) in the file
        # header. This is the single biggest unverified assumption in this
        # file: every other sequence in this repo calls
        # TX0.setup_receive(ReceiverPoints, ...) because every other
        # sequence calls Receiver1() exactly once per scan. This one calls
        # it N times per scan above, so the buffer is sized for all N
        # echoes' worth of samples -- VALIDATE this actually delivers
        # NumEchoes*ReceiverPoints samples, echo-major, before trusting any
        # data from this sequence (CallBack1D.process_data() above will log
        # a loud error if the LENGTH doesn't match, but cannot detect
        # echoes silently overwriting each other in a correctly-sized
        # buffer).
        TX0.setup_receive(points*N,32,scan=seqScans,circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX0.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))
        comms.log("NumEchoes %s, Encode Time %s" % (N, P.GradientOnTime))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 11/08/26 (v1.0, draft) - Initial version. Multi-echo (CPMG)
#      adaptation of LWG_1D-Image-Echo-GradAfterRefocus_H.py (itself a
#      diffusion-compensated variant of LWG_1D-Image-Echo_H.py v3.1),
#      written to answer whether the effective-T2 gradient seen near a
#      liquid-liquid interface in octanol-water imaging data is really T2
#      or is diffusion contamination from the single-echo sequences'
#      inherent (if small) diffusion weighting -- see chat discussion,
#      2026-08-11. NOT YET RUN ON HARDWARE -- see "BEFORE YOU RUN THIS" at
#      the top of this file, especially item (1) (multi-Receiver1()-per-
#      scan is unprecedented in this repo and unverified against the real
#      framework behaviour) and item (2) (inter-echo branch-sync padding
#      derived algebraically, not checked against real timing).
#
# -----------------------------------------------------------------------------
