#-------------------------------------------------------------------------------
# Name:        LWG_1D-Image-Echo-GradAfterRefocus_H.py
# Purpose:     Variant of LWG_1D-Image-Echo_H.py (v3.1) that moves ALL
#              gradient activity to the SECOND TAU period (after the
#              refocusing pulse) instead of splitting it: a dephase lobe
#              before refocusing (exploiting the 180's sign flip) + a
#              readout lobe after. Written as an experimental A/B test
#              against the original, NOT a validated replacement -- see
#              "IMPLICATIONS" below and the changelog before running.
#
# Author:      Claude for L. Gordon (DTU), branched from
#              LWG_1D-Image-Echo_H.py v3.1 (CM/LWG base, Claude revisions)
#
# Created:     06/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.3 (experimental)
#
# WHAT CHANGED VS v3.1:
#   v3.1 (original):  90 -- [dephase gradient, +G1] -- TAU -- 180 -- TAU --
#                      [readout gradient, +G1] -- echo/acquisition
#                      The dephase lobe sits BEFORE the refocusing pulse.
#                      Same polarity (+G1) as the readout lobe is used for
#                      BOTH, because the 180 pulse itself inverts the
#                      accumulated phase -- so a same-sign pre-180 lobe acts
#                      as a NEGATIVE (dephase) lobe once you get to the
#                      post-180 frame. This is the standard trick used to
#                      align a gradient echo (k=0) with a spin echo.
#
#   This file:         90 -- TAU (no gradient) -- 180 -- [dephase gradient,
#                      -G1] -- [readout gradient, +G1] -- echo/acquisition
#                      The dephase lobe is moved to immediately AFTER the
#                      refocusing pulse, with genuinely OPPOSITE polarity
#                      (there is no 180 left to exploit for a sign flip --
#                      both lobes now happen on the same side of it). The
#                      readout lobe itself -- amplitude, duration, and its
#                      timing relative to the acquisition window -- is
#                      UNCHANGED from v3.1, so the already-validated
#                      "echo centred in the acquisition window" timing is
#                      untouched; only the source of the -kmax starting
#                      offset moves from before-180 to after-180.
#
# IMPLICATIONS (see also chat discussion -- summarising here for anyone
# reading this file cold):
#
#   T2 / echo time (TE): UNCHANGED. TE = 2*TAU is set entirely by the RF
#   branch (90 -- TAU -- 180 -- TAU -- acquisition), which this file does
#   not touch at all. Moving the gradient's position within the TAU-TAU
#   window does not change when the spin echo happens or how much T2 decay
#   has occurred by the time you acquire it.
#
#   Diffusion weighting: THIS is what actually changes, and is the main
#   reason to try this variant. In v3.1, the dephase lobe (before 180) and
#   the readout lobe (after 180) are separated by roughly TAU and straddle
#   the refocusing pulse with the SAME nominal polarity -- which, after the
#   180 flips the sign of the first one, is exactly the Stejskal-Tanner
#   (PGSE) diffusion-weighting geometry: two gradient pulses of effective
#   opposite sign, separated by a diffusion time Delta ~ TAU, refocused by
#   a 180 in between. That geometry inevitably imparts a small but real
#   b-value on top of the intended spatial encoding -- so v3.1's profile
#   intensity is, strictly, a little diffusion-weighted as well as
#   T2-weighted, and a sample with a higher self-diffusion coefficient
#   (e.g. free water vs. something more restricted) will show somewhat
#   more signal loss than its spin density alone would predict.
#
#   In this file, BOTH gradient lobes sit on the same (post-180) side of
#   the refocusing pulse, back-to-back with no refocusing pulse between
#   them. That is a diffusion-COMPENSATED (bipolar) geometry, not a
#   diffusion-weighting one -- the separation between the two lobes' centres
#   is now on the order of the lobes' own duration (microseconds-to-low-
#   milliseconds), not TAU, so the effective b-value collapses to something
#   much smaller (b scales with roughly the square of the lobe separation).
#   The trade-off: this profile is "purer" spin-density/T2-weighted data,
#   less confounded by the sample's diffusion coefficient -- useful if
#   you're comparing regions/samples with different D and don't want that
#   difference bleeding into apparent intensity.
#
#   Practical costs of the swap:
#     - The post-180 gradient branch now has to fit BOTH lobes into the
#       second TAU instead of one, so there is less slack before TAU is
#       too short for your GradientOnTime/RampTime/GradSettle combination
#       -- watch for the safe_delay() FATAL TIMING ERROR on this branch
#       and lengthen TAU if you hit it (this compounds with, doesn't
#       replace, the "increase TAU to move the echo away from pulse
#       bleed" advice from earlier).
#     - Two back-to-back OPPOSITE-polarity lobes (dephase then readout)
#       means the gradient amplifier has to reverse polarity quickly
#       instead of ramping to zero and idling -- if your amplifier's
#       eddy-current settling is polarity-direction-dependent, you may
#       need MORE settle time between these two lobes than GradSettle
#       currently gives it (see the GradSettle re-use note in run()
#       below); this is exactly the kind of thing that could produce a
#       new artefact of its own, so validate this file's profile shape
#       (symmetric, positive, centred) before trusting it quantitatively,
#       the same way the v3.1 rollout was validated.
#     - Estimated gradient duty cycle is unchanged in total on-time, but
#       is now concentrated entirely in the second TAU rather than spread
#       across both -- check the MaxGradDuty warning if you shorten RD.
#
# Basis:       Directly derived from LWG_1D-Image-Echo_H.py v3.1 -- only
#              the gradient branch of run() differs; the RF branch,
#              Parameters block, helper functions and CallBack1D are
#              carried over unchanged (this is deliberate: keeping the RF/
#              acquisition timing byte-for-byte identical is what makes
#              the T2/TE-unchanged claim above actually true, rather than
#              just asserted).
#
# NOTE ON FILE ENCODING: the X-Pulse pulse-sequence compiler requires
# Windows CRLF line endings (see LWG_1D-Image-Echo_H.py's own changelog
# item 10g) -- this file must be saved/kept as CRLF.
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
#   Callback to handle one-dimensional sequences
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
            self.comms.log("seqTime LWG_1D-Image-Echo-GradAfterRefocus_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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

def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)

    RefocusPulseWidth = P.P90*4 + 4 + P.TXEnableTime

    def lobe_elapsed(delta):
        return P.PreGrad + 2.0*P.RampTime + delta + P.GradSettle

    # The pulse/gradient event in run() is TWO PARALLEL branches (gradient
    # timing vs RF timing) inside a single 'with parallel:' block -- its
    # real elapsed time is the MAX of the two branches, not their sum (see
    # xpulse-pulse-programmes skill sec.4). Each branch total below is
    # derived by summing every Delay()/duration argument in that exact
    # branch of run(), term for term -- this variant's gradient branch
    # differs from LWG_1D-Image-Echo_H.py's (dephase lobe moved after the
    # refocusing pulse), but its own lobe_elapsed(GOT) cancels out exactly
    # the same way the encode-lobe term did there.
    T_grad = (P.P90 + P.TXEnableTime + 3)
    T_grad += (P.Tau - RefocusPulseWidth/2.0)
    T_grad += RefocusPulseWidth
    T_grad += lobe_elapsed(P.GradientOnTime)
    T_grad += ((P.Tau - P.PreGrad) - lobe_elapsed(P.GradientOnTime))
    T_grad += lobe_elapsed(P.GradientOnTime*2 + P.RampTime*2)

    T_rf = (1 + P.TXEnableTime + P.P90)
    T_rf += (P.Tau - (P.P90+3+P.TXEnableTime)/2 - RefocusPulseWidth/2.0)
    T_rf += (1 + P.TXEnableTime + P.P90 + 1 + P.P90*3)
    T_rf += (P.Tau - RefocusPulseWidth/2.0 + P.PreGrad - ReceiverFilter.dead_time)
    T_rf += 2*P.TXEnableTime
    T_rf += P.Dead1
    T_rf += ReceiverFilter.dead_time
    T_rf += ReceiverTime

    # RD was previously OMITTED here entirely -- by far the dominant term
    # -- and the whole gradient/RF-parallel structure was collapsed into a
    # bare P.P90*5 approximation that accounted for neither gradients nor
    # acquisition. Receiver1FilterFlush's 200us is the only fixed cost
    # outside the parallel block.
    t_scanTime = P.RecycleDelay + 200.0 + max(T_grad, T_rf)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("1D MRI / image-echo (frequency-encoded profile) sequence on "
                "the {H/F} channel -- GradAfterRefocus variant.\nSpin-echo "
                "(90-TAU-[90-270 composite 180]-TAU-echo) with BOTH the "
                "dephase and readout gradient lobes placed after the "
                "refocusing pulse (bipolar, back-to-back), instead of "
                "splitting them before/after as in LWG_1D-Image-Echo_H.py. "
                "TE (T2 weighting) is unchanged; diffusion weighting is "
                "greatly reduced since the two lobes no longer straddle the "
                "refocusing pulse with a ~TAU separation. Gradient axis "
                "selectable via GradAxis (x/y/z).")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,D70,D73,D71,D75,G1,GradAxis,Probe,Dead1"

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
    in run() below, NOT here. gradient_value's SIGN sets the lobe's
    polarity -- this file calls it with -P.G1 for the post-refocus dephase
    lobe and +P.G1 for the readout lobe, unlike v3.1 which only ever uses
    +P.G1 (relying on the refocusing pulse to flip the dephase lobe's
    effective sign instead)."""
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
    TXEnableTime + 4*length90 -- every Delay() around this call in run()
    accounts for that true duration via RefocusPulseWidth."""
    Channel1SetBasePhase(1, phase90)
    Transmit1BlankingOn(1)
    Delay(txenabletime)
    Transmit1(length90)
    Channel1SetBasePhase(1, phase270)
    Transmit1(length90*3)
    Transmit1BlankingOff(1)


def safe_delay(value, label, comms):
    """Delay() wrapper for every timing-formula-derived (as opposed to
    fixed) delay in the sequence. If TAU/GradientOnTime/RampTime/PreGrad
    don't leave enough room, the raw formula would go negative -- which can
    silently desynchronise the read gradient from the acquisition window
    ('leaking' signal/artefacts outside the intended echo) instead of
    failing loudly. This turns that failure mode into a specific,
    immediate error instead.

    NOTE for this variant: the second-TAU wait now has to make room for
    the dephase lobe as well as the readout lobe (v3.1 only had the
    readout lobe there), so this is more likely to go negative than in
    v3.1 for the same TAU/GradientOnTime -- if you hit this on the
    'gradient second-TAU wait (after dephase)' label, increase TAU."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU and/or PreGrad, or decrease "
               "GradientOnTime/RampTime/GradSettle, then retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger before the first pulse event.
    OFF by default -- see LWG_1D-Image-Echo_H.py's docstring for the full
    rationale (unchanged here)."""
    if int(P.UseMainsLock) == 0:
        return
    triggers = {1: ExternalTrigger1, 2: ExternalTrigger2, 3: ExternalTrigger3}
    ch = int(P.MainsLockChannel)
    if ch not in triggers:
        comms.log("WARNING: MainsLockChannel={0} is not 1, 2, or 3 -- mains "
                  "lock trigger skipped.".format(ch))
        return
    triggers[ch]()


def estimate_duty_cycles(P, rf_on_time, grad_on_time, comms):
    """Log estimated RF and gradient duty cycles for this scan and warn if
    they exceed the (user-adjustable) MaxRFDuty / MaxGradDuty parameters.
    Same conservative-placeholder caveat as v3.1 applies."""
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR
    grad_duty = grad_on_time / TR

    comms.log("Estimated duty cycle -- RF: {0:.3%} (limit {1:.1%}), "
              "Gradient [{2}-axis]: {3:.3%} (limit {4:.1%})".format(
              rf_duty, P.MaxRFDuty, P.Axis, grad_duty, P.MaxGradDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD or shorter pulses/TXEnable, "
                  "or confirm with Oxford Instruments that this is within the "
                  "transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} on the "
                  "{1}-axis coil exceeds MaxGradDuty ({2:.2%}). Consider a "
                  "longer RD or shorter GradientOnTime, or confirm the rated "
                  "duty cycle for the X-Pulse gradient amplifier/coil with "
                  "Oxford Instruments -- this limit is a conservative "
                  "placeholder, not a vendor spec."
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
    Sequence = Parameter("Sequence", "LWG_1D-Image-Echo-GradAfterRefocus_H", ParameterTypes.String, "Sequence Name")

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
    # gradient-strength knob used for both the dephase and readout lobes.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet [HFX(default, calibrated)/Low Gamma(not yet calibrated)]")
    ProbeList = Parameter("ProbeList", "HFX,Low Gamma", ParameterTypes.String, "Probe Options")

    # Sequence specific
    GradientOnTime = Parameter("D71", 1000.0, ParameterTypes.Double, "Gradient Duration [&#956;s] -- dephase lobe plateau (readout = 2*this+2*RampTime)")
    RampTime = Parameter("D70", 100.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s] -- used for both lobes incl. polarity reversal")
    Tau = Parameter("TAU", 10000, ParameterTypes.Int32, "Echo Time [&#956;s]")
    PreGrad = Parameter("D75", 100.0, ParameterTypes.Double, "Pre-Gradient Time [&#956;s]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")

    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rails (conservative, user-adjustable placeholders)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Phases
    PH1 = Parameter("PH1", "0,180,90,270,180,0,270,90", ParameterTypes.String, "H/F 180&#176; Pulse Phase")
    PH2 = Parameter("PH2", "0",ParameterTypes.String, "H/F 180&#176; Pulse Phase")
    PH3 = Parameter("PH3", "180",ParameterTypes.String, "H/F 180&#176; Pulse Phase")
    PHRX = Parameter("PHRX", "0,180,270,90,180,0,90,270", ParameterTypes.String, "Acquisition Phase")

def run(comms):

    P = Parameters

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

    # Total acquisition-window duration -- Receiver1()'s 'duration' arg.
    ReceiverTime = DW*(points+1)

    times = np.arange(0, points*DW, DW) / 1.0e6

    # Duration of the refocusing element (composite 90-270, via
    # ninety270()): TXEnableTime + 4*P90. Identical to v3.1.
    RefocusPulseWidth = P.P90*4 + 4 + P.TXEnableTime

    # Elapsed time of a single apply_gradient() lobe call, given its own
    # plateau duration 'delta' -- used to keep the gradient-branch timing
    # formulas below self-documenting instead of repeating the same sum.
    def lobe_elapsed(delta):
        return P.PreGrad + 2.0*P.RampTime + delta + P.GradSettle

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
                "jc_ns=1"
                    ]
        return jc_data

    comms.log("Sequence Time: {0}".format(time_calculation(Parameters)))
    seqAcqu = CallBack1D(Parameters, comms, ReceiverFilter, jcamp_meta())

    def RecvCallback(acqData, scan):
        seqAcqu.process_data(scan, acqData)

    # ---- Read-gradient coverage check (warning only, not a hard stop) ----
    # Identical check to v3.1: the readout lobe's own plateau (unchanged
    # here) should span at least the acquisition window.
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

    # ---- Duty-cycle warnings (computed once, outside the scan loop) ------
    # Total gradient on-time is the same two lobes as v3.1 (one of duration
    # GradientOnTime, one of duration GradientOnTime*2+RampTime*2) -- just
    # both now on the post-refocus side.
    rf_on_time = (P.P90 + P.TXEnableTime) + (RefocusPulseWidth)
    grad_on_time = lobe_elapsed(P.GradientOnTime) + lobe_elapsed(P.GradientOnTime*2+P.RampTime*2)
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

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
                    # 90 delay -- mirrors the RF branch's P1 elapsed time,
                    # same as v3.1.
                    Delay((P.P90+P.TXEnableTime+3))
                    # No gradient during the first TAU in this variant --
                    # wait out the FULL first-TAU period (v3.1 spent part
                    # of it on the encode lobe; here that time is simply
                    # idle, and the lobe reappears after the refocusing
                    # pulse below instead).
                    safe_delay(P.Tau - RefocusPulseWidth/2.0,
                               "gradient first-TAU wait (grad-after-refocus variant)", comms)
                    # Wait out the refocusing pulse itself (this branch
                    # doesn't transmit it, just needs to stay in step).
                    safe_delay(RefocusPulseWidth,
                               "gradient refocus-width wait", comms)
                    # Dephase gradient -- NOW placed immediately after the
                    # refocusing pulse, with explicit NEGATIVE polarity
                    # (there's no more 180 ahead of it to flip the sign
                    # for us). Same GradientOnTime/RampTime/GradSettle as
                    # v3.1's pre-180 encode lobe, so the -kmax offset it
                    # produces is identical in magnitude.
                    apply_gradient(P.Axis, P.GradientOnTime, P.RampTime, -P.G1, P.GradSettle, P.PreGrad)
                    # Wait out the remainder of the second TAU -- same
                    # total (Tau - PreGrad) as v3.1's "second-TAU wait",
                    # minus the time the dephase lobe above just spent, so
                    # the readout lobe below still starts at the exact
                    # same moment relative to TAU as it did in v3.1 (i.e.
                    # the already-validated acquisition-window centring is
                    # untouched).
                    safe_delay((P.Tau - P.PreGrad) - lobe_elapsed(P.GradientOnTime),
                               "gradient second-TAU wait (after dephase)", comms)
                    # Read Gradient -- unchanged from v3.1: same duration,
                    # same +G1 polarity, same position relative to TAU.
                    apply_gradient(P.Axis, P.GradientOnTime*2+P.RampTime*2, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                with sequential:
                    # P1
                    pulse(P.P90, ph["PH1"], P.TXEnableTime)
                    # TAU
                    safe_delay(P.Tau - (P.P90+3+P.TXEnableTime)/2 - RefocusPulseWidth/2.0,
                               "RF first-TAU wait", comms) # first tau delay
                    # P2
                    ninety270(P.P90, ph["PH2"], ph["PH3"], P.TXEnableTime)
                    # TAU -- ReceiverFilter.dead_time is reserved out of this
                    # wait and paid back explicitly (as its own Delay, right
                    # before Receiver1 below), exactly as in v3.1. This RF
                    # branch is BYTE-FOR-BYTE IDENTICAL to v3.1's -- TE and
                    # T2 weighting are therefore unchanged by this variant.
                    safe_delay(P.Tau - RefocusPulseWidth/2.0 + 1.*P.PreGrad - ReceiverFilter.dead_time,
                               "RF second-TAU wait", comms)
                    # ACQU
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
        comms.log("Dephase/Encode Time %s" % (P.GradientOnTime))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 06/08/26 - Branched from LWG_1D-Image-Echo_H.py v3.1 as an
#      experimental A/B test. Only the gradient branch of run() changed:
#      the dephase lobe moved from before the refocusing pulse (same
#      polarity as the readout lobe, relying on the 180 to flip its
#      effective sign) to immediately after it (explicit opposite
#      polarity, back-to-back with the readout lobe). RF branch, Parameters,
#      helper functions, and CallBack1D are otherwise unchanged from v3.1.
#      See the file header for the full T2-vs-diffusion-weighting rationale
#      and the practical timing/duty-cycle trade-offs of the swap.
#      NOT YET VALIDATED ON HARDWARE -- before trusting results
#      quantitatively, confirm (as was done for v3.1) that the profile
#      comes out symmetric, positive after phasing, and centred in the
#      acquisition window, and watch for a new artefact at the dephase/
#      readout polarity-reversal boundary if GradSettle turns out to be
#      too short for your amplifier's eddy-current behaviour there.
# 2. Claude - 10/08/26 - FIXED time_calculation(): the reported 'Sequence
#    Time' was a bare P.P90*5 + P.Tau*2 + points*DW guess that omitted
#    RecycleDelay (RD, by far the dominant term) AND every gradient timing
#    parameter entirely (found while auditing every pp file's
#    time_calculation() after you reported the reported time never
#    reflecting real parameters). Rewrote to derive the actual per-scan
#    time from run()'s real structure: MAX of the gradient/RF parallel
#    branches (RF branch identical to LWG_1D-Image-Echo_H.py's; gradient
#    branch reflects this variant's post-refocus dephase-lobe placement),
#    each summed term for term, + RD + Receiver1FilterFlush. Verified via
#    the mock harness (2091054.96us at default Parameters -- differs from
#    LWG_1D-Image-Echo_H.py's 2091462.96us by exactly the expected 408us,
#    i.e. 102us/scan x 4 scans, matching the gradient-branch structural
#    difference between the two variants).
#
# 3. Claude - 18/08/26 - Added a ProbeList Parameter ("HFX,Low Gamma")
#    paired with Probe (same convention as an existing GradAxis/
#    GradAxisList pairing elsewhere in this repo, e.g.
#    LWG_CPMG-Image-Echo_H.py) -- this is how this framework exposes a
#    real dropdown (a Parameter plus a comma-separated "...List"
#    companion), not a distinct enum ParameterType. Probe's description
#    now reads "Low Gamma" (matching your wording) instead of "LOWGAMMA"
#    -- report_probe_gradient()'s existing normalization (strip spaces/
#    hyphens, uppercase) already maps that to the unchanged MAXGRAD_TABLE
#    key "LOWGAMMA", so no table changes were needed. Purely a Parameter-
#    panel change -- no effect on timing/hardware calls.
#
# 4. Claude - 18/08/26 - Added Probe and Dead1 to sequence_basic() at
#    your request, so both already-existing Parameters show up in
#    SpinFlow's quick "basic parameters" panel rather than requiring the
#    full parameter list to find/set them. No change to the Parameters
#    themselves or to timing/hardware calls.
#
# -----------------------------------------------------------------------------
