#-------------------------------------------------------------------------------
# Name:        LWG_1D-Image-UTE_X.py
# Purpose:     One-dimensional MRI / short-echo-time ('quasi-UTE') profile
#              sequence for short-T2/T2* species (e.g. quadrupolar nuclei
#              such as 23Na, or semi-solids) on the Oxford Instruments
#              X-Pulse Broadband Benchtop NMR Spectrometer (1H/19F channel
#              in this _H version). Axis-selectable (x/y/z).
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     04/08/2026
# Revised:     06/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.2
#
# Design notes -- READ THIS, this sequence is structurally different from
# LWG_1D-Image-Echo_H.py, not just a parameter change:
#
#  - This is NOT a spin echo. There is no refocusing pulse. The read
#    gradient is ramped up and settled BEFORE the excitation pulse, the
#    pulse fires with the gradient already flat, and acquisition starts as
#    soon as hardware-fixed dead times (TXEnableTime/blanking switching,
#    Dead1 probe ring-down, ReceiverFilter.dead_time filter settling) allow
#    -- there is no TAU to wait out, so TE here is essentially the sum of
#    those fixed, mostly-hardware-determined delays. That's what makes it
#    'ultra-short' relative to the ms-scale TE of the spin-echo profile
#    sequences.
#  - This acquires ONE-SIDED k-space (from k=0, at the pulse, out to +kmax,
#    not -kmax to +kmax like the spin-echo profile sequences). This is
#    standard for UTE/ZTE-style acquisition: you cannot dephase to -kmax
#    before the pulse without adding exactly the delay you are trying to
#    eliminate. In practice this is equivalent to acquiring an FID under a
#    read gradient and treating its magnitude Fourier transform as the
#    projection -- normal, accepted practice for magnitude-mode short-T2
#    imaging, but if you need a fully symmetric/phased profile you will
#    need to acquire twice with the gradient sign flipped and combine, or
#    add a genuine ramp-sampling/half-echo reconstruction step in
#    post-processing. Not attempted here.
#  - This spectrometer is a general-purpose benchtop system, not
#    UTE-optimised hardware -- the achievable minimum TE is set by
#    TXEnableTime, P90/2, a couple of microseconds of fixed instruction
#    overhead, Dead1 (probe ring-down), and ReceiverFilter.dead_time
#    (which depends on the Filter/bandwidth you choose -- wider Filter =
#    shorter dead_time = shorter TE). run() computes and comms.logs the
#    resulting 'MinimumTE' every time you run this, so you can see exactly
#    what you're getting and what to shorten (shorter P90, shorter
#    TXEnableTime if your hardware allows it, and/or a wider Filter) if you
#    need it shorter still. This will likely be tens to a couple hundred
#    microseconds on this hardware, not the sub-10us some dedicated
#    UTE/solid-state probes achieve -- still a large improvement over a
#    ms-scale spin-echo TE for many short-T2/T2* species.
#  - No refocusing pulse also means no phase cycling is doing anything
#    useful beyond simple CYCLOPS-style receiver cycling -- kept simple
#    (2-step) accordingly.
#
#
# X-CHANNEL CALIBRATION -- IMPORTANT: this is a mechanical H->X port (all
# Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1). SF
# now defaults to 15.01 MHz as a placeholder -- P90/TXAmplitude/RA below are
# still the 1H-calibrated values from the _H file and are almost certainly
# WRONG for whatever nucleus you actually put on the X channel. This
# matters more than usual here: for a short-T2 quadrupolar nucleus like
# 23Na, get the hard-pulse calibration (P90, TXAmplitude) and probe tuning
# right first, since the whole point of this sequence is to not waste any
# of a signal that may already be decaying within tens to hundreds of us.
# Recalibrate SF, P90, TXAmplitude (RFA0) and RA for your target nucleus
# before running.
#
# Changes/Modifications: At end of file.
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
            self.comms.log("seqTime LWG_1D-Image-UTE_X Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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

    # Matches run()'s actual per-scan timing exactly (previous formula
    # already included RD -- unlike most other pp files in this repo, see
    # changelog for the wider audit -- but omitted PreGrad/RampTime/
    # GradSettle, a small ~0.1% underestimate given RD already dominates).
    # The pulse/gradient event is a 'with parallel:' block whose two
    # branches (gradient hold vs RF+acquisition) are constructed in run()
    # to have matching total elapsed time via Hold/TrailingPad -- the
    # gradient branch is used here since it's very slightly the longer of
    # the two by construction (its own PreGrad/2*RampTime/GradSettle
    # overhead isn't offset by anything on the RF side). RD +
    # Receiver2FilterFlush's 200us are the fixed costs outside the
    # parallel block.
    PulseElapsed = P.P90 + P.TXEnableTime + 3
    PostPulseOverhead = 2
    Hold = PulseElapsed + PostPulseOverhead + P.Dead1 + ReceiverFilter.dead_time + ReceiverTime
    GradBranch = P.PreGrad + Hold + 2*P.RampTime + P.GradSettle + 2

    t_scanTime = P.RecycleDelay + 200.0 + GradBranch
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("1D MRI short-echo-time ('quasi-UTE') profile sequence on "
                "the {X} channel, for short-T2/T2* species.\nNo "
                "refocusing pulse -- the read gradient is on and settled "
                "BEFORE excitation, and acquisition begins as soon as fixed "
                "hardware dead times allow. Acquires one-sided k-space "
                "(k=0 at the pulse, out to +kmax). Gradient axis selectable "
                "via GradAxis (x/y/z).")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,D70,D73,D71,D75,G1,GradAxis,P90"

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
    """Ramp a gradient axis up to gradient_value, hold flat for delta, ramp
    back to zero, then settle. gradient_value is the final logical
    (normalised, -1..1) amplitude -- axis calibration (FPX/FPY/FPZ) is
    applied exactly once, by GradientMatrix() in run() below, NOT here."""
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
    """Apply a short, non-selective (broadband) hard pulse."""
    Channel2SetBasePhase(1, phase)
    Transmit2BlankingOn(1)
    Delay(txenabletime)
    Transmit2(length)
    Transmit2BlankingOff(1)


def safe_delay(value, label, comms):
    """Delay() wrapper with a clear, specific error if a computed delay
    would be negative, instead of a silent/confusing failure."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Check RampTime/GradSettle/PreGrad/Dead1, then "
               "retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger (ExternalTrigger[n]()) before the
    first pulse event, per Pulse Sequence Programming User Manual 01-U-049
    section 3.7.19. The manual's trigger-channel table only lists MQC+
    (mains lock = ExternalTrigger2()) and MQR (mains lock =
    ExternalTrigger3()) -- no X-Pulse row is given, so the channel is a
    Parameter for you to confirm/correct rather than hard-coded.

    OFF (UseMainsLock=0) by default. Reasoning: a mains-lock trigger waits
    for the next AC zero-crossing before firing, adding a non-deterministic
    delay (0 up to one mains half-cycle -- up to ~10 ms @ 50 Hz / ~8.3 ms @
    60 Hz) before the first event of EVERY scan. This is especially
    unwelcome here: the whole point of this sequence is a short, TIGHT,
    reproducible TE, and a variable pre-pulse delay would undermine that
    directly. Set UseMainsLock=1 (and confirm MainsLockChannel) only if you
    determine your X-Pulse needs mains-synchronous triggering.
    """
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

    IMPORTANT: MaxRFDuty and MaxGradDuty are conservative, editable
    placeholders (5% / 10%), NOT a vendor-confirmed rating -- no public duty
    cycle specification for the X-Pulse 60 MHz RF/gradient amplifiers was
    available when this was written. Confirm the real limits with Oxford
    Instruments. Note the gradient is held on for the FULL scan (pulse +
    dead time + acquisition), not just a short lobe, so gradient duty cycle
    is worth watching closely here, especially with short RD/high NS.
    """
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR
    grad_duty = grad_on_time / TR

    comms.log("Estimated duty cycle -- RF: {0:.3%} (limit {1:.1%}), "
              "Gradient [{2}-axis]: {3:.3%} (limit {4:.1%})".format(
              rf_duty, P.MaxRFDuty, P.Axis, grad_duty, P.MaxGradDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD or shorter P90, or "
                  "confirm with Oxford Instruments that this is within the "
                  "transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} on the "
                  "{1}-axis coil exceeds MaxGradDuty ({2:.2%}). The gradient "
                  "is held on for the whole pulse+acquisition window in this "
                  "sequence, so this is easy to hit with short RD/large NP -- "
                  "consider a longer RD, fewer NP, or a lower G1, or confirm "
                  "the rated duty cycle for the X-Pulse gradient amplifier/"
                  "coil with Oxford Instruments -- this limit is a "
                  "conservative placeholder, not a vendor spec."
                  .format(grad_duty, P.Axis, P.MaxGradDuty))


# Per-axis max gradient strength (G/cm) for probes that can be fitted to
# this magnet. There is no way to detect which probe is mounted from
# software, so this is a MANUAL selector (Probe Parameter below) -- the
# operator must set it to match what is actually on the magnet. KEEP THIS
# TABLE IN SYNC with leonmr/xpulse_imaging.py's own GRADIENT_CALIBRATION
# dict (duplicated here rather than imported, since pulse programmes can't
# import external Python modules -- same house convention used for
# generate_shape()/safe_delay()/mains_lock_trigger() across this family).
MAXGRAD_TABLE = {
    "HFX": {"x": 11.879, "y": 11.978, "z": 57.915},   # G/cm, measured/averaged 06/08/2026
    "LOWGAMMA": {"x": None, "y": None, "z": None},     # NOT YET CALIBRATED
}

def report_probe_gradient(P, comms, axis=None):
    """Log the max gradient strength (G/cm) for the currently-selected Probe
    and gradient axis. Since Probe and GradAxis are both Parameters, this
    also gets auto-recorded in the resulting JCAMP file's SpinFlow
    parameter block -- letting leonmr/xpulse_imaging.py convert Hz to mm
    later without the operator having to remember/re-enter which probe was
    fitted for a given experiment. Warns (does not fail) if the fitted
    probe's axis has no calibration yet."""
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
    Sequence = Parameter("Sequence", "LWG_1D-Image-UTE_X", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 15.01, ParameterTypes.Double, "X Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "X Freq Offset [Hz]")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "X TX Freq Offset [ppm]")

    ReceiverPoints = Parameter("NP", 1024, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 34, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    # Wide Filter = short ReceiverFilter.dead_time = shorter achievable TE.
    Filter = Parameter("Filter", "100000", ParameterTypes.String, "f2 Spectral Window [Hz] -- wider = shorter dead_time = shorter TE")

    # Scans
    NumScans = Parameter("NS", 4, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 0, ParameterTypes.Int32,"Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 20.0, ParameterTypes.Double, "X Probe Ringdown Time [&#956;s] -- part of TE floor")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s] -- part of TE floor")
    RecycleDelay = Parameter("RD", 500000, ParameterTypes.Int32, "Relaxation Delay [s]",
                              RD, min=100000, max=2000000000)

    # Hard pulse -- short, broadband, non-selective. A smaller flip angle
    # (shorter P90 and/or lower TXAmplitude) both shortens TE slightly and
    # allows faster repetition; tune to your T1/SNR needs.
    P90 = Parameter("P90", 9.58, ParameterTypes.Double, "X 90&#176; (or smaller) Pulse Width [&#956;s] -- TE includes P90/2")
    TXAmplitude = Parameter("RFA0", 0.4, ParameterTypes.Double,
                            "X TX Power [0.0&#8230;1.0]", RFA, min=0.0, max=1.0)

    # Gradients -- FPX/FPY/FPZ are the ONLY place per-axis calibration is
    # applied (via GradientMatrix() in run()); G1 is the logical gradient-
    # strength knob. Unlike the spin-echo profile sequences, there is only
    # ONE gradient event here: it is ramped up and settled before the
    # pulse, held flat through the pulse + dead time + acquisition, then
    # ramped down.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    # Probe fitted to the magnet -- no automatic detection is possible, so
    # this must be set MANUALLY to match what's actually mounted. Used only
    # for logging/downstream Hz-to-mm conversion (see report_probe_gradient()
    # above and leonmr/xpulse_imaging.py) -- has no effect on this pp's own
    # timing/hardware calls.
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet [HFX(default, calibrated)/LOWGAMMA(not yet calibrated)]")

    RampTime = Parameter("D70", 100.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s]")
    PreGrad = Parameter("D75", 100.0, ParameterTypes.Double, "Pre-Gradient Time [&#956;s] -- before pulse, does not add to TE")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")

    # Mains-lock trigger -- OFF by default (see mains_lock_trigger()
    # docstring: especially undesirable here, since it would add a variable
    # pre-pulse delay directly to this sequence's short, tight TE). Channel
    # is unconfirmed for X-Pulse specifically.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rails (see estimate_duty_cycles() docstring: these are
    # conservative, user-adjustable placeholders, not vendor-confirmed specs)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Phases -- simple 2-step CYCLOPS-style cycling (no refocusing pulse to
    # phase-cycle here)
    PH1 = Parameter("PH1", "0,180", ParameterTypes.String, "X Pulse Phase")
    PHRX = Parameter("PHRX", "0,180", ParameterTypes.String, "Acquisition Phase")

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
    Channel2SetFrequency(10,Frequency)
    Channel2RestartSynth(10)

    ReceiverFilter = Filter(P.Filter)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell

    # Total acquisition-window duration -- Receiver2()'s 'duration' argument
    # is the TOTAL window (points*dwell), NOT the per-point dwell.
    ReceiverTime = DW*(points+1)

    times = np.arange(0, points*DW, DW) / 1.0e6

    # ---- Timeline -----------------------------------------------------
    # pulse() elapsed time: Channel2SetBasePhase(1)+Transmit2BlankingOn(1)
    # +TXEnableTime+P90+Transmit2BlankingOff(1) ~= P90+TXEnableTime+3.
    PulseElapsed = P.P90 + P.TXEnableTime + 3
    # Fixed instruction cost of Channel2SetBasePhase(1,0)+Receiver2Phase(1,..)
    # between the end of the pulse and the start of the dead-time waits.
    PostPulseOverhead = 2
    # Gradient plateau (hold) duration: from the end of the ramp-up to the
    # end of acquisition -- i.e. covers the pulse, the fixed dead times, and
    # the acquisition itself, so the gradient never moves during any of it.
    Hold = PulseElapsed + PostPulseOverhead + P.Dead1 + ReceiverFilter.dead_time + ReceiverTime
    # Trailing pad so the RF/acquisition branch's total duration exactly
    # matches the gradient branch's total duration (which has an extra
    # ramp-down + settle tail after Hold) -- required for the two branches
    # inside 'with parallel:' to stay synchronised.
    TrailingPad = P.RampTime + P.GradSettle + 1

    # Achievable TE (pulse-centre to acquisition-start) at the current
    # parameters -- logged so you can see exactly what you're getting.
    MinimumTE = P.P90/2.0 + PostPulseOverhead + P.Dead1 + ReceiverFilter.dead_time
    comms.log("MinimumTE (pulse-centre to acquisition-start) = {0:.2f} us "
              "-- to shorten further, reduce P90, TXEnableTime, and/or "
              "Dead1, or widen Filter (lowers ReceiverFilter.dead_time)."
              .format(MinimumTE))

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

    with sequential:

        # Mains-lock trigger, if enabled -- must come before the first pulse
        # event (manual 3.7.19). OFF by default; see mains_lock_trigger().
        mains_lock_trigger(P, comms)

        Transmit2SelectPort(1,1)
        Transmit2LPEnable(1,0)
        Delay(10000) # This is required for changing the relay state from tune mode
        Transmit2SetScale(5, P.TXAmplitude)

        Receiver2Preamp(128, P.ReceiverAttenuation)
        Receiver2Filter(200, ReceiverFilter)

        # Gradient Setup -- the ONLY place FPX/FPY/FPZ scale the gradients
        GradientMatrix(200,Matrix.T)

        Phases.Reset()

    # ---- Probe/gradient-calibration report (once, outside the scan loop) -
    report_probe_gradient(P, comms)

    # ---- Duty-cycle warnings (computed once, outside the scan loop) ------
    rf_on_time = P.P90 + P.TXEnableTime
    grad_on_time = P.PreGrad + 2*P.RampTime + Hold + P.GradSettle
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver2FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            # Sanity check only (Hold/TrailingPad are always >=0 by
            # construction given non-negative Parameters, but this keeps the
            # same safe-guard convention as the other sequences in this
            # family and will catch a mistyped negative Parameter loudly).
            if Hold < 0 or TrailingPad < 0:
                msg = ("FATAL TIMING ERROR: Hold={0:.2f} us, TrailingPad="
                       "{1:.2f} us -- one of these is negative, check "
                       "P90/TXEnable/Dead1/RampTime/GradSettle.").format(Hold, TrailingPad)
                comms.log(msg)
                raise ValueError(msg)

            with parallel:
                with sequential:
                    # Ramp the read gradient up and let it settle BEFORE the
                    # pulse, hold it flat through pulse+dead-time+acquisition,
                    # then ramp down. This is what makes TE short: there is
                    # no gradient switching event between the pulse and the
                    # start of acquisition.
                    apply_gradient(P.Axis, Hold, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                with sequential:
                    # Wait for the gradient to ramp up and reach its plateau
                    Delay(P.PreGrad + P.RampTime)
                    # Short, broadband excitation pulse (gradient already on
                    # and flat -- this is what makes this a UTE-style
                    # excitation rather than a slice-selective one: it
                    # excites everything within the coil at once, projected
                    # by the already-active read gradient)
                    pulse(P.P90, ph["PH1"], P.TXEnableTime)
                    # Minimal fixed dead time -- this, plus P90/2, IS TE
                    Channel2SetBasePhase(1,0)
                    Receiver2Phase(1, ph["PHRX"])
                    Delay(P.Dead1)
                    Delay(ReceiverFilter.dead_time)
                    Receiver2(ReceiverTime, points)
                    # Pad so this branch's total matches the gradient
                    # branch's total (which has a trailing ramp-down+settle)
                    Delay(TrailingPad)

        Delay(9.0e4)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        TX1.setup_receive(points,32,scan=seqScans,circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX1.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 04/08/26 - Initial version. Built directly on the working
#    gradient/hard-pulse helper functions from LWG_1D-Image-Echo_H.py
#    (apply_gradient, get_gradient_functions, PhasesManager/CallBack1D
#    scaffolding, safe_delay guard, duty-cycle warnings, dead-time
#    handling), but with a completely different pulse sequence: no
#    refocusing pulse, gradient on and settled before the excitation pulse
#    rather than switched during TAU, acquisition starting as soon as fixed
#    hardware dead times allow. Gradient/RF branch timing matched by exact
#    construction (Hold/TrailingPad defined so both branches' totals are
#    equal by arithmetic, not by a separately hand-verified formula on each
#    side -- see comments above Hold/TrailingPad in run()).
# 2. Claude - 04/08/26 - MAINS LOCK: added UseMainsLock/MainsLockChannel
#    parameters and mains_lock_trigger(), called at the very start of the
#    first 'with sequential:' block. Defaults to OFF -- see
#    mains_lock_trigger() docstring; particularly important for THIS
#    sequence since mains triggering would add a variable delay directly
#    onto its short, tight TE. MainsLockChannel is a Parameter, not
#    hard-coded, since the manual doesn't list an X-Pulse row.
# 3. Claude - 06/08/26 - PROBE SELECTOR: added a Probe Parameter ('HFX'
#    default, or 'LOWGAMMA') plus MAXGRAD_TABLE/report_probe_gradient(),
#    since there is no way to detect which probe is fitted from software.
#    report_probe_gradient() logs the selected probe's max gradient
#    strength (G/cm) for the current GradAxis once per run, and -- because
#    Probe is a Parameter -- it is auto-recorded in the resulting JCAMP
#    file's SpinFlow block, so leonmr/xpulse_imaging.py can convert the
#    acquired Hz axis to mm without the operator re-entering which probe
#    was used. MAXGRAD_TABLE currently has real calibration numbers for
#    'HFX' (x=11.879, y=11.978, z=57.915 G/cm, measured/averaged) and
#    placeholder None entries for 'LOWGAMMA' (not yet calibrated) -- keep
#    this table in sync with leonmr/xpulse_imaging.py's own
#    GRADIENT_CALIBRATION dict if either is updated.
# 4. Claude - 10/08/26 - REFINED time_calculation(): mirrors
#    LWG_1D-Image-UTE_H.py's identical fix -- this file already included
#    RD (unlike most other pp files in this repo), but omitted PreGrad/
#    RampTime/GradSettle, a small (~0.1%) systematic underestimate. Now
#    derives the exact gradient-branch total using the same Hold/
#    PulseElapsed logic already defined in run(), + RD +
#    Receiver2FilterFlush. Verified via the mock harness against
#    LWG_1D-Image-UTE_H.py's result at matching default Parameters
#    (2006766.32us, identical).
#
# -----------------------------------------------------------------------------
