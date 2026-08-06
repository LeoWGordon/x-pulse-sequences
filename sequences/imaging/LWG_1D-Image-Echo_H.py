#-------------------------------------------------------------------------------
# Name:        LWG_1D-Image-Echo_H.py
# Purpose:     One-dimensional MRI / image-echo (frequency-encoded profile)
#              sequence for the Oxford Instruments X-Pulse Broadband Benchtop
#              NMR Spectrometer (1H/19F channel). Axis-selectable (x/y/z).
#
# Author:      CM / LWG (base sequence); revised by Claude for L. Gordon (DTU)
#
# Created:     05/08/2020
# Revised:     04/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     3.0
#
# Basis:       v3.0 is a NEAR-VERBATIM copy of LWG_diffprof_2_H.py, which is
#              confirmed to compile and run on this X-Pulse (unlike the
#              official vendor 1D-profile-{x,y,z}_H.py sequences, which
#              apparently only work inside Application Developer, not here).
#              A v2.0 rewrite of this file that derived the read-gradient
#              timing from NP*Filter (matching the vendor reference's
#              algebra) did NOT compile -- most likely because it called
#              Receiver1(dwell, points) instead of the correct
#              Receiver1(total_acquisition_time, points), see changelog
#              item (i). Rather than re-guess, this version starts from your
#              own proven-working sequence and only makes small, additive
#              changes on top of it (see changelog).
#
# Changes/Modifications: At end of file -- READ THIS before running on a
#              sample.
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
            self.comms.log("seqTime LWG_1D-Image-Echo_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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

    # Rough estimate only (used purely for the comms.log 'Sequence Time'
    # message) -- kept close to LWG_diffprof_2_H.py's own approximation.
    t_scanTime = (P.P90*5 + P.Tau*2 + points*DW)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("1D MRI / image-echo (frequency-encoded profile) sequence on "
                "the {H/F} channel.\nSpin-echo (90-TAU-[90-270 composite "
                "180]-TAU-echo) with a dephase gradient before the "
                "refocusing pulse and a matched readout gradient spanning "
                "the acquisition window. Gradient axis selectable via "
                "GradAxis (x/y/z).")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,D70,D73,D71,D75,G1,GradAxis"

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
    in run() below, NOT here (see changelog item (a))."""
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
    accounts for that true duration (see changelog item (b) /
    RefocusPulseWidth)."""
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
    immediate error instead."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU and/or PreGrad, or decrease "
               "GradientOnTime/RampTime/GradSettle, then retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger (ExternalTrigger[n]()) before the
    first pulse event, per Pulse Sequence Programming User Manual 01-U-049
    section 3.7.19 ("...mains-lock instructions..."). The manual's trigger-
    channel table only lists MQC+ (mains lock = ExternalTrigger2()) and MQR
    (mains lock = ExternalTrigger3()) -- no X-Pulse row is given, so the
    channel is a Parameter for you to confirm/correct rather than hard-coded.

    OFF (UseMainsLock=0) by default for every sequence in this imaging
    family, per your instruction. Reasoning: a mains-lock trigger waits for
    the next AC zero-crossing before firing, adding a non-deterministic
    delay (0 up to one mains half-cycle -- up to ~10 ms @ 50 Hz / ~8.3 ms @
    60 Hz) before the very first RF/gradient event of EVERY scan. That is
    harmless for ordinary spectroscopy but works against exactly what this
    sequence family is trying to guarantee: tight, reproducible scan-to-scan
    timing for gradient-echo imaging. Set UseMainsLock=1 (and confirm
    MainsLockChannel) if you determine your X-Pulse needs mains-synchronous
    triggering for interference suppression.
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
    Instruments (or your instrument's service documentation) and adjust
    accordingly, especially before running long, low-RD, high-NS
    experiments unattended.
    """
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


@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_1D-Image-Echo_H", ParameterTypes.String, "Sequence Name")
#   Basic = Parameter("Basic", "NS,RD,NP,Filter,D71,D74,G1", ParameterTypes.String, "List of basic parameters")

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

    # Decoupling

    # Gradients -- FPX/FPY/FPZ are the ONLY place per-axis calibration is
    # applied (via GradientMatrix() in run()); G1 is the single logical
    # gradient-strength knob used for both encode and readout lobes (see
    # changelog item (a) -- LWG_diffprof_2_H.py applied FPX/FPY/FPZ a
    # second time via 'G1*gradscaler', which double-scales the current).
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    # Soft/Selective pulses

    # Special acquisition

    # Sequence specific
    GradientOnTime = Parameter("D71", 1000.0, ParameterTypes.Double, "Gradient Duration [&#956;s]")
    RampTime = Parameter("D70", 100.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s]")
    Tau = Parameter("TAU", 10000, ParameterTypes.Int32, "Echo Time [&#956;s]")
    PreGrad = Parameter("D75", 100.0, ParameterTypes.Double, "Pre-Gradient Time [&#956;s]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")

    # Mains-lock trigger -- OFF by default for imaging (see
    # mains_lock_trigger() docstring for why). Channel is unconfirmed for
    # X-Pulse specifically; 2 (MQC+'s mains-lock channel) is used as a
    # placeholder default, only relevant if you turn this on.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse -- see manual 3.7.19]")

    # Duty-cycle guard rails (see estimate_duty_cycles() docstring: these are
    # conservative, user-adjustable placeholders, not vendor-confirmed specs)
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

    # Total acquisition-window duration, i.e. what Receiver1() actually
    # needs as its 'duration' argument (NOT the per-point dwell -- see
    # changelog item (i), this was the compile-breaking bug in v2.0).
    ReceiverTime = DW*(points+1)

    times = np.arange(0, points*DW, DW) / 1.0e6

    # Duration of the refocusing element used below (composite 90-270,
    # via ninety270()): TXEnableTime + 4*P90. Used consistently in every
    # surrounding wait so the echo stays centred on TAU regardless of P90.
    RefocusPulseWidth = P.P90*4 + 4 + P.TXEnableTime

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
    # The readout/decode gradient plateau (2*GradientOnTime, ramped over
    # RampTime on each side) should span at least the digitised acquisition
    # window (ReceiverTime), otherwise part of the FID gets sampled while
    # the gradient is off or still ramping -- a source of truncation/
    # 'leakage' artefacts in the profile. GradientOnTime/RampTime are now
    # free parameters (not derived from NP/Filter), so this is a warning
    # you can act on rather than a hard failure.
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

        # Mains-lock trigger, if enabled -- must come before the first pulse
        # event (manual 3.7.19). OFF by default; see mains_lock_trigger().
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

    # ---- Duty-cycle warnings (computed once, outside the scan loop) ------
    rf_on_time = (P.P90 + P.TXEnableTime) + (RefocusPulseWidth)
    grad_on_time = (P.PreGrad + 2*P.RampTime + P.GradientOnTime + P.GradSettle
                    + P.PreGrad + 2*P.RampTime + (P.GradientOnTime*2+P.RampTime*2) + P.GradSettle)
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
                    # 90 delay
                    Delay((P.P90+P.TXEnableTime+3))
                    # Encode Gradient
                    apply_gradient(P.Axis, P.GradientOnTime, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                    # TAU - gradient time
                    safe_delay(P.Tau - (P.GradientOnTime + 2*P.RampTime + P.GradSettle + 2) - RefocusPulseWidth/2.0,
                               "gradient first-TAU wait", comms) # first tau delay
                    # 180 delay
                    safe_delay(RefocusPulseWidth,
                               "gradient refocus-width wait", comms) # length of the 90-270 pulse and initialisation times
                    # TAU
                    safe_delay(P.Tau-P.PreGrad, "gradient second-TAU wait", comms) # second tau delay
                    # Read Gradient
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
                    # before Receiver1 below) rather than left out entirely.
                    # This branch's TOTAL elapsed time is unchanged, so the
                    # gradient branch does not need to be touched.
                    safe_delay(P.Tau - RefocusPulseWidth/2.0 + 1.*P.PreGrad - ReceiverFilter.dead_time,
                               "RF second-TAU wait", comms)
                    # ACQU
                    Channel1SetBasePhase(P.TXEnableTime,0)
                    Receiver1Phase(P.TXEnableTime, ph["PHRX"])
                    # Dead1: probe ring-down time (hardware/probe-specific,
                    # independent of Filter).
                    Delay(P.Dead1)
                    # ReceiverFilter.dead_time: digital-filter settling time
                    # (varies with the Filter bandwidth you choose). Without
                    # this, Receiver1() can start capturing before the
                    # receiver chain has actually settled, so the first part
                    # of the acquired 'signal' is really filter transient /
                    # ring-down, not the NMR signal -- exactly the kind of
                    # thing that shows up as spurious 'leakage' in the FT.
                    Delay(ReceiverFilter.dead_time)
                    Receiver1(ReceiverTime, points)

        Delay(9.0e4)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        TX0.setup_receive(points,32,scan=seqScans,circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX0.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))
        comms.log("Encode Time %s" % (P.GradientOnTime))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. CM  - 05/08/20 - Initial Version for X-Pulse
# 2. AS  - 25/05/21 - "Preview" instruction added
# 3. RJB - 01/10/21 - Updated List of Basic Parameters
# 4. AS  - 22/12/21 - Updated to run RD>=0.1s
# 5. RJB - 27/02/23 - addition of 'def sequence_description():' and 'def sequence_basic():'
# 6. LWG - 10/03/26 - Updated with consolidated functions for gradients, pulses, and phases.
# 7. LWG - 17/03/26 - Updated to a 1D profile experiment using a 90-tau-90-270-acq sequence.
# 8. LWG - ??/??/26 - "This version I will change the timing according to my older sequence"
#                      (-> LWG_diffprof_2_H.py, confirmed to compile on this X-Pulse).
# 9. Claude - 04/08/26 (v2.0, SUPERSEDED) - Rebuilt around timing derived from
#      NP*Filter (matching the official 1D-profile-{x,y,z}_H.py algebra).
#      Reported as not compiling -- most likely cause identified below (i).
#      Kept for reference in git history/chat, not carried forward here.
# 10. Claude - 04/08/26 (v3.0, THIS VERSION) - Rebuilt starting from
#      LWG_diffprof_2_H.py (your confirmed-working sequence) instead, making
#      only small, additive changes on top of it:
#      a) GRADIENT DOUBLE-SCALING: LWG_diffprof_2_H.py computes
#         'gradscaler' (= FPX/FPY/FPZ depending on GradAxis) and passes
#         'P.G1*gradscaler' as the gradient amplitude into apply_gradient(),
#         while GradientMatrix() (built from the same FPX/FPY/FPZ) ALSO
#         scales that amplitude again in hardware. With the shipped
#         defaults (FPX=FPY=FPZ=1.0) this is invisible, but once you
#         calibrate FPZ to a real value (e.g. ~0.17-0.34, see the
#         'optimised parameters' table in your imaging docs) the true
#         current becomes G1*FPZ^2 instead of G1*FPZ. Removed the
#         '*gradscaler' multiply; G1 is now the sole logical gradient
#         knob and FPX/FPY/FPZ calibration is applied exactly once, by
#         GradientMatrix(). (This is the same class of bug patched
#         separately in the vendor's 1D-profile-z_H.py, where the
#         SlewRate calculation was hard-coded to FPX instead of FPZ.)
#      b) READ-GRADIENT COVERAGE CHECK: GradientOnTime (D71) and RampTime
#         (D70) are free parameters, independent of NP/Filter, so nothing
#         stops you from setting an acquisition window (NP*dwell) longer
#         than the readout gradient's flat plateau -- part of the FID would
#         then be sampled with the gradient off/ramping (truncation/
#         'leakage'). Added a comms.log WARNING (not a hard stop, since
#         these are legitimately independent, user-tunable parameters) that
#         compares the two and tells you which parameter to change.
#      c) GUARD RAILS: wrapped every TAU/PreGrad/RampTime-derived Delay()
#         in safe_delay(), which raises a clear, specific error (naming the
#         exact wait) instead of a silent/confusing failure if TAU is too
#         short for your PreGrad/GradientOnTime/RampTime combination.
#      d) DUTY-CYCLE WARNINGS: added estimate_duty_cycles(), logging
#         estimated RF and per-axis gradient duty cycle and warning if they
#         exceed MaxRFDuty/MaxGradDuty (new parameters, conservative
#         user-adjustable placeholders -- no public duty-cycle spec for the
#         X-Pulse 60 MHz RF/gradient amplifiers was available when this was
#         written; confirm real limits with Oxford Instruments).
#      e) Everything else (Parameters defaults, apply_gradient()'s pre_grad
#         argument, pulse()'s Channel1SetBasePhase(1,...), the composite
#         90-270 refocusing pulse and its RefocusPulseWidth-based timing,
#         Receiver1(ReceiverTime, points), PhasesManager spelling) is
#         unchanged from LWG_diffprof_2_H.py.
#      f) SUSPECTED (BUT NOT THE ACTUAL) ROOT CAUSE OF THE v2.0 COMPILE
#         FAILURE: v2.0 called Receiver1(DW, points) -- DW being the
#         per-point dwell, not the total acquisition duration -- which is
#         still fixed here (see Receiver1(ReceiverTime, points) below), but
#         turned out not to be why v2.0 failed to compile; see (g).
#      g) ACTUAL ROOT CAUSE OF THE COMPILE FAILURE: the file had plain
#         Unix LF line endings. The X-Pulse pulse-sequence compiler only
#         accepts Windows CRLF line endings (confirmed: every vendor/
#         working .py in 'Default pps' and your own LWG_diffprof_2_H.py use
#         CRLF; the broken v1.0/v2.0 files used LF). This file is saved
#         with CRLF endings -- if you ever re-save it from an editor that
#         defaults to LF, convert back to CRLF before loading it.
#      h) DEAD TIME: previously only Dead1 (probe ring-down) was waited out
#         before Receiver1(). Added an explicit Delay(ReceiverFilter.dead_time)
#         (the digital filter's own settling time, which depends on the
#         Filter/bandwidth parameter) immediately before Receiver1(), so
#         Receiver1() only starts once BOTH the probe ring-down AND the
#         filter's settling transient are over -- otherwise the first part
#         of the 'signal' acquired is filter transient, not the NMR signal.
#         This time is reserved out of (not added on top of) the preceding
#         TAU wait, so the RF branch's total duration -- and therefore its
#         sync with the gradient branch -- is unchanged.
# 11. Claude - 04/08/26 - MAINS LOCK: added UseMainsLock/MainsLockChannel
#      parameters and mains_lock_trigger(), called at the very start of the
#      first 'with sequential:' block. Defaults to OFF (UseMainsLock=0) for
#      this and every other sequence in the imaging family, per your
#      instruction -- see mains_lock_trigger() docstring for the reasoning
#      (mains triggering adds a non-deterministic pre-scan delay that works
#      against reproducible gradient-echo timing). The manual doesn't list
#      an X-Pulse row for the mains-lock trigger channel (only MQC+/MQR), so
#      MainsLockChannel is left as a Parameter, not hard-coded -- confirm
#      the right channel before ever setting UseMainsLock=1.
#
# -----------------------------------------------------------------------------
