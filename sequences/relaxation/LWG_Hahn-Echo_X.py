#-------------------------------------------------------------------------------
# Name:        LWG_Hahn-Echo_X.py
# Purpose:     Basic Hahn (spin) echo for T2 relaxometry: hard-90 - TAU -
#              hard-180 - TAU - echo, single acquisition per scan, NO
#              gradients/imaging/slice-selection. X-Pulse Broadband Benchtop
#              NMR Spectrometer, on the X (2nd RF) channel -- any nucleus
#              your X probe/coil supports (e.g. 19F, 23Na, 31P...).
#
#              INTENDED USE: run this sequence as a quasi-2D / VC-list sweep
#              of TAU (e.g. 200, 400, 800, 1600, ... us) across SEPARATE
#              scans/experiments -- exactly the "diffprof.lwg-style single-
#              echo TAU array" pattern already used elsewhere in this repo
#              (see LWG_1D-Image-Echo-GradAfterRefocus_H.py's design notes).
#              Echo intensity vs. TE = 2*TAU (approximately -- see
#              time_calculation()/comms.log for the exact predicted value)
#              gives a T2 decay curve via a mono-exponential fit. Compare
#              against LWG_CPMG_X.py (same channel, same phase-cycle/pulse
#              helpers) -- see "HAHN ECHO vs. CPMG" below for why the two
#              give genuinely different measurements, not just two ways to
#              get the same number.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     17/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.0 (draft) -- NOT YET RUN ON HARDWARE. Per this repo's
#              README: no AI-assisted sequence is used for real data
#              collection until a named person has reviewed and validated it
#              on the instrument.
# status: draft
#
# X-CHANNEL CALIBRATION -- IMPORTANT: this is a mechanical H->X port (all
# Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1, per
# the vendor's own H/X sequence pairs, e.g. PGSE_H.py vs PGSE_X.py). SF now
# defaults to 15.01 MHz as a placeholder -- P90/P180/TXAmplitude(RFA0)/
# ReceiverAttenuation(RA) below are still the 1H-calibrated values from
# LWG_Hahn-Echo_H.py and are almost certainly WRONG for whatever nucleus you
# actually put on the X channel. Recalibrate SF, P90, P180, RFA0, and RA for
# your target nucleus (e.g. via a simple pulse-acquire nutation curve) before
# running this on a sample.
#
# Design notes:
#  - This is the SIMPLEST possible echo: one hard 90, one hard 180, one
#    acquisition. No composite refocusing pulse, no gradients, no imaging.
#    Deliberately built as a close sibling of LWG_CPMG_X.py -- both files
#    share the same pulse()/safe_delay()/mains_lock_trigger() helpers, the
#    same _first_gap()/_subsequent_gap() timing formulas, and the same
#    8-step Meiboom-Gill-style phase cycle (PH1/PH2/PHRX), all copied
#    verbatim from the vendor's own CPMG_H.py (Created 10/09/2013, Author
#    AS, supplied by L. Gordon 11/08/2026 -- see
#    sequences/imaging/LWG_CPMG-Image-Echo_H.py's design notes for the full
#    provenance, and LWG_Hahn-Echo_H.py for the original 1H-channel port).
#    This sequence is in fact EXACTLY CPMG_H.py's/LWG_CPMG_X.py's own NECH=1
#    case: one excitation pulse, one gap (_first_gap), one refocusing pulse,
#    one post-refocus wait, one acquisition -- verified algebraically below
#    (time_calculation() reduces to the identical formula LWG_CPMG_X.py's
#    own time_calculation() produces at NECH=1).
#  - POST-REFOCUS WAIT: uses the vendor's own plain-spectroscopy formula,
#    (Tau - P180/2 - dead_time + group_delay - 3), NOT the gradient-echo-
#    aligned formula used by the imaging family (LWG_1D-Image-Echo-
#    GradAfterRefocus_H.py etc.) -- there is no gradient here, so there is
#    no gradient echo to centre; this formula instead compensates for the
#    receiver filter's own group delay so the acquired FID/echo starts
#    cleanly. See LWG_CPMG-Image-Echo_H.py's "WHY THE FINAL ECHO DOESN'T
#    REUSE CPMG_H.py's OWN PRE-ACQUISITION FORMULA" for the full reasoning
#    (written from the imaging sequence's point of view, i.e. why it does
#    NOT reuse this formula -- this file is the case where the vendor's
#    own formula DOES apply, unmodified).
#  - HAHN ECHO vs. CPMG -- why build both: a single Hahn echo with TAU
#    itself swept directly TE-couples molecular diffusion into the decay
#    (bigger TAU -> longer TE -> more diffusion attenuation on top of true
#    T2 decay, via the same b~TAU^3-ish mechanism as a PGSE experiment,
#    even with G1=0 -- background/susceptibility gradients alone are enough
#    at long TAU). CPMG instead keeps TAU SHORT and FIXED and increases TE
#    by adding more refocusing pulses (NECH), which refocuses diffusion
#    attenuation far more effectively (echo spacing, not total TE, sets the
#    diffusion sensitivity). Comparing a Hahn-echo TAU-array T2 against a
#    CPMG NECH-array T2 on the SAME sample is a direct, practical way to
#    see whether your measured "T2" is diffusion-contaminated -- if the two
#    disagree (Hahn echo apparently shorter/faster-decaying), diffusion is
#    the likely explanation, not a genuinely shorter T2.
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
            self.comms.log("seqTime LWG_Hahn-Echo_X Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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


def safe_delay(value, label, comms):
    """Delay() wrapper -- raises a clear, specific error instead of a
    silent/confusing failure if a computed delay would be negative."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU, or decrease P90/P180/TXEnable, "
               "then retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger (ExternalTrigger[n]()) before the
    first pulse event, per Pulse Sequence Programming User Manual 01-U-049
    section 3.7.19. OFF (UseMainsLock=0) by default."""
    if int(P.UseMainsLock) == 0:
        return
    triggers = {1: ExternalTrigger1, 2: ExternalTrigger2, 3: ExternalTrigger3}
    ch = int(P.MainsLockChannel)
    if ch not in triggers:
        comms.log("WARNING: MainsLockChannel={0} is not 1, 2, or 3 -- mains "
                  "lock trigger skipped.".format(ch))
        return
    triggers[ch]()


def estimate_rf_duty_cycle(P, rf_on_time, comms):
    """Log estimated RF duty cycle and warn if it exceeds MaxRFDuty
    (conservative, editable placeholder -- confirm with Oxford Instruments).
    No gradients in this sequence, so RF is the only duty-cycle concern."""
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR
    comms.log("Estimated RF duty cycle: {0:.3%} (limit {1:.1%})".format(rf_duty, P.MaxRFDuty))
    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD or shorter P90/P180, or "
                  "confirm with Oxford Instruments that this is within the "
                  "transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, P.MaxRFDuty))


def pulse(length, phase, txenabletime):
    """Apply a hard pulse of given length and phase. Verbatim from the
    vendor's CPMG_H.py -- see LWG_CPMG-Image-Echo_H.py's design notes."""
    Channel2SetBasePhase(1, phase)
    Transmit2BlankingOn(1)
    Delay(txenabletime)
    Transmit2(length)
    Transmit2BlankingOff(1)


def _first_gap(P):
    """Delay between P1 (excitation) and the refocusing pulse. Verbatim
    from the vendor's CPMG_H.py."""
    return P.Tau - P.TXEnableTime - ((P.P90+P.P180)/2.0) - 12.0


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)

    # Exact term-by-term derivation of run()'s per-scan body. RD is split
    # across Delay(RD-9e4) + a trailing Delay(9e4), summing to exactly RD;
    # Receiver2FilterFlush(200,...) contributes a fixed 200us. Each pulse()
    # call costs TXEnableTime+pulse_width+3 (three 1us duration-arg
    # instructions: Channel2SetBasePhase(1,.)+Transmit2BlankingOn(1)+
    # Transmit2BlankingOff(1)). The post-P1 gap (_first_gap) and the
    # post-refocus wait (Tau-P180/2-dead_time+group_delay-3, verbatim from
    # the vendor's CPMG_H.py) algebraically cancel TXEnableTime and P180
    # entirely, leaving a clean 2*Tau term (TE=2*Tau, by construction) plus
    # a residual P90/2 offset from the excitation pulse -- worked by hand;
    # cross-checked numerically against LWG_CPMG_X.py's own
    # time_calculation() at NECH=1, which must (and does) give the
    # identical result, since this sequence IS that file's NECH=1 case.
    t_pulse_event = (2.0*P.Tau + P.TXEnableTime + P.P90/2.0
                      + ReceiverFilter.group_delay + P.Dead1 + ReceiverTime - 7.0)
    t_scanTime = P.RecycleDelay + 200.0 + t_pulse_event
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Basic Hahn (spin) echo for T2 relaxometry on the X "
                "channel: hard-90 - TAU - hard-180 - TAU - echo, single "
                "acquisition, NO gradients/imaging. ARRAY TAU across "
                "separate scans to build a T2 decay curve (echo intensity "
                "vs. TE=2*TAU). Compare against LWG_CPMG_X.py to check for "
                "diffusion contamination at long TAU. X-CHANNEL: recalibrate "
                "SF/P90/P180/RFA0/RA for your target nucleus first.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,P90,P180"

    return basic

@ParameterBlock
class Parameters:

    # Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_Hahn-Echo_X", ParameterTypes.String, "Sequence Name")

    # General acquisition -- SF defaults to a PLACEHOLDER (15.01 MHz); see
    # X-CHANNEL CALIBRATION note at the top of this file.
    FrequencyBase = Parameter("SF", 15.01, ParameterTypes.Double, "X Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "X Freq Offset [Hz]")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "X TX Freq Offset [ppm]")

    ReceiverPoints = Parameter("NP", 1024, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 34, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "100000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans -- NumScans (and DS) should be a MULTIPLE OF 8 to complete the
    # Meiboom-Gill phase cycle below (PH1/PH2/PHRX are 8-step lists).
    NumScans = Parameter("NS", 8, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 8, ParameterTypes.Int32,"Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 100.0, ParameterTypes.Double, "X Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 2000000, ParameterTypes.Int32, "Relaxation Delay [s] -- SET >= 5x your sample's T1",
                              RD, min=100000, max=2000000000)

    # Hard pulses -- plain 90/180, HP port, no shaping. STILL 1H-CALIBRATED
    # DEFAULTS -- see X-CHANNEL CALIBRATION note at the top of this file.
    P90 = Parameter("P90", 9.58, ParameterTypes.Double, "X 90&#176; Pulse Width [&#956;s]")
    P180 = Parameter("P180", 19.16, ParameterTypes.Double, "X 180&#176; Pulse Width [&#956;s]")
    TXAmplitude = Parameter("RFA0", 0.4, ParameterTypes.Double,
                            "X TX Power [0.0&#8230;1.0], SAME for P90+P180", RFA, min=0.0, max=1.0)

    # *** ARRAY THIS PARAMETER (across separate scans/experiments) to
    # measure T2 -- echo intensity vs. TE=2*TAU gives a mono-exponential
    # decay curve. Keep TAU >= (TXEnable+(P90+P180)/2+12)/1 so _first_gap
    # stays positive (safe_delay() will raise a clear error otherwise). ***
    Tau = Parameter("TAU", 2000.0, ParameterTypes.Double,
                     "&#964; Delay [&#956;s], P1-to-P2 and P2-to-acq (TE=2*TAU) -- *** ARRAY THIS for T2 ***")

    # Duty-cycle guard rail (RF only -- no gradients in this sequence)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Mains-lock trigger -- OFF by default (see mains_lock_trigger()).
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock Trigger Channel [1-3, unconfirmed for X-Pulse]")

    # Phases -- 8-step Meiboom-Gill-style cycle, copied verbatim from the
    # vendor's CPMG_H.py (same cycle LWG_CPMG_X.py uses, for direct
    # comparability). Attribute name = short code = PhasesManager dict key,
    # ALL THREE MUST MATCH EXACTLY -- see the pulse-programme-parameters
    # skill's phase-cycling rule.
    PH1 = Parameter("PH1", "0,0,180,180,90,90,270,270", ParameterTypes.String, "X 90&#176; Pulse Phase")
    PH2 = Parameter("PH2", "90,270,90,270,0,180,0,180", ParameterTypes.String, "X 180&#176; Pulse Phase")
    PHRX = Parameter("PHRX", "0,0,180,180,90,90,270,270", ParameterTypes.String, "Acquisition Phase")

def run(comms):

    P = Parameters

    Phases = PhasesManager(P)
    Phases.Reset()

    Frequency = P.FrequencyBase + (P.FrequencyOffset*1.0e-6) + (P.FrequencyBase*P.TxPPM*1.0e-6)
    Channel2SetFrequency(10,Frequency)
    Channel2RestartSynth(10)

    ReceiverFilter = Filter(P.Filter)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)

    times = np.arange(0, points*DW, DW) / 1.0e6

    TE = 2.0*P.Tau

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
                "jc_te={0}".format(TE),
                    ]
        return jc_data

    comms.log("Sequence Time: {0}".format(time_calculation(Parameters)))
    comms.log("Hahn echo: TAU={0}us, approx TE=2*TAU={1}us ({2:.3f} ms) -- "
              "real hardware overhead adds a small amount on top.".format(P.Tau, TE, TE/1000.0))
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

        Phases.Reset()

    # ---- Duty-cycle warning (RF only, no gradients here) ------------------
    rf_on_time = P.P90 + P.P180
    estimate_rf_duty_cycle(P, rf_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver2FilterFlush(200, ReceiverFilter)
            # RD -- quiet wait before the first pulse.
            Delay(P.RecycleDelay-9.0e4)

            # P1 -- excitation 90
            pulse(P.P90, ph["PH1"], P.TXEnableTime)

            # Gap to the refocusing pulse.
            safe_delay(_first_gap(P), "excitation-to-refocus gap", comms)

            # P2 -- refocusing 180
            pulse(P.P180, ph["PH2"], P.TXEnableTime)

            # Post-refocus wait -- verbatim from the vendor's CPMG_H.py,
            # tuned for a plain spectroscopy acquisition (compensates the
            # receiver filter's own group delay so the echo/FID starts
            # cleanly; NOT the gradient-echo-aligned formula the imaging
            # family uses -- there is no gradient here).
            safe_delay(P.Tau - P.P180/2.0 - ReceiverFilter.dead_time + ReceiverFilter.group_delay - 3,
                       "post-refocus wait", comms)

            # ACQU
            Channel2SetBasePhase(1,0)
            Receiver2Phase(1, ph["PHRX"])
            Delay(P.Dead1)
            Delay(ReceiverFilter.dead_time)
            Receiver2(ReceiverTime, points)

        Delay(9.0e4)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        TX1.setup_receive(points,32,scan=seqScans,circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX1.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 17/08/26 - Initial version. Mechanical X-channel port of
#    LWG_Hahn-Echo_H.py (Channel1/Transmit1/Receiver1/TX0 ->
#    Channel2/Transmit2/Receiver2/TX1, per the vendor's own H/X sequence
#    pairs) -- see that file's changelog for the full derivation/
#    verification of the shared timing formulas. SF defaults to a
#    placeholder 15.01 MHz; P90/P180/RFA0/RA remain the 1H-calibrated
#    values and MUST be recalibrated for your target nucleus (see
#    X-CHANNEL CALIBRATION note at the top of this file). Verified via the
#    mock harness. DRAFT, not yet run on hardware.
#
# -----------------------------------------------------------------------------
