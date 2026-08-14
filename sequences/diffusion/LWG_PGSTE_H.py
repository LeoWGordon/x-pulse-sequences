#-------------------------------------------------------------------------------
# Name:        LWG_PGSTE_H.py
# Purpose:     Pulsed-Gradient Stimulated Echo (PGSTE / Tanner stimulated
#              echo) for translational-diffusion measurement: hard-90 - TAU
#              (gradient) - hard-90 - TM - hard-90 - TAU (gradient) - echo.
#              X-Pulse Broadband Benchtop NMR Spectrometer (1H/19F channel).
#              NO IMAGING / NO SLICE SELECTION -- single-shot, non-localised
#              diffusion measurement (spectroscopic use, or as a bench test
#              of a b-value/gradient calibration before building a localised
#              diffusion-imaging sequence on top of it).
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     06/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.4
#
# Design notes:
#  - Structurally a NEAR-VERBATIM port of the vendor's PGSE_H.py (three-pulse
#    variant), with the refocusing 180 pulse replaced by a second hard-90
#    (storage pulse), a mixing delay TM inserted, and a third hard-90
#    (readout/restore pulse) added before the second gradient + acquisition.
#    Reuses PGSE_H.py's GradientMatrix/Gradient3/ramp-and-settle pattern
#    verbatim for both gradient pulses, since that pattern is confirmed
#    working on this instrument (see PGSE_H.py); only the RF/phase/delay
#    skeleton around it changes.
#  - WHY A STIMULATED ECHO (vs. the simpler PGSE two-pulse spin echo): during
#    TM the coherence is stored along +/-Z as longitudinal magnetization, so
#    it decays with T1 rather than T2 during the mixing period. For samples
#    with T2 << T1 (common for restricted/slow diffusion, viscous fluids,
#    or samples with fast transverse relaxation), this lets you probe much
#    longer diffusion times (Delta) than a spin echo can, at the cost of
#    losing half the signal (only the cosine-modulated component survives
#    the second 90, i.e. an intrinsic factor of ~2 versus PGSE, on top of
#    T1 decay during TM instead of no decay).
#  - GRADIENT POLARITY: both diffusion-encoding gradient pulses use the SAME
#    sign and magnitude (G1), matching PGSE_H.py's convention -- but the
#    underlying reason is different. In PGSE, the intervening 180 inverts
#    the accumulated phase, so equal-sign gradients refocus it. Here there
#    is no inversion pulse; instead, the second hard-90 (P2, the "storage"
#    pulse) only stores the COSINE component of the phase accumulated
#    during the first gradient (cos(phi) = cos(-phi), i.e. an even
#    function), so the SIGN of the phase accumulated by the first gradient
#    is effectively discarded at storage. The third hard-90 (P3) restores
#    that cosine-modulated Z-magnetization to the transverse plane, where
#    it then acquires phase from the SECOND gradient of the SAME sign;
#    stationary spins refocus at the same total phase (0) as in PGSE. This
#    is the standard Tanner (1970) / Stejskal-Tanner PGSTE convention --
#    do not flip G1's sign between the two gradient pulses.
#  - TAU is used SYMMETRICALLY: it is the same total delay (encode-gradient
#    period) on BOTH sides of TM -- P1-to-P2 and P3-to-acquisition are both
#    length TAU, with an identical PreGrad/Gradient/GradSettle sub-structure
#    on each side (see D2a/D2b and D3a/D3b below). This is required for the
#    two encode/decode gradient pulses to have identical duration (delta) --
#    an asymmetric implementation would break the basic Stejskal-Tanner
#    attenuation formula, which assumes matched delta on both sides.
#    Nominal diffusion time (gradient-centre to gradient-centre):
#        Delta ~= TM + TAU   (see time_calculation() / run() for the exact
#        sub-delay algebra; this is an approximate description for anyone
#        reasoning about b-values, not the literal code path).
#    Approximate b-value (Stejskal &amp; Tanner 1965): b = (gamma*G1_abs*delta)^2
#    * (Delta - delta/3), where delta = GradientOnTime + 2*RampTime (the
#    effective encoding gradient duration) and G1_abs is your gradient in
#    real units (T/m) -- G1 here is only a DIMENSIONLESS -1.0..1.0 DAC
#    fraction; convert using your own gradient-coil calibration (mT/m per
#    unit G1) before using this formula quantitatively.
#  - NO SPOILER GRADIENT during TM. The Z-storage pathway (what you want)
#    survives TM regardless, but unwanted transverse coherences that leak
#    through P2 (e.g. from imperfect 90 calibration) do NOT get spoiled
#    here -- they are instead suppressed only by the phase cycle below. If
#    you see artefacts that track RFA0 calibration error rather than TM,
#    consider adding a homospoil/spoiler gradient pulse immediately after
#    P2 (before TM) in a future revision.
#  - PHASE CYCLE: taken DIRECTLY from Bruker's vendor-supplied "diffSte"
#    pulse programme (TopSpin 3.5pl7, lists/pp/diffSte -- a standard,
#    field-tested PGSTE phase cycle), converting Bruker's quarter-cycle
#    (0=0,1=90,2=180,3=270) ph1/ph2/ph3/ph31 tables to degrees:
#      Bruker ph1  = 0 0 0 0 2 2 2 2 1 1 1 1 3 3 3 3  (16 steps, P1/excite)
#      Bruker ph2  = 1 3 0 2                           (4 steps, P2/store)
#      Bruker ph3  = 1 3 0 2                           (4 steps, P3/restore)
#      Bruker ph31 = 0 0 2 2 2 2 0 0 3 3 1 1 1 1 3 3   (16 steps, receiver)
#    -> PH1/PHRX below are 16-element degree lists, PH2/PH3 are 4-element
#    degree lists. All four are driven from a SINGLE PhasesManager(P) (per
#    house convention -- see LWG_Selective-Echo_H.py/LWG_1D-Image-Echo-
#    Selective_H.py for the same pattern), rather than one PhaseListContainer
#    per pulse: one Phases.Reset()/Phases.Incd() call per scan returns a
#    dict keyed by Parameter short code (ph["PH1"], ph["PH2"], ph["PH3"],
#    ph["PHRX"]). Each entry is assumed to cycle independently, modulo its
#    OWN length, matching Bruker's own differing-length ph-list semantics
#    (PH2/PH3, length 4, repeat 4x over one full PH1/PHRX, length 16, cycle)
#    -- this assumption is based on this family's other PhasesManager-based
#    pps (e.g. LWG_Selective-Echo_H.py, whose PH1="0,180" and PH2="0" are
#    themselves different lengths and work correctly), not on inspecting
#    PhasesManager's own source, since the framework internals aren't
#    available outside the instrument -- verify against the Pulse Sequence
#    Programming manual if you rely on this for quantitative artefact
#    suppression. NumScans (and DS, for the dummy scans to also complete
#    whole cycles) should be a MULTIPLE OF 16 to complete this cycle; NumScans
#    default below is set to 16 accordingly.
#  - Like PGSE_H.py, this is a HARD-PULSE (non-selective) sequence -- no
#    shaped/selective-pulse machinery is used here. If you need slice- or
#    chemical-shift-selective PGSTE, port the shaped_pulse()/generate_shape()
#    machinery from LWG_Selective-Echo_H.py on top of this file rather than
#    rewriting the gradient/phase-cycle skeleton from scratch.
#  - Mains-lock trigger defaults OFF (UseMainsLock=0), per house convention
#    -- see mains_lock_trigger() below.
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
            self.comms.log("seqTime LWG_PGSTE_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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
               "(negative). Increase TAU (or TM), or decrease P90/PreGrad/"
               "GradientOnTime/RampTime/GradSettle, then retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger (ExternalTrigger[n]()) before the
    first pulse event, per Pulse Sequence Programming User Manual 01-U-049
    section 3.7.19. OFF (UseMainsLock=0) by default -- adds a non-
    deterministic 0-to-one-mains-half-cycle delay before the first event of
    every scan, which is undesirable unless your instrument specifically
    needs mains-synchronous triggering."""
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
    """Log estimated RF and gradient duty cycle for this scan and warn if
    either exceeds a (user-adjustable) threshold Parameter.

    IMPORTANT: MaxRFDuty/MaxGradDuty are conservative, editable placeholders,
    NOT vendor-confirmed ratings -- no public duty-cycle spec for the
    X-Pulse RF/gradient amplifiers has been available; confirm real limits
    with Oxford Instruments before relying on this for unattended runs."""
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR
    grad_duty = grad_on_time / TR
    comms.log("Estimated duty cycle -- RF: {0:.3%} (limit {1:.1%}), Gradient: {2:.3%} (limit {3:.1%})"
              .format(rf_duty, P.MaxRFDuty, grad_duty, P.MaxGradDuty))
    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD or shorter P90, or confirm "
                  "with Oxford Instruments that this is within the "
                  "transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} exceeds "
                  "MaxGradDuty ({1:.2%}). Consider a longer RD or shorter "
                  "GradientOnTime, or confirm with Oxford Instruments that "
                  "this is within the gradient amplifier's rated duty cycle "
                  "before running unattended.".format(grad_duty, P.MaxGradDuty))


# Per-axis max gradient strength (G/cm) for probes that can be fitted to
# this magnet. There is no way to detect which probe is mounted from
# software, so this is a MANUAL selector (Probe Parameter below) -- the
# operator must set it to match what is actually on the magnet. KEEP THIS
# TABLE IN SYNC with leonmr/xpulse_imaging.py's own GRADIENT_CALIBRATION
# dict (duplicated here rather than imported, since pulse programmes can't
# import external Python modules -- same house convention used for
# safe_delay()/mains_lock_trigger() across this family). This sequence's
# diffusion gradient is hard-wired to the z-axis (Matrix = diag(0,0,1) in
# run()), so there is no GradAxis Parameter here -- axis is always 'z'.
MAXGRAD_TABLE = {
    "HFX": {"x": 11.879, "y": 11.978, "z": 57.915},   # G/cm, measured/averaged 06/08/2026
    "LOWGAMMA": {"x": None, "y": None, "z": None},     # NOT YET CALIBRATED
}

def report_probe_gradient(P, comms, axis="z"):
    """Log the max gradient strength (G/cm) for the currently-selected Probe
    on this sequence's (fixed z-axis) diffusion gradient. Since Probe is a
    Parameter, this also gets auto-recorded in the resulting JCAMP file's
    SpinFlow parameter block -- letting leonmr/xpulse_imaging.py convert Hz
    to mm later without the operator having to remember/re-enter which
    probe was fitted for a given experiment. Warns (does not fail) if the
    fitted probe's z-axis has no calibration yet."""
    probe_key = str(P.Probe).upper().replace('/', '').replace('-', '').replace(' ', '').replace('_', '')
    axis_key = str(axis).strip().lower()
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


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    DW = ReceiverFilter.dwell

    # Matches run()'s actual per-scan timing exactly (this file already
    # included RD -- unlike most other pp files in this repo, see
    # changelog for the wider audit -- but the previous D3/PreGrad/
    # GradWidth bookkeeping double-counted the two TAU periods' explicit
    # gradient delays against D3's own "rest of TAU" subtraction, a real
    # ~100-150us/scan error, small next to RD but non-zero). Purely
    # sequential (no 'with parallel:' here), so this is just every
    # Delay()/duration argument in the per-scan block summed exactly, then
    # algebraically simplified: the three P90 pulses and TXEnableTime
    # delays partially cancel against the TM/D3a/D3b subtractions (each of
    # which nets out a P90/2 or TXEnableTime term), leaving TM - P90 +
    # 2*TXEnableTime net; the two TAU periods contribute 2*Tau net (their
    # own PreGrad/GradWidth/GradSettle terms cancel exactly against D3a/
    # D3b's subtraction of the same terms); 208us is the sum of every
    # small fixed-duration instruction (Channel1SetBasePhase/
    # Transmit1Blanking/Gradient3SlewRate/Receiver1FilterFlush/etc.).
    t_scanTime = (P.RecycleDelay + 2*P.Tau + P.TM - P.P90 + 2*P.TXEnableTime
                  + ReceiverFilter.group_delay + P.ReceiverPoints*DW + 208.0)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Pulsed-Gradient Stimulated Echo (PGSTE / Tanner stimulated "
                "echo) on the {H/F} channel: hard-90 - TAU(gradient) - "
                "hard-90 - TM - hard-90 - TAU(gradient) - echo. Used for "
                "measurement of diffusion constants where T2 << T1 makes a "
                "spin-echo PGSE sequence (PGSE_H.py) impractical for the "
                "diffusion time (Delta) you need, at the cost of half the "
                "signal and T1 (not T2) decay during the TM mixing delay. "
                "Non-localised (no imaging/slice-selection).")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,TM,G1"

    return basic

@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_PGSTE_H", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 59.7, ParameterTypes.Double, "H/F Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "H/F Freq Offset [Hz]")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "H/F TX Freq Offset [ppm]")

    ReceiverPoints = Parameter("NP", 32768, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 30, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "5000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans -- NumScans (and DS) should be a MULTIPLE OF 16 to complete the
    # Bruker-derived phase cycle below (PH1/PHRX are 16-step lists; PH2/PH3
    # are 4-step lists that repeat 4x within one PH1/PHRX cycle).
    NumScans = Parameter("NS", 16, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 16, ParameterTypes.Int32,"Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 0.0, ParameterTypes.Double, "H/F Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 2000000, ParameterTypes.Int32, "Relaxation Delay [s]",
                              RD, min=100000, max=2000000000)

    # Hard pulses -- three identical 90s (excitation, storage, restore), no
    # 180 in this sequence.
    P90 = Parameter("P90", 10.0, ParameterTypes.Double, "H/F 90&#176; Pulse Width [&#956;s]")
    TXAmplitude = Parameter("RFA0", 0.4, ParameterTypes.Double,
                            "H/F TX Power [0.0&#8230;1.0]", RFA, min=0.0, max=1.0)

    # Gradients -- same G1 amplitude used for BOTH the encode and decode
    # gradient pulses (see design notes: DO NOT flip sign between them).
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Diffusion Gradient Strength [-1.0&#8230;1.0]")
    GradientOnTime = Parameter("D71", 4000.0, ParameterTypes.Double, "Gradient Duration (&#948;, plateau only) [&#956;s]")
    RampTime = Parameter("D70", 500.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s]")
    PreGrad = Parameter("D75", 1000.0, ParameterTypes.Double, "Pre-Gradient Time (delay from pulse to gradient start) [&#956;s]")
    # Probe fitted to the magnet -- no automatic detection is possible, so
    # this must be set MANUALLY to match what's actually mounted. Used only
    # for logging/downstream Hz-to-mm conversion (see report_probe_gradient()
    # above and leonmr/xpulse_imaging.py) -- has no effect on this pp's own
    # timing/hardware calls.
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet ['HFX'=calibrated, 'LOWGAMMA'=NOT YET CALIBRATED]")

    # Sequence-specific timing
    Tau = Parameter("TAU", 20000.0, ParameterTypes.Double,
                     "TAU [&#956;s], P1-to-P2 and P3-to-acq (symmetric) -- "
                     "must exceed PreGrad+2*RampTime+GradientOnTime+GradSettle+P90")
    TM = Parameter("TM", 50000.0, ParameterTypes.Double,
                    "Mixing Time [&#956;s], P2-to-P3 delay -- Delta &#8776; TM + TAU")

    # Duty-cycle guard rails
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Mains-lock trigger -- OFF by default (see mains_lock_trigger()).
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock Trigger Channel [1-3, unconfirmed for X-Pulse]")

    # Phases -- taken directly from Bruker's diffSte pulse programme
    # (ph1/ph2/ph3/ph31), converted from quarter-cycle to degrees. See
    # design notes for the full derivation. PH1/PHRX: 16 steps. PH2/PH3:
    # 4 steps (repeat 4x per PH1/PHRX cycle). Attribute name = short code
    # = PhasesManager dict key, ALL THREE MUST MATCH EXACTLY (e.g.
    # PH1/PH1/ph["PH1"]) -- see the pulse-programme-parameters skill's
    # phase-cycling rule.
    PH1 = Parameter("PH1","0,0,0,0,180,180,180,180,90,90,90,90,270,270,270,270", ParameterTypes.String,"P1 (excitation) RF-Pulse Phase [Bruker diffSte ph1, 16-step]")
    PH2 = Parameter("PH2","90,270,0,180", ParameterTypes.String,"P2 (storage) RF-Pulse Phase [Bruker diffSte ph2, 4-step]")
    PH3 = Parameter("PH3","90,270,0,180", ParameterTypes.String,"P3 (restore) RF-Pulse Phase [Bruker diffSte ph3, 4-step]")
    PHRX = Parameter("PHRX","0,0,180,180,180,180,0,0,270,270,90,90,90,90,270,270", ParameterTypes.String,"Acquisition Phase [Bruker diffSte ph31, 16-step]")

def run(comms):

    P = Parameters

    Matrix = np.array([
            [0,  0,  0],
            [0  ,0,  0],
            [0  ,0,1.0]
            ])

    Frequency = P.FrequencyBase + (P.FrequencyOffset*1.0e-6) + (P.FrequencyBase*P.TxPPM*1.0e-6)
    Channel1SetFrequency(10,Frequency)
    Channel1RestartSynth(10)

    ReceiverFilter = Filter(P.Filter)
    ReceiverTime = ReceiverFilter.dwell*(P.ReceiverPoints+0.99)

    Phases = PhasesManager(P)
    Phases.Reset()

    points = P.ReceiverPoints
    times = np.arange(0, points*ReceiverFilter.dwell, ReceiverFilter.dwell) / 1.0e6

    GradWidth = 2*P.RampTime + P.GradientOnTime
    GradientSlewRate = abs(P.G1/P.RampTime)

    # "Rest of TAU" after the gradient/settle and pulse tails, on EACH side
    # (identical algebra both sides, by design -- see design notes on why
    # TAU must be symmetric). Same -2/-7 fudge-constant convention as
    # PGSE_H.py's D3 (matches that file's proven-working overhead).
    D3a = P.Tau - P.PreGrad - GradWidth - P.GradSettle - (P.P90/2.) - (P.P90/2.) - 7
    D3b = D3a

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

        Transmit1SelectPort(1,1)
        Transmit1LPEnable(1,0)
        Delay(10000) # This is required for changing the relay state from tune mode
        Transmit1SetScale(5, P.TXAmplitude)

        Receiver1Preamp(128, P.ReceiverAttenuation)
        Receiver1Filter(200, ReceiverFilter)

        # Gradient Setup
        GradientMatrix(200,Matrix.T)

        Phases.Reset()

    # ---- Probe/gradient-calibration report ---------------------------------
    report_probe_gradient(P, comms)

    # ---- Duty-cycle warning ------------------------------------------------
    rf_on_time = 3*P.P90
    grad_on_time = 2*GradWidth
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Channel1SetBasePhase(10, ph["PH1"])
            Receiver1FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            # P1 -- excitation 90
            Transmit1BlankingOn(1)
            Delay(P.TXEnableTime)
            Transmit1(P.P90)
            Transmit1BlankingOff(1)

            # First TAU (encode side): PreGrad, gradient, settle, rest
            safe_delay(P.PreGrad-2-(P.P90/2.), "pre-grad-1 wait", comms)
            Gradient3SlewRate(1.0, GradientSlewRate)
            Gradient3(P.RampTime, P.G1)
            Delay(P.GradientOnTime)
            Gradient3(P.RampTime, 0)
            Gradient3(1,0) # Turn Gradient off explicitly for safety
            Delay(P.GradSettle)
            safe_delay(D3a, "first-TAU rest wait", comms)

            # P2 -- storage 90 (stores cos-modulated phase along Z)
            Channel1SetBasePhase(3, ph["PH2"])
            Transmit1BlankingOn(1)
            Delay(P.TXEnableTime)
            Transmit1(P.P90)
            Transmit1BlankingOff(1)

            # TM -- mixing delay, T1 (not T2) relaxation while stored on Z.
            # No spoiler gradient here -- see design notes.
            safe_delay(P.TM - P.TXEnableTime - (P.P90/2.) - (P.P90/2.), "TM mixing wait", comms)

            # P3 -- restore 90 (returns stored Z-magnetization to transverse)
            Channel1SetBasePhase(3, ph["PH3"])
            Transmit1BlankingOn(1)
            Delay(P.TXEnableTime)
            Transmit1(P.P90)
            Transmit1BlankingOff(1)

            # Second TAU (decode side): PreGrad, gradient (SAME sign as
            # above -- see design notes), settle, rest (minus receiver
            # dead-time/group-delay, paid back explicitly at ACQU).
            safe_delay(P.PreGrad-2-(P.P90/2.), "pre-grad-2 wait", comms)
            Gradient3SlewRate(1.0, GradientSlewRate)
            Gradient3(P.RampTime, P.G1)
            Delay(P.GradientOnTime)
            Gradient3(P.RampTime, 0)
            Gradient3(1,0) # Turn Gradient off explicitly for safety
            Delay(P.GradSettle)
            safe_delay(D3b-ReceiverFilter.dead_time-2+ReceiverFilter.group_delay, "second-TAU rest wait", comms)

            # ACQU
            Channel1SetBasePhase(1,0)
            Receiver1Phase(1, ph["PHRX"])
            Delay(ReceiverFilter.dead_time)
            Receiver1(P.ReceiverPoints*ReceiverFilter.dwell, P.ReceiverPoints)

        Delay(9.0e4)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        TX0.setup_receive(points,32,scan=seqScans,circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX0.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 06/08/26 - Initial version. Built at your request as a test of
#    the repo's auto-save-to-correct-folder workflow: a Pulsed-Gradient
#    Stimulated Echo (PGSTE), i.e. the vendor's PGSE_H.py three-pulse
#    diffusion sequence with the refocusing 180 replaced by a second hard-90
#    (storage) + mixing delay TM + third hard-90 (restore), so that T1 (not
#    T2) governs decay during the diffusion time. Same GradientMatrix/
#    Gradient3 ramp-and-settle pattern as PGSE_H.py (proven working on this
#    instrument), reused verbatim for both the encode and decode gradient
#    pulses -- see design notes at top of file for why both use the SAME
#    gradient sign (unlike a naive "just like PGSE" guess might assume,
#    since there's no 180 here to invert the phase). Added house-convention
#    safe_delay()/mains_lock_trigger()/estimate_duty_cycles() helpers (see
#    LWG_Selective-Echo_H.py for the same pattern) on top of the vendor
#    PGSE_H.py skeleton this was ported from. Phase cycle was an ad hoc
#    4-step EXORCYCLE-style cycle on P1/P2 only (flagged as not literature-
#    confirmed).
# 2. Claude - 06/08/26 - PHASE CYCLE REPLACED with the real cycle from
#    Bruker's own "diffSte" pulse programme (TopSpin 3.5pl7 lists/pp/
#    diffSte, supplied by you), converted from Bruker's quarter-cycle
#    (0/1/2/3 = 0/90/180/270 deg) ph1/ph2/ph3/ph31 tables -- see design
#    notes for the full conversion. PH1/PHRX are now 16-step lists, PH2/PH3
#    4-step lists (each Parameter's phase list cycles independently modulo
#    its own length, exactly as Bruker phase lists of differing length do
#    -- no change needed to the run()/PhaseListContainer mechanism, only
#    the default strings). NumScans/DS defaults changed from 4/2 to 16/16
#    to complete this longer cycle.
# 3. Claude - 06/08/26 - PHASE HANDLING SWITCHED from four independent
#    PhaseListContainer objects (one per pulse, PGSE_H.py's pattern) to a
#    single PhasesManager(P) (this family's house convention, per
#    LWG_Selective-Echo_H.py/LWG_1D-Image-Echo-Selective_H.py -- one
#    Phases.Reset()/Phases.Incd() per scan, phases read back via
#    ph["PH1"]/ph["PH2"]/ph["PH3"]/ph["PHRX"]), for consistency with the
#    rest of the repo. No change to the Bruker-derived phase VALUES from
#    item 2 above, or to the resulting per-scan phase sequence -- only how
#    the running phase state is tracked and incremented.
# 4. Claude - 06/08/26 - PROBE SELECTOR: added a Probe Parameter ('HFX'
#    default, or 'LOWGAMMA') plus MAXGRAD_TABLE/report_probe_gradient(),
#    since there is no way to detect which probe is fitted from software.
#    report_probe_gradient() logs the selected probe's max gradient
#    strength (G/cm) for this sequence's fixed z-axis diffusion gradient
#    once per run, and -- because Probe is a Parameter -- it is
#    auto-recorded in the resulting JCAMP file's SpinFlow block, so
#    leonmr/xpulse_imaging.py can convert the acquired Hz axis to mm
#    without the operator re-entering which probe was used. MAXGRAD_TABLE
#    currently has real calibration numbers for 'HFX' (x=11.879, y=11.978,
#    z=57.915 G/cm, measured/averaged) and placeholder None entries for
#    'LOWGAMMA' (not yet calibrated) -- keep this table in sync with
#    leonmr/xpulse_imaging.py's own GRADIENT_CALIBRATION dict if either is
#    updated.
# 6. Claude - 10/08/26 - FIXED time_calculation(): this file already
#    included RD (unlike most other pp files in this repo -- see the
#    wider audit prompted by you reporting the reported time never
#    reflecting real parameters), but the D3/PreGrad/GradWidth bookkeeping
#    double-counted the two TAU periods' explicit gradient delays against
#    D3's own "rest of TAU" subtraction -- a real ~100-150us/scan error,
#    small next to RD but non-zero. Rewrote as an exact term-for-term sum
#    of every Delay()/duration argument in the (purely sequential, no
#    'with parallel:') per-scan block, algebraically simplified to
#    RD + 2*Tau + TM - P90 + 2*TXEnableTime + group_delay + points*DW +
#    208 (fixed small-instruction overhead). Verified via the mock
#    harness (67936352.0us at default Parameters, matching an independent
#    symbolic derivation exactly).
# 7. Claude - 14/08/26 - FIXED a latent "KeyError: 'PH1'" crash (found via
#    the hardened mock harness, after an identical crash was confirmed on
#    real hardware in the newer LWG_PGSTE-WET_H.py, which was forked from
#    this file and inherited the same bug). Root cause: phase-cycle
#    Parameters were declared as P1Phase = Parameter("PH1", ...) -- Python
#    attribute name "P1Phase" != SpinFlow short code "PH1" != the
#    ph["PH1"] key used at Channel1SetBasePhase(10, ph["PH1"]). Every
#    CONFIRMED-WORKING file in this repo (LWG_1D-Image-Echo_H.py and all
#    the imaging/selective sequences) instead declares these with
#    attribute name == short code == PhasesManager dict key, all
#    identical (PH1 = Parameter("PH1", ...)). FIXED by renaming the
#    attributes themselves from P1Phase/P2Phase/P3Phase/RXPhase to
#    PH1/PH2/PH3/PHRX (matching the short code exactly) -- run()'s
#    PhasesManager(P)/Phases.Incd()/ph["PH1"] logic was already correct
#    and unchanged. Verified via the mock harness (previously failed with
#    this exact KeyError; now completes the full scan loop).
#
# -----------------------------------------------------------------------------
