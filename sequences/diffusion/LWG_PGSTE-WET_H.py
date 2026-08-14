#-------------------------------------------------------------------------------
# Name:        LWG_PGSTE-WET_H.py
# Purpose:     Pulsed-Gradient Stimulated Echo (PGSTE / Tanner stimulated
#              echo) for translational-diffusion measurement, with WET
#              (Water suppression Enhanced through T1 effects) solvent
#              suppression prepended to every scan: hard-90 - TAU(gradient)
#              - hard-90 - TM - hard-90 - TAU(gradient) - echo. X-Pulse
#              Broadband Benchtop NMR Spectrometer (1H/19F channel).
#              NO IMAGING / NO SLICE SELECTION -- single-shot, non-localised
#              diffusion measurement of a solute in a strongly-dominant
#              solvent (e.g. water), where without suppression the solvent
#              signal would swamp the receiver dynamic range / dominate the
#              spectrum.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     10/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.0
# status:      draft
#
# Design notes:
#  - THIS FILE vs LWG_PGSTE_H.py: identical PGSTE diffusion timing/gradient/
#    phase-cycle skeleton (see that file's own design notes for the full
#    Tanner-PGSTE derivation, gradient-polarity reasoning, and b-value
#    formula -- not repeated here). The ONLY addition is the WET module
#    (wet_suppression() below), called once per scan, immediately after the
#    RD wait and before the excitation pulse (P1). Set WetOn=0 to recover
#    LWG_PGSTE_H.py's exact behaviour (zero added RF/gradient/time when
#    off) without needing the separate, non-WET file.
#  - WET MODULE: the standard 4-pulse WET solvent-suppression sequence
#    (Ogg, R. J., Kingsley, P. B. &amp; Taylor, J. S., "WET, a T1- and
#    B1-insensitive water-suppression method for in vivo localized 1H NMR
#    spectroscopy", J. Magn. Reson. B 104, 1-10, 1994). Four selectively-
#    shaped pulses tip the TARGET resonance (WetOffsetPPM, typically water)
#    by successive flip angles 81.4/101.4/69.3/161.0 degrees -- a set
#    numerically optimised to be insensitive to both T1 (across a range of
#    values) and B1 (miscalibration) -- each immediately followed by a
#    spoiler gradient of DECREASING strength (WetSpoil1..4, default
#    0.80/0.40/0.20/0.10, halving each time): a later spoiler only has to
#    dephase whatever residual transverse magnetisation an imperfectly-
#    calibrated preceding pulse left behind, not re-spoil the full
#    magnetisation again. This exact flip-angle set and decreasing-spoiler
#    pattern matches this instrument's own vendor sequence, WET-PRE_H.py
#    (Default pps folder) -- WET-PRE_H.py additionally layers on a DANTE
#    presaturation train ahead of WET and uses the vendor's native
#    Transmit1ShapedPulse(width, phase, shapename) hardware call (referencing
#    a named shape in the instrument's on-board shape library, e.g.
#    "GAUSS201"); this file deliberately uses ONLY plain WET (no DANTE, a
#    separate/optional presaturation technique, not what was asked for) and
#    synthesises the Gaussian shape ON THE FLY (generate_shape()/
#    shaped_pulse(), this repo's house convention throughout -- see
#    xpulse-pulse-programmes skill sec.6) rather than depending on a named
#    shape already loaded onto the instrument, so this file is fully self-
#    contained/auditable like every other pp in this repo.
#  - WET PULSE BANDWIDTH / R_GAUSSIAN_90: a plain truncated Gaussian shape
#    (this repo's generate_shape('GAUSSIAN', ...): tau in [-2.5,2.5],
#    amp=exp(-0.5*tau^2)) was Bloch-simulated (leonmr/bloch_sim.py) at a
#    90-degree rotation, 5000us duration: measured FWHM excitation
#    bandwidth = 440 Hz, giving R_GAUSSIAN_90 = FWHM_Hz * duration_s = 2.2
#    -- a duration-independent property of this shape (see
#    LWG_Slice-Selective-PulseAcquire_H.py's R_EBURP2_90 docstring for why
#    R is duration-independent). WetPulseWidth's default (20000us) targets
#    ~110 Hz FWHM (2.2e6/20000), narrow enough to leave nearby resonances
#    on a 60 MHz instrument's ~500Hz total 1H window largely untouched (see
#    the Bloch-simulated EBURP2 selectivity analysis worked through for the
#    frequency-selective pp family) while still being short enough to keep
#    the WET module's own duty-cycle/duration overhead modest. Widen/narrow
#    via WetPulseWidth if your water peak's own linewidth or its separation
#    from nearby resonances calls for it -- R_GAUSSIAN_90 scales linearly
#    (bandwidth_Hz = R_GAUSSIAN_90*1.0e6/WetPulseWidth), same as every
#    other shape in this repo.
#  - WET FLIP-ANGLE AMPLITUDES (WetAmp1..4): computed via the SAME hard-90-
#    calibration + shape-integration-factor formula used throughout this
#    repo (see LWG_Selective-PulseAcquire_H.py's changelog for the full
#    derivation), scaled by each pulse's own target rotation (81.4/101.4/
#    69.3/161.0, vs. 90 for the hard-pulse reference). Cross-checked
#    against a FULL Bloch simulation (not just the linear/small-tip-angle
#    approximation the formula itself uses): achieved flip angles were
#    within 0.1 degree of target for all four pulses, including the
#    161-degree one -- confirms the linear approximation holds well even
#    at large flip angles for this shape/duration. Still only a STARTING
#    POINT: unlike a normal hard-pulse or single-shape calibration (nutate
#    and look for a signal max/null), WET calibration should be done by
#    directly minimising the RESIDUAL water signal while sweeping each
#    WetAmp in turn (or scaling all four together first, then fine-tuning
#    individually) -- a more direct, model-free calibration than nutation,
#    and the standard way WET is calibrated in practice.
#  - WET TARGET FREQUENCY (WetOffsetPPM): set to your water peak's ppm
#    position relative to SF+O1 -- SEPARATE from TxPPM (the main sequence's
#    excitation/acquisition reference), exactly as this instrument's own
#    WET-PRE_H.py keeps its WET/DANTE target frequency (TxPPMP) separate
#    from the main pulse frequency (TxPPM). shaped_pulse() handles the
#    temporary frequency shift/restore internally (same mechanism as every
#    other selective pp in this repo) -- no extra Channel1SetFrequency
#    bookkeeping needed in run() itself.
#  - LP/HP PORT SWITCHING: WET pulses are long and low-power, so -- per
#    house convention (xpulse-pulse-programmes skill sec.6, matching the
#    vendor's own WET-PRE_H.py) -- they run on the LOW-POWER (LP) TX port,
#    while the main PGSTE hard pulses need the HIGH-POWER (HP) port. This
#    means TWO relay switches per scan (LP-in before WET, HP-in after) that
#    LWG_PGSTE_H.py (permanently on HP) doesn't need -- each carries the
#    usual 10000us relay-settle delay (see shaped_pulse()/LP-port sections
#    elsewhere in this repo for the same value/reasoning), both counted in
#    time_calculation() and estimate_duty_cycles() below. Set WetOn=0 to
#    skip the WET block AND both port switches entirely (wet_suppression()
#    returns immediately, before touching the TX port).
#  - Everything else (Tau/TM symmetry, gradient-polarity convention, Bruker
#    diffSte-derived phase cycle, Probe selector/MAXGRAD_TABLE, RD-correct
#    time_calculation(), mains-lock-off default, no-spoiler-during-TM
#    caveat) is carried over UNCHANGED from LWG_PGSTE_H.py -- see that
#    file's design notes for the full reasoning on each.
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
            self.comms.log("seqTime LWG_PGSTE-WET_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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
                  "({1:.2%}). Consider a longer RD or shorter P90/WetPulseWidth, "
                  "or confirm with Oxford Instruments that this is within the "
                  "transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} exceeds "
                  "MaxGradDuty ({1:.2%}). Consider a longer RD or shorter "
                  "GradientOnTime/WetGradTime, or confirm with Oxford "
                  "Instruments that this is within the gradient amplifier's "
                  "rated duty cycle before running unattended.".format(grad_duty, P.MaxGradDuty))


# Per-axis max gradient strength (G/cm) for probes that can be fitted to
# this magnet. There is no way to detect which probe is mounted from
# software, so this is a MANUAL selector (Probe Parameter below) -- the
# operator must set it to match what is actually on the magnet. KEEP THIS
# TABLE IN SYNC with leonmr/xpulse_imaging.py's own GRADIENT_CALIBRATION
# dict (duplicated here rather than imported, since pulse programmes can't
# import external Python modules -- same house convention used for
# safe_delay()/mains_lock_trigger() across this family). This sequence's
# diffusion gradient (and the WET spoilers) are hard-wired to the z-axis
# (Matrix = diag(0,0,1) in run()), so there is no GradAxis Parameter here --
# axis is always 'z'.
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


# ---- WET on-the-fly shape synthesis ---------------------------------------
# Bandwidth-time product for a plain truncated GAUSSIAN 90-degree pulse (see
# design notes at top of file for the Bloch-simulation derivation:
# leonmr/bloch_sim.py, FWHM=440Hz at 5000us -> R=2.2). Used only to size
# WetPulseWidth's default for a target excitation bandwidth.
R_GAUSSIAN_90 = 2.2

def generate_shape(shape_name, n_steps):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- see LWG_Selective-PulseAcquire_H.py for the full family
    docstring; duplicated here so this file stays self-contained. Only
    GAUSSIAN/SINC are supported -- WET doesn't need the BURP-family
    coefficient machinery used by the frequency-selective pp family."""
    shape_name = shape_name.upper()

    if shape_name == 'GAUSSIAN':
        tau = np.linspace(-2.5, 2.5, n_steps)
        amp = np.exp(-0.5*tau**2)
        phase = np.zeros(n_steps)

    elif shape_name == 'SINC':
        tau = np.linspace(-5.0, 5.0, n_steps)
        with np.errstate(divide='ignore', invalid='ignore'):
            raw = np.where(tau == 0, 1.0, np.sin(np.pi*tau)/(np.pi*tau))
        window = 0.54 + 0.46*np.cos(np.pi*tau/5.0)
        raw = raw*window
        amp = np.abs(raw)
        phase = np.where(raw < 0, 180.0, 0.0)

    else:
        raise ValueError("Unknown WET shape '{0}'. Known: GAUSSIAN, SINC."
                          .format(shape_name))

    peak = np.max(np.abs(amp))
    if peak > 0:
        amp = amp/peak
    return amp, phase


SHAPED_PULSE_FIXED_OVERHEAD = 14
# Exact fixed (non-parallel) instruction cost of shaped_pulse(), in us -- see
# LWG_1D-Image-Echo-Selective_H.py's identical constant for the full
# instruction-by-instruction derivation.

def shaped_pulse(duration, phase, shape, amplitude, txenabletime,
                  base_frequency, pulse_offset):
    """Apply a shaped RF pulse, synthesised on the fly -- identical
    mechanism to LWG_Selective-PulseAcquire_H.py's shaped_pulse() (no
    burp_coeffs_A/B args here, since only GAUSSIAN/SINC are supported)."""
    n_steps = int(round(duration))
    if n_steps < 4:
        raise ValueError("shaped_pulse duration too short for step-wise "
                          "synthesis (need >=4 us so the shape has enough "
                          "points to be meaningful): got {0} us"
                          .format(duration))
    amp_profile, phase_profile = generate_shape(shape, n_steps)

    Channel1SetFrequency(1, base_frequency + pulse_offset*1.0e-6)
    Transmit1SetScale(5, amplitude)
    Channel1SetBasePhase(5, phase)
    Transmit1BlankingOn(1)
    Delay(txenabletime)
    with parallel:
        with sequential:
            Transmit1(float(n_steps))
        with sequential:
            for a in amp_profile:
                Transmit1SetScale(1.0, amplitude*float(a))
        with sequential:
            for p in phase_profile:
                Channel1SetBasePhase(1.0, phase + float(p))
    Transmit1BlankingOff(1)
    Channel1SetFrequency(1, base_frequency)


def wet_suppression(P, comms, Frequency):
    """Emit one WET (Water suppression Enhanced through T1 effects,
    Ogg/Kingsley/Taylor, JMR B 104, 1-10, 1994) solvent-suppression module:
    4 selectively-shaped pulses at the standard 81.4/101.4/69.3/161.0 deg
    flip-angle set, each immediately followed by a spoiler gradient of
    DECREASING strength -- see design notes at top of file. Call ONCE per
    scan, before the main excitation pulse; emits hardware instructions
    directly (no return value). Switches TX to the LOW-POWER port for the
    duration and back to HIGH-POWER afterward (ready for the sequence's own
    hard pulses) -- returns immediately, before touching the TX port at
    all, if WetOn=0."""
    if int(P.WetOn) == 0:
        return

    WetOffsetHz = P.WetOffsetPPM * P.FrequencyBase

    Transmit1SelectPort(1, 0)
    Transmit1LPEnable(1, 1)
    Delay(10000)  # relay settle -- LP-port house convention, see design notes

    for amp, spoil in ((P.WetAmp1, P.WetSpoil1), (P.WetAmp2, P.WetSpoil2),
                        (P.WetAmp3, P.WetSpoil3), (P.WetAmp4, P.WetSpoil4)):
        shaped_pulse(P.WetPulseWidth, P.WetPhase, P.WetShape, amp,
                     P.TXEnableTime, Frequency, WetOffsetHz)
        Gradient3SlewRate(1.0, abs(spoil/P.WetRampTime))
        Gradient3(P.WetRampTime, spoil)
        Delay(P.WetGradTime)
        Gradient3(P.WetRampTime, 0)
        Gradient3(1, 0)  # Turn Gradient off explicitly for safety
        Delay(P.WetGradSettle)
        safe_delay(P.WetInterDelay, "WET inter-pulse wait", comms)

    Transmit1SelectPort(1, 1)
    Transmit1LPEnable(1, 0)
    Delay(10000)  # relay settle back to HP for the hard pulses below


# Precomputed 10/08/2026 for the standard WET flip-angle set (81.4/101.4/
# 69.3/161.0 deg) via the same hard-90-calibration + shape-integration-
# factor formula used throughout this repo (generate_shape('GAUSSIAN',
# 20000) -> integration factor 0.4951; RFAsh_i = HardAmplitude*HardWidth*
# flip_i/(WetPulseWidth*90*IntegFactor)/LPMaxFraction) -- see design notes
# at top of file for the full derivation and Bloch-simulation cross-check.
# Deliberately NOT evaluated at import time -- SpinFlow executes this whole
# file just to LOAD it into the parameter panel, so any on-the-fly shape
# synthesis at module scope runs on EVERY load, not just every scan (see
# the SpinFlow-load-failure fix applied across this repo 10/08/2026).
_DEFAULT_WET_AMP1 = 0.0035    # 81.4 deg
_DEFAULT_WET_AMP2 = 0.00436   # 101.4 deg
_DEFAULT_WET_AMP3 = 0.00298   # 69.3 deg
_DEFAULT_WET_AMP4 = 0.00692   # 161.0 deg


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    DW = ReceiverFilter.dwell

    # Matches run()'s actual per-scan timing exactly -- see LWG_PGSTE_H.py's
    # changelog for the base-formula derivation (RD + 2*Tau + TM - P90 +
    # 2*TXEnableTime + group_delay + points*DW + 208), unchanged here. The
    # WET block adds its own per-scan time on TOP of that (it runs once,
    # after RD, before P1 -- see design notes): 4x(shaped pulse + spoiler
    # gradient + inter-pulse wait), i.e.
    # 4*(WetPulseWidth+TXEnableTime+SHAPED_PULSE_FIXED_OVERHEAD +
    #    2*WetRampTime+WetGradTime+WetGradSettle+WetInterDelay+2), plus the
    # two 10000us LP/HP port-switch relay-settle delays. Zero when
    # WetOn=0, matching wet_suppression()'s early return.
    if int(P.WetOn) != 0:
        wet_time = (4*(P.WetPulseWidth + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD
                        + 2*P.WetRampTime + P.WetGradTime + P.WetGradSettle
                        + P.WetInterDelay + 2.0)
                    + 20000.0)
    else:
        wet_time = 0.0

    t_scanTime = (P.RecycleDelay + 2*P.Tau + P.TM - P.P90 + 2*P.TXEnableTime
                  + ReceiverFilter.group_delay + P.ReceiverPoints*DW + 208.0
                  + wet_time)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Pulsed-Gradient Stimulated Echo (PGSTE / Tanner stimulated "
                "echo) on the {H/F} channel, with a WET water-suppression "
                "module (4 selective pulses + decreasing spoiler gradients, "
                "targeting WetOffsetPPM) prepended to every scan: [WET] - "
                "hard-90 - TAU(gradient) - hard-90 - TM - hard-90 - "
                "TAU(gradient) - echo. Used for measurement of diffusion "
                "constants where T2 &lt;&lt; T1 makes a spin-echo PGSE sequence "
                "impractical for the diffusion time (Delta) you need, at "
                "the cost of half the signal and T1 (not T2) decay during "
                "the TM mixing delay. Non-localised (no imaging/slice-"
                "selection). Set WetOn=0 to disable suppression.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,TM,G1,WetOn,WetOffsetPPM"

    return basic

@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_PGSTE-WET_H", ParameterTypes.String, "Sequence Name")

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
    RecycleDelay = Parameter("RD", 2000000, ParameterTypes.Int32, "Relaxation Delay [s] -- quiet wait BEFORE the WET module starts",
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

    # ---- WET (Water suppression Enhanced through T1 effects) -------------
    # Standard 4-pulse WET module (Ogg, Kingsley &amp; Freeman 1994 -- see
    # design notes), prepended to every scan before P1. Set WetOn=0 to
    # disable entirely (zero added time/RF/gradient when off).
    WetOn = Parameter("WetOn", 1, ParameterTypes.Int32, "Enable WET water suppression [1=On(default), 0=Off]")
    WetOffsetPPM = Parameter("WetOffsetPPM", 4.7, ParameterTypes.Double, "*** WET target [ppm] rel. SF+O1 -- SET to your water peak ***")
    WetShape = Parameter("WetShape", "GAUSSIAN", ParameterTypes.String, "WET Pulse Shape [GAUSSIAN(default)/SINC]")
    WetPhase = Parameter("WetPhase", 0.0, ParameterTypes.Double, "WET Pulse Phase [&#176;] -- same for all 4, spoiled after each so phase doesn't matter")
    WetPulseWidth = Parameter("WetP90sh", 20000.0, ParameterTypes.Double, "WET Shaped Pulse Width [&#956;s] -- sets excitation BW (~110Hz FWHM default, R_GAUSSIAN_90=2.2)")
    WetAmp1 = Parameter("WetAmp1", _DEFAULT_WET_AMP1, ParameterTypes.Double,
                        "WET Pulse 1 TX Power [0&#8230;1 of LP max] -- 81.4&#176; flip, CALIBRATE by nulling water", RFA, min=0.0, max=1.0)
    WetAmp2 = Parameter("WetAmp2", _DEFAULT_WET_AMP2, ParameterTypes.Double,
                        "WET Pulse 2 TX Power [0&#8230;1 of LP max] -- 101.4&#176; flip, CALIBRATE by nulling water", RFA, min=0.0, max=1.0)
    WetAmp3 = Parameter("WetAmp3", _DEFAULT_WET_AMP3, ParameterTypes.Double,
                        "WET Pulse 3 TX Power [0&#8230;1 of LP max] -- 69.3&#176; flip, CALIBRATE by nulling water", RFA, min=0.0, max=1.0)
    WetAmp4 = Parameter("WetAmp4", _DEFAULT_WET_AMP4, ParameterTypes.Double,
                        "WET Pulse 4 TX Power [0&#8230;1 of LP max] -- 161.0&#176; flip, CALIBRATE by nulling water", RFA, min=0.0, max=1.0)
    WetSpoil1 = Parameter("WetSpoil1", 0.80, ParameterTypes.Double, "WET Spoiler 1 Gradient Strength [-1.0&#8230;1.0] -- after pulse 1")
    WetSpoil2 = Parameter("WetSpoil2", 0.40, ParameterTypes.Double, "WET Spoiler 2 Gradient Strength [-1.0&#8230;1.0] -- after pulse 2")
    WetSpoil3 = Parameter("WetSpoil3", 0.20, ParameterTypes.Double, "WET Spoiler 3 Gradient Strength [-1.0&#8230;1.0] -- after pulse 3")
    WetSpoil4 = Parameter("WetSpoil4", 0.10, ParameterTypes.Double, "WET Spoiler 4 Gradient Strength [-1.0&#8230;1.0] -- after pulse 4")
    WetRampTime = Parameter("WetD70", 500.0, ParameterTypes.Double, "WET Spoiler Gradient Ramp Time [&#956;s]")
    WetGradTime = Parameter("WetD71", 2000.0, ParameterTypes.Double, "WET Spoiler Gradient Duration [&#956;s]")
    WetGradSettle = Parameter("WetD73", 200.0, ParameterTypes.Double, "WET Spoiler Gradient Settling Duration [&#956;s]")
    WetInterDelay = Parameter("WetD1", 2000.0, ParameterTypes.Double, "WET Inter-Pulse Delay (after each spoiler settles) [&#956;s]")

    # Duty-cycle guard rails
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Mains-lock trigger -- OFF by default (see mains_lock_trigger()).
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock Trigger Channel [1-3, unconfirmed for X-Pulse]")

    # Phases -- taken directly from Bruker's diffSte pulse programme
    # (ph1/ph2/ph3/ph31), converted from quarter-cycle to degrees. See
    # LWG_PGSTE_H.py's design notes for the full derivation. PH1/PHRX: 16
    # steps. PH2/PH3: 4 steps (repeat 4x per PH1/PHRX cycle). The WET
    # module's own WetPhase (above) is fixed, not cycled -- see its
    # description.
    P1Phase = Parameter("PH1","0,0,0,0,180,180,180,180,90,90,90,90,270,270,270,270", ParameterTypes.String,"P1 (excitation) RF-Pulse Phase [Bruker diffSte ph1, 16-step]")
    P2Phase = Parameter("PH2","90,270,0,180", ParameterTypes.String,"P2 (storage) RF-Pulse Phase [Bruker diffSte ph2, 4-step]")
    P3Phase = Parameter("PH3","90,270,0,180", ParameterTypes.String,"P3 (restore) RF-Pulse Phase [Bruker diffSte ph3, 4-step]")
    RXPhase = Parameter("PHRX","0,0,180,180,180,180,0,0,270,270,90,90,90,90,270,270", ParameterTypes.String,"Acquisition Phase [Bruker diffSte ph31, 16-step]")

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
    # (identical algebra both sides, by design -- see LWG_PGSTE_H.py's
    # design notes on why TAU must be symmetric).
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
                "jc_ns=1",
                "jc_wetOn={0}".format(int(P.WetOn)),
                "jc_wetOffsetPPM={0}".format(P.WetOffsetPPM),
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

    # ---- WET settings report -------------------------------------------
    if int(P.WetOn) != 0:
        comms.log("WET water suppression ON -- targeting WetOffsetPPM={0}ppm "
                  "(={1:.2f} Hz rel. SF), WetPulseWidth={2}us (~{3:.1f}Hz "
                  "FWHM excitation bandwidth), flip angles 81.4/101.4/69.3/"
                  "161.0deg via WetAmp1-4 (CALIBRATE by nulling water)."
                  .format(P.WetOffsetPPM, P.WetOffsetPPM*P.FrequencyBase,
                          P.WetPulseWidth, R_GAUSSIAN_90*1.0e6/P.WetPulseWidth))
    else:
        comms.log("WET water suppression OFF (WetOn=0).")

    # ---- Duty-cycle warning ------------------------------------------------
    if int(P.WetOn) != 0:
        wet_rf_on = 4*(P.WetPulseWidth + P.TXEnableTime)
        wet_grad_on = 4*(2*P.WetRampTime + P.WetGradTime)
    else:
        wet_rf_on = 0.0
        wet_grad_on = 0.0
    rf_on_time = 3*P.P90 + wet_rf_on
    grad_on_time = 2*GradWidth + wet_grad_on
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Channel1SetBasePhase(10, ph["PH1"])
            Receiver1FilterFlush(200, ReceiverFilter)
            # RD -- quiet wait BEFORE the WET module starts (see design notes)
            Delay(P.RecycleDelay-9.0e4)

            # ---- WET water suppression (once per scan, before P1) --------
            wet_suppression(P, comms, Frequency)

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
            # No spoiler gradient here -- see LWG_PGSTE_H.py design notes.
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
# 1. Claude - 10/08/26 - Initial version. Built at your request as a WET-
#    enabled variant of LWG_PGSTE_H.py: identical PGSTE diffusion timing/
#    gradient/phase-cycle skeleton (see LWG_PGSTE_H.py's own changelog for
#    its full history), with a standard 4-pulse WET water-suppression
#    module (Ogg, Kingsley &amp; Taylor 1994 -- matching this instrument's own
#    vendor WET-PRE_H.py flip angles/spoiler pattern, but using on-the-fly
#    GAUSSIAN shape synthesis instead of a named on-board shape, and
#    without WET-PRE_H.py's separate DANTE presaturation train) prepended
#    to every scan, before P1. WetAmp1-4 defaults computed via the
#    established hard-90-calibration + shape-integration-factor formula,
#    cross-validated against a full Bloch simulation (leonmr/bloch_sim.py)
#    to confirm the linear approximation holds within 0.1deg even at the
#    161deg flip. R_GAUSSIAN_90=2.2 similarly Bloch-derived (see design
#    notes). WetOn=0 recovers LWG_PGSTE_H.py's exact behaviour with zero
#    added time/RF/gradient. time_calculation()/estimate_duty_cycles()
#    updated to include the WET block's own contribution (including the
#    two LP/HP port-switch relay delays), verified via the mock harness.
#
# -----------------------------------------------------------------------------
