#-------------------------------------------------------------------------------
# Name:        LWG_PGSE-WET_H.py
# Purpose:     Pulsed-Gradient Spin Echo (PGSE / Stejskal-Tanner) for
#              translational-diffusion measurement, with WET (Water
#              suppression Enhanced through T1 effects) solvent suppression
#              prepended to every scan: hard-90 - TAU(gradient) - hard-180
#              - TAU(gradient) - echo. X-Pulse Broadband Benchtop NMR
#              Spectrometer (1H/19F channel). NO IMAGING / NO SLICE
#              SELECTION -- single-shot, non-localised diffusion
#              measurement of a solute in a strongly-dominant solvent (e.g.
#              water), where without suppression the solvent signal would
#              swamp the receiver dynamic range / dominate the spectrum.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     10/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.2
# status:      draft
#
# Design notes:
#  - THIS FILE vs the vendor's PGSE_H.py (Default pps folder) / this repo's
#    own LWG_PGSTE-WET_H.py: a near-verbatim port of PGSE_H.py's two-pulse
#    (90-180) spin-echo diffusion skeleton -- same GradientMatrix/
#    Gradient3/ramp-and-settle pattern, same "SAME gradient sign on both
#    sides" convention (the intervening 180 inverts the phase accumulated
#    by the first gradient, so equal-sign gradients refocus it -- standard
#    Stejskal-Tanner PGSE) -- with house-convention helpers (safe_delay(),
#    mains_lock_trigger(), estimate_duty_cycles(), Probe selector/
#    MAXGRAD_TABLE, PhasesManager, RD-correct term-by-term time_calculation
#    derivation) added on top, matching LWG_PGSTE-WET_H.py, PLUS the SAME
#    WET water-suppression module as that file (see its own design notes
#    for the full WET derivation/literature reference/Bloch-validation --
#    not repeated here; wet_suppression()/generate_shape()/shaped_pulse()
#    below are an EXACT duplicate).
#  - PGSE vs PGSTE (this repo's LWG_PGSTE-WET_H.py): PGSE stores the
#    dephased coherence in the TRANSVERSE plane throughout (via the 180
#    refocusing pulse), so it decays with T2 during the whole diffusion
#    time Delta -- simpler (one fewer pulse, no mixing-time T1 leakage
#    concern, twice the signal of PGSTE for the same conditions) but
#    unusable once Delta needs to exceed T2 (e.g. slow/restricted
#    diffusion, viscous samples, or samples with fast transverse
#    relaxation) -- use LWG_PGSTE-WET_H.py in that case instead.
#  - TAU is used SYMMETRICALLY, exactly as in LWG_PGSTE-WET_H.py: the same
#    total delay (P1-centre-to-P2-centre, and P2-centre-to-acquisition) on
#    BOTH sides of the refocusing pulse, with matched PreGrad/Gradient/
#    GradSettle sub-structure each side (D3a/D3b below) -- required for the
#    two gradient pulses to have identical duration (delta), which the
#    basic Stejskal-Tanner attenuation formula assumes. Echo time TE =
#    2*TAU; diffusion time (gradient-centre to gradient-centre) Delta ~=
#    TAU (see time_calculation()/run() for the exact sub-delay algebra --
#    an approximate description for anyone reasoning about b-values, not
#    the literal code path, exactly as noted in LWG_PGSTE-WET_H.py).
#    Approximate b-value (Stejskal &amp; Tanner 1965): b =
#    (gamma*G1_abs*delta)^2 * (Delta - delta/3), where delta =
#    GradientOnTime + 2*RampTime and G1_abs is your gradient in real units
#    (T/m) -- G1 here is only a DIMENSIONLESS -1.0..1.0 DAC fraction;
#    convert using your own gradient-coil calibration before using this
#    formula quantitatively.
#  - WET MODULE, PORT SWITCHING, WET PULSE BANDWIDTH/AMPLITUDE DERIVATION:
#    identical to LWG_PGSTE-WET_H.py -- see that file's design notes for
#    the full literature reference (Ogg, Kingsley &amp; Taylor 1994),
#    R_GAUSSIAN_90 Bloch-simulation derivation, WetAmp1-4 calibration
#    formula + Bloch cross-check, and LP/HP port-switching reasoning. Set
#    WetOn=0 to disable entirely (zero added time/RF/gradient when off).
#  - Like PGSE_H.py, this is a HARD-PULSE (non-selective) sequence for the
#    diffusion-encoding P1/P2 -- only the WET module uses shaped/selective
#    pulses. If you need slice- or chemical-shift-selective diffusion
#    encoding, port shaped_pulse()/generate_shape() onto P1/P2 as well
#    (see LWG_Selective-Echo_H.py for the pattern) rather than rewriting
#    the gradient/timing skeleton from scratch.
#  - Mains-lock trigger defaults OFF (UseMainsLock=0), per house convention
#    -- see mains_lock_trigger() below.
#  - PHASE CYCLE: a straightforward 2-step EXORCYCLE-style cycle on the
#    refocusing pulse (PH2=90,270) with the receiver following (PHRX=0,180)
#    to select the correctly-refocused echo and reject artefacts from an
#    imperfectly-calibrated 180 -- NOT taken from a specific vendor
#    reference pulse programme (unlike LWG_PGSTE-WET_H.py's Bruker-diffSte-
#    derived cycle) since PGSE_H.py itself only cycles PH2 the same way
#    (P.P2Phase="90,270") with a FIXED PH1/PHRX; extended here to also
#    cycle PHRX so imperfect-180 artefacts are actively rejected rather
#    than just averaged. NumScans should be a multiple of 2.
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
            self.comms.log("seqTime LWG_PGSE-WET_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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
               "(negative). Increase TAU, or decrease P90/P180/PreGrad/"
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
                  "({1:.2%}). Consider a longer RD or shorter P90/P180/"
                  "WetPulseWidth, or confirm with Oxford Instruments that "
                  "this is within the transmitter's rated duty cycle before "
                  "running unattended.".format(rf_duty, P.MaxRFDuty))
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


# ---- WET on-the-fly shape synthesis (EXACT duplicate of
# LWG_PGSTE-WET_H.py -- see that file's design notes for the full
# literature reference and derivation) ---------------------------------
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
    Ogg/Kingsley/Taylor, JMR B 104, 1-10, 1994) solvent-suppression module
    -- EXACT duplicate of LWG_PGSTE-WET_H.py's identical function, see that
    file's design notes for the full derivation. Call ONCE per scan, before
    the main excitation pulse; emits hardware instructions directly (no
    return value). Switches TX to the LOW-POWER port for the duration and
    back to HIGH-POWER afterward -- returns immediately, before touching
    the TX port at all, if WetOn=0."""
    if int(P.WetOn) == 0:
        return

    WetOffsetHz = P.WetOffsetPPM * P.FrequencyBase

    Transmit1SelectPort(1, 0)
    Transmit1LPEnable(1, 1)
    Delay(10000)  # relay settle -- LP-port house convention, see design notes

    for amp, spoil in ((P.WetAmp1, P.WetSpoil1), (P.WetAmp2, P.WetSpoil2),
                        (P.WetAmp3, P.WetSpoil3), (P.WetAmp4, P.WetSpoil4)):
        shaped_pulse(P.WetPulseWidth, P.WetAngle, P.WetShape, amp,
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
# 69.3/161.0 deg) -- EXACT duplicate of LWG_PGSTE-WET_H.py's identical
# constants; see that file's design notes for the full derivation and
# Bloch-simulation cross-check. Deliberately NOT evaluated at import time
# -- see the SpinFlow-load-failure fix applied across this repo 10/08/2026.
_DEFAULT_WET_AMP1 = 0.0035    # 81.4 deg
_DEFAULT_WET_AMP2 = 0.00436   # 101.4 deg
_DEFAULT_WET_AMP3 = 0.00298   # 69.3 deg
_DEFAULT_WET_AMP4 = 0.00692   # 161.0 deg


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    DW = ReceiverFilter.dwell
    ReceiverTime = P.ReceiverPoints*DW

    # Matches run()'s actual per-scan timing exactly (derived by summing
    # every Delay()/duration argument in the per-scan block, term for
    # term, then algebraically simplified via symbolic expansion -- same
    # method used throughout this repo's timing audit, see xpulse-pulse-
    # programmes skill sec.4). Purely sequential (no 'with parallel:' in
    # the diffusion-encoding part), so no MAX-of-branches step needed here.
    # RD + 2*Tau survive in full (as in LWG_PGSTE-WET_H.py); P90/P180 each
    # contribute only HALF their width net (the other half cancels against
    # the D3a/D3b "rest of TAU" subtractions, exactly as in
    # LWG_PGSTE_H.py/LWG_PGSTE-WET_H.py, generalised here to P90!=P180).
    # The WET block adds its own per-scan time on TOP (it runs once, after
    # RD, before P1) -- identical formula to LWG_PGSTE-WET_H.py, zero when
    # WetOn=0.
    if int(P.WetOn) != 0:
        wet_time = (4*(P.WetPulseWidth + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD
                        + 2*P.WetRampTime + P.WetGradTime + P.WetGradSettle
                        + P.WetInterDelay + 2.0)
                    + 20000.0)
    else:
        wet_time = 0.0

    t_scanTime = (P.RecycleDelay + 2*P.Tau + 2*P.TXEnableTime
                  - P.P90/2.0 - P.P180/2.0
                  + ReceiverFilter.group_delay + ReceiverTime + 203.0
                  + wet_time)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Pulsed-Gradient Spin Echo (PGSE / Stejskal-Tanner) on the "
                "{H/F} channel, with a WET water-suppression module (4 "
                "selective pulses + decreasing spoiler gradients, "
                "targeting WetOffsetPPM) prepended to every scan: [WET] - "
                "hard-90 - TAU(gradient) - hard-180 - TAU(gradient) - "
                "echo. Used for measurement of diffusion constants where "
                "Delta stays well within T2 -- see LWG_PGSTE-WET_H.py for "
                "the T2&lt;&lt;T1 stimulated-echo alternative. Non-localised "
                "(no imaging/slice-selection). Set WetOn=0 to disable "
                "suppression.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,G1,WetOn,WetOffsetPPM"

    return basic

@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_PGSE-WET_H", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 59.7, ParameterTypes.Double, "H/F Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "H/F Freq Offset [Hz]")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "H/F TX Freq Offset [ppm]")

    ReceiverPoints = Parameter("NP", 32768, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 30, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "5000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans -- NumScans (and DS) should be a MULTIPLE OF 2 to complete the
    # EXORCYCLE-style phase cycle below.
    NumScans = Parameter("NS", 8, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 4, ParameterTypes.Int32,"Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 0.0, ParameterTypes.Double, "H/F Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 2000000, ParameterTypes.Int32, "Relaxation Delay [s] -- quiet wait BEFORE the WET module starts",
                              RD, min=100000, max=2000000000)

    # Hard pulses
    P90 = Parameter("P90", 10.0, ParameterTypes.Double, "H/F 90&#176; Pulse Width [&#956;s]")
    P180 = Parameter("P180", 20.0, ParameterTypes.Double, "H/F 180&#176; Pulse Width [&#956;s]")
    TXAmplitude = Parameter("RFA0", 0.4, ParameterTypes.Double,
                            "H/F TX Power [0.0&#8230;1.0]", RFA, min=0.0, max=1.0)

    # Gradients -- same G1 amplitude used for BOTH gradient pulses (see
    # design notes: DO NOT flip sign between them -- the 180 does the
    # inverting for you).
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
                     "TAU [&#956;s], P1-to-P2 and P2-to-acq (symmetric); "
                     "TE=2*TAU, Delta&#8776;TAU -- must exceed "
                     "PreGrad+2*RampTime+GradientOnTime+GradSettle+P90/2+P180/2")

    # ---- WET (Water suppression Enhanced through T1 effects) -------------
    # Standard 4-pulse WET module (Ogg, Kingsley &amp; Taylor 1994 -- see
    # design notes), prepended to every scan before P1. Set WetOn=0 to
    # disable entirely (zero added time/RF/gradient when off).
    WetOn = Parameter("WetOn", 1, ParameterTypes.Int32, "Enable WET water suppression [1=On(default), 0=Off]")
    WetOffsetPPM = Parameter("WetOffsetPPM", 4.7, ParameterTypes.Double, "*** WET target [ppm] rel. SF+O1 -- SET to your water peak ***")
    WetShape = Parameter("WetShape", "GAUSSIAN", ParameterTypes.String, "WET Pulse Shape [GAUSSIAN(default)/SINC]")
    # NOTE: attribute name deliberately does NOT end in "Phase" -- SpinFlow's
    # PhasesManager auto-discovers ANY Parameter whose attribute name ends
    # with "Phase" (matching P1Phase/P2Phase/RXPhase above) and
    # unconditionally tries to parse its value as a comma-separated integer
    # phase-cycle list, regardless of declared type. A first version of this
    # file named this WetPhase and crashed on real hardware with
    # "AttributeError: 'float' object has no attribute 'split'" inside
    # PhasesManager.__init__ -- see changelog.
    WetAngle = Parameter("WETANG", 0.0, ParameterTypes.Double, "WET Pulse Phase [&#176;] -- same for all 4, spoiled after each so phase doesn't matter")
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

    # Phases -- simple EXORCYCLE-style 2-step cycle (see design notes). The
    # WET module's own WetAngle (above) is fixed, not cycled.
    P1Phase = Parameter("PH1","0", ParameterTypes.String,"P1 (excitation) RF-Pulse Phase")
    P2Phase = Parameter("PH2","90,270", ParameterTypes.String,"P2 (refocus) RF-Pulse Phase [EXORCYCLE, 2-step]")
    RXPhase = Parameter("PHRX","0,180", ParameterTypes.String,"Acquisition Phase [2-step, matches PH2]")

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

    # Individual PhaseListContainer per phase Parameter, matching the
    # vendor's own confirmed-working pattern (PGSE_H.py/WET-PRE_H.py,
    # Default pps folder) -- NOT PhasesManager(P)/ph["PH1"]-style dict
    # indexing, which crashed on real hardware with "KeyError: 'PH1'"
    # (PhasesManager's returned dict is evidently NOT keyed by the raw
    # short code) -- see changelog. Incremented phase values are named
    # PH1/PH2/PHRX (matching the short codes exactly), since those are
    # what get passed straight into Channel1SetBasePhase()/Receiver1Phase().
    P1Phase = PhaseListContainer(P.P1Phase)
    P2Phase = PhaseListContainer(P.P2Phase)
    RXPhase = PhaseListContainer(P.RXPhase)

    points = P.ReceiverPoints
    times = np.arange(0, points*ReceiverFilter.dwell, ReceiverFilter.dwell) / 1.0e6

    GradWidth = 2*P.RampTime + P.GradientOnTime
    GradientSlewRate = abs(P.G1/P.RampTime)

    # "Rest of TAU" after the gradient/settle and pulse tails, on EACH side
    # -- same pattern as LWG_PGSTE_H.py's D3a/D3b, generalised here since
    # P1(90) and P2(180) have different widths (each side subtracts its
    # OWN adjacent pulse's half-width, not a shared P90/2 twice).
    D3a = P.Tau - P.PreGrad - GradWidth - P.GradSettle - (P.P90/2.) - (P.P180/2.) - 7
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

        P1Phase.Reset()
        P2Phase.Reset()
        RXPhase.Reset()

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
    rf_on_time = (P.P90 + P.P180) + wet_rf_on
    grad_on_time = 2*GradWidth + wet_grad_on
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                P1Phase.Reset()
                P2Phase.Reset()
                RXPhase.Reset()

            PH1 = P1Phase.Inc()
            PH2 = P2Phase.Inc()
            PHRX = RXPhase.Inc()

            Channel1SetBasePhase(10, PH1)
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

            # First TAU: PreGrad, gradient, settle, rest
            safe_delay(P.PreGrad-2-(P.P90/2.), "pre-grad-1 wait", comms)
            Gradient3SlewRate(1.0, GradientSlewRate)
            Gradient3(P.RampTime, P.G1)
            Delay(P.GradientOnTime)
            Gradient3(P.RampTime, 0)
            Gradient3(1,0) # Turn Gradient off explicitly for safety
            Delay(P.GradSettle)
            safe_delay(D3a, "first-TAU rest wait", comms)

            # P2 -- refocusing 180 (SAME gradient sign follows -- see design notes)
            Channel1SetBasePhase(3, PH2)
            Transmit1BlankingOn(1)
            Delay(P.TXEnableTime)
            Transmit1(P.P180)
            Transmit1BlankingOff(1)

            # Second TAU: PreGrad, gradient (SAME sign as above), settle,
            # rest (minus receiver dead-time/group-delay, paid back
            # explicitly at ACQU).
            safe_delay(P.PreGrad-2-(P.P180/2.), "pre-grad-2 wait", comms)
            Gradient3SlewRate(1.0, GradientSlewRate)
            Gradient3(P.RampTime, P.G1)
            Delay(P.GradientOnTime)
            Gradient3(P.RampTime, 0)
            Gradient3(1,0) # Turn Gradient off explicitly for safety
            Delay(P.GradSettle)
            safe_delay(D3b-ReceiverFilter.dead_time-2+ReceiverFilter.group_delay, "second-TAU rest wait", comms)

            # ACQU
            Channel1SetBasePhase(1,0)
            Receiver1Phase(1, PHRX)
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
# 1. Claude - 10/08/26 - Initial version. Built at your request alongside
#    LWG_PGSTE-WET_H.py: a WET-enabled two-pulse (90-180) PGSE diffusion
#    sequence, structurally a near-verbatim port of the vendor's PGSE_H.py
#    (Default pps folder) -- same GradientMatrix/Gradient3/ramp-and-settle
#    pattern and SAME-sign-both-gradients convention, confirmed working on
#    this instrument -- with this repo's house-convention helpers
#    (safe_delay()/mains_lock_trigger()/estimate_duty_cycles()/Probe
#    selector/PhasesManager) and an exact term-by-term (not the vendor's
#    own approximate) time_calculation() derivation added on top, matching
#    LWG_PGSTE-WET_H.py's approach. WET module (wet_suppression()/
#    generate_shape()/shaped_pulse()/R_GAUSSIAN_90/_DEFAULT_WET_AMP1-4) is
#    an EXACT duplicate of LWG_PGSTE-WET_H.py's -- see that file's design
#    notes/changelog for the full WET literature reference, R_GAUSSIAN_90
#    Bloch-simulation derivation, and WetAmp1-4 calibration-formula +
#    Bloch cross-check. WetOn=0 gives a plain (unsuppressed) PGSE
#    equivalent to a house-convention-ised port of vendor PGSE_H.py, with
#    zero added time/RF/gradient. Verified via the mock harness, cross-
#    checked symbolically (sympy) against LWG_PGSTE-WET_H.py's own
#    derivation method for consistency.
# 2. Claude - 14/08/26 - FIXED a real crash on real hardware: mirrors
#    LWG_PGSTE-WET_H.py's identical fix -- the WET pulse-phase Parameter
#    was originally named WetPhase (attribute AND short code), which
#    collided with SpinFlow's PhasesManager auto-discovery of any
#    Parameter attribute ending in "Phase" (matching P1Phase/P2Phase/
#    RXPhase, genuinely comma-separated lists), crashing with
#    "AttributeError: 'float' object has no attribute 'split'" inside
#    PhasesManager.__init__ before any WET-specific code even ran (you
#    supplied the full instrument traceback, which identified the exact
#    mechanism). Renamed to WetAngle (short code WETANG). Also hardened
#    the mock hardware harness's PhasesManager to reproduce this exact
#    check locally, so this class of bug is caught before deployment.
# 3. Claude - 14/08/26 - FIXED a second, distinct crash on real hardware:
#    mirrors LWG_PGSTE-WET_H.py's identical fix -- "KeyError: 'PH1'" at
#    Channel1SetBasePhase(10, ph["PH1"]), where ph = Phases.Incd() and
#    Phases = PhasesManager(P). PhasesManager's returned dict is evidently
#    NOT keyed by the raw short code "PH1"/"PH2"/"PHRX" the way this file
#    (and this repo generally) assumed. Switched to individual
#    PhaseListContainer objects per phase Parameter (P1Phase/P2Phase/
#    RXPhase), matching the vendor's own confirmed-working pattern
#    (PGSE_H.py/WET-PRE_H.py, Default pps folder) instead of going through
#    PhasesManager at all. Incremented per-scan phase values are named
#    PH1/PH2/PHRX -- matching the short codes exactly. Verified via the
#    mock harness (no exception, full scan loop completes).
#
# -----------------------------------------------------------------------------
