#-------------------------------------------------------------------------------
# Name:        LWG_Selective-Echo_X.py
# Purpose:     Frequency-selective spin echo, NO IMAGING / NO GRADIENTS --
#              shaped-90 - TAU - shaped-180 - TAU - acquire on a single
#              targeted resonance. This is a CALIBRATION tool: use it to
#              find the correct P90sh/P180sh/RFAsh0/RFAsh1/PulseOffset
#              before trusting LWG_1D-Image-Echo-Selective_X.py's selective
#              image. X-Pulse Broadband Benchtop NMR Spectrometer (X
#              channel) -- mechanical port of LWG_Selective-Echo_H.py, see
#              X-CHANNEL CALIBRATION note below.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     04/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     2.0
#
# X-CHANNEL CALIBRATION -- IMPORTANT: this is a mechanical H->X port (all
# Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1). SF
# now defaults to 15.01 MHz as a placeholder -- P90sh/P180sh/RFAsh0/RFAsh1/
# RA below are still the 1H-calibrated values from LWG_Selective-Echo_H.py
# and are almost certainly WRONG for whatever nucleus you actually put on
# the X channel. This file exists precisely so you can (re)calibrate those
# values for your target nucleus before trusting LWG_1D-Image-Echo-
# Selective_X.py.
#
# Design notes:
#  - This is the RF-only "core" of LWG_1D-Image-Echo-Selective_H.py with
#    every gradient (encode, read, apply_gradient, GradientMatrix, axis
#    selection) removed -- just the two shaped pulses, the two TAU delays,
#    and acquisition, all fully sequential (no 'with parallel:' needed at
#    all, since there's nothing left to run concurrently). Shares the exact
#    same generate_shape()/shaped_pulse() machinery (on-the-fly EBURP1/
#    REBURP/GAUSSIAN/SINC synthesis, LOW-POWER port, dedicated PulseOffset
#    Parameter, mains-lock-off default) as the imaging version, so whatever
#    you calibrate here (P90sh, P180sh, RFAsh0, RFAsh1, PulseOffset) can be
#    copied straight across.
#  - HOW TO CALIBRATE: acquire a plain (non-selective, hard-pulse) spectrum
#    of your sample first, so you know where your target resonance sits.
#    Then here:
#      1) Set PulseOffset to the Hz offset of your target peak (relative to
#         SF+O1 -- see the PulseOffset Parameter description).
#      2) With RefocusShape temporarily set to 'GAUSSIAN' or a very low
#         RFAsh1 (i.e. effectively no refocusing), sweep RFAsh0 (or P90sh)
#         and watch the acquired FID/spectrum amplitude to find the
#         EXCITATION 90 condition (maximum signal for a simple pulse-
#         acquire; for a true nutation curve go past 90 to see the signal
#         fall back through zero at 180, which pins down the calibration
#         more precisely than a single maximum).
#      3) Once P90sh/RFAsh0 (or PulseOffset) are set, sweep RFAsh1 (or
#         P180sh) with the excitation pulse fixed, and look for the
#         REFOCUSED ECHO amplitude to maximise -- that is your 180
#         condition for the refocusing shape.
#      4) Carry the calibrated PulseOffset/P90sh/P180sh/RFAsh0/RFAsh1
#         across to LWG_1D-Image-Echo-Selective_X.py's identically-named
#         Parameters.
#    See also the markdown calibration-procedure document (dB-to-relative-
#    power conversion + a more detailed nutation-calibration writeup) that
#    should be sitting alongside this file.
#  - Selective pulses are long (ms, not us) -- see the design notes in
#    LWG_1D-Image-Echo-Selective_X.py for why, and for the on-the-fly
#    shape-synthesis explanation (identical here).
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
            self.comms.log("seqTime LWG_Selective-Echo_X Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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

    t_scanTime = (P.P90sh + P.P180sh + P.Tau*2 + points*DW)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Frequency-selective spin echo (shaped-90-TAU-shaped-180-"
                "TAU-acquire), NO gradients/imaging, on the {H/F} channel. "
                "Calibration tool for LWG_1D-Image-Echo-Selective_H.py: use "
                "PulseOffset to target a resonance and sweep P90sh/P180sh/"
                "RFAsh0/RFAsh1 to find the true 90/180 condition before "
                "combining the selective pulses with imaging gradients.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,P90sh,P180sh,RFAsh0,RFAsh1,PulseOffset"

    return basic

def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical to LWG_1D-Image-Echo-Selective_H.py's version (see
    that file for the full docstring); duplicated here so this file stays
    fully self-contained.

    shape_name: 'GAUSSIAN', 'SINC', or 'EBURP1'/'REBURP'/'BURP' (needs
    burp_coeffs_A -- see ExBurpCoeffsA/B, RefocusBurpCoeffsA/B Parameters)."""
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

    elif shape_name in ('EBURP1', 'REBURP', 'BURP'):
        if not burp_coeffs_A:
            raise ValueError(
                "shape '{0}' requested but no Fourier A-coefficients were "
                "supplied. EBURP1/REBURP are defined in Geen, H. &amp; "
                "Freeman, R., J. Magn. Reson. 93, 93-141 (1991) as a "
                "truncated Fourier series -- enter coefficients as "
                "ExBurpCoeffsA/B or RefocusBurpCoeffsA/B, or use 'GAUSSIAN' / "
                "'SINC' instead, which need no coefficients."
                .format(shape_name))
        tau = np.linspace(0.0, 1.0, n_steps, endpoint=False)
        b1 = np.full(n_steps, burp_coeffs_A[0], dtype=float)
        for n in range(1, max(len(burp_coeffs_A), len(burp_coeffs_B))):
            An = burp_coeffs_A[n] if n < len(burp_coeffs_A) else 0.0
            Bn = burp_coeffs_B[n] if n < len(burp_coeffs_B) else 0.0
            b1 = b1 + An*np.cos(2*np.pi*n*tau) + Bn*np.sin(2*np.pi*n*tau)
        amp = np.abs(b1)
        phase = np.where(b1 < 0, 180.0, 0.0)

    else:
        raise ValueError("Unknown on-the-fly shape '{0}'. Known: GAUSSIAN, "
                          "SINC, EBURP1/REBURP/BURP (needs BurpCoeffsA/B)."
                          .format(shape_name))

    peak = np.max(np.abs(amp))
    if peak > 0:
        amp = amp/peak
    return amp, phase


SHAPED_PULSE_FIXED_OVERHEAD = 14
# See LWG_1D-Image-Echo-Selective_H.py's shaped_pulse() for the exact
# instruction-cost derivation of this constant (identical here).

def shaped_pulse(duration, phase, shape, amplitude, txenabletime,
                  base_frequency, pulse_offset,
                  burp_coeffs_A=None, burp_coeffs_B=None):
    """Apply a chemical-shift-selective shaped RF pulse, synthesised on the
    fly -- identical mechanism/parameters to LWG_1D-Image-Echo-
    Selective_H.py's shaped_pulse(); see that file for the full docstring."""
    n_steps = int(round(duration))
    if n_steps < 4:
        raise ValueError("shaped_pulse duration too short for step-wise "
                          "synthesis (need >=4 us so the shape has enough "
                          "points to be meaningful): got {0} us"
                          .format(duration))
    amp_profile, phase_profile = generate_shape(shape, n_steps, burp_coeffs_A, burp_coeffs_B)

    Channel2SetFrequency(1, base_frequency + pulse_offset)
    Transmit2SetScale(5, amplitude)
    Channel2SetBasePhase(5, phase)
    Transmit2BlankingOn(1)
    Delay(txenabletime)
    with parallel:
        with sequential:
            Transmit2(float(n_steps))
        with sequential:
            for a in amp_profile:
                Transmit2SetScale(1.0, amplitude*float(a))
        with sequential:
            for p in phase_profile:
                Channel2SetBasePhase(1.0, phase + float(p))
    Transmit2BlankingOff(1)
    Channel2SetFrequency(1, base_frequency)


def safe_delay(value, label, comms):
    """Delay() wrapper -- raises a clear, specific error instead of a
    silent/confusing failure if a computed delay would be negative."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU, or decrease P90sh/P180sh, then "
               "retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger (ExternalTrigger[n]()) before the
    first pulse event, per Pulse Sequence Programming User Manual 01-U-049
    section 3.7.19. OFF (UseMainsLock=0) by default -- see the identical
    function in LWG_1D-Image-Echo-Selective_H.py for the full reasoning."""
    if int(P.UseMainsLock) == 0:
        return
    triggers = {1: ExternalTrigger1, 2: ExternalTrigger2, 3: ExternalTrigger3}
    ch = int(P.MainsLockChannel)
    if ch not in triggers:
        comms.log("WARNING: MainsLockChannel={0} is not 1, 2, or 3 -- mains "
                  "lock trigger skipped.".format(ch))
        return
    triggers[ch]()


def parse_burp_coeffs(s):
    """Parse a comma-separated 'A0,A1,A2,...' string Parameter into a list
    of floats. Returns [] for an empty/whitespace-only string."""
    s = (s or "").strip()
    if not s:
        return []
    return [float(x) for x in s.split(",") if x.strip() != ""]


def db_to_relative_scale(dB, reference_relative_scale, convention='amplitude'):
    """Convert a dB-referenced RF power/amplitude figure (the units older
    NMR literature -- including Geen &amp; Freeman 1991 -- typically reports
    shaped-pulse power in) into an X-Pulse linear relative-scale value
    (0.0..1.0, the units Transmit2SetScale()/RFAsh0/RFAsh1 use).

    dB is interpreted as an ATTENUATION relative to reference_relative_scale
    (i.e. dB=0 means "the same as the reference"; positive dB means LESS
    power/amplitude than the reference, matching standard attenuator-dB
    convention). reference_relative_scale is YOUR OWN calibrated relative-
    scale value for a known, real flip angle (e.g. RFA0 for your normal
    hard 90 pulse) -- this function cannot do anything useful without a
    real calibration point to anchor to; it does not know your probe/coil
    characteristics.

    convention:
      'amplitude' (default) -- dB defined on B1/voltage amplitude, the
                   usual convention for RF pulse "power" in NMR (flip angle
                   is linear in B1, so a 6 dB attenuation halves B1 and
                   quarters the deposited power). Conversion: 10**(-dB/20).
      'power'    -- dB defined on deposited power directly (as in many RF
                   component/attenuator datasheets, W or mW figures).
                   Conversion: 10**(-dB/10).
    CHECK WHICH CONVENTION THE PAPER (OR WHATEVER YOU'RE CONVERTING FROM)
    ACTUALLY USES before trusting a result -- getting this wrong is exactly
    a factor-of-2-in-dB (i.e. sqrt()) error in the resulting flip angle.
    """
    if convention == 'amplitude':
        divisor = 20.0
    elif convention == 'power':
        divisor = 10.0
    else:
        raise ValueError("convention must be 'amplitude' or 'power', got %r" % (convention,))
    return reference_relative_scale * (10.0 ** (-dB / divisor))


def db_to_lp_relative_scale(dB, hp_reference_relative_scale, lp_max_fraction, convention='amplitude'):
    """As db_to_relative_scale() above, but for the common situation here:
    your dB figure and reference calibration point are both referenced to
    the HIGH-power port, but the shaped pulse itself runs on the LOW-power
    port (see 'Switch to LOW-power TX port' in run()). This re-expresses
    the result as a scale relative to LP's OWN (smaller) full scale, i.e.
    what you'd actually put in RFAsh0/RFAsh1.

    lp_max_fraction: your measured LPMaxFraction (LP full scale as a
    fraction of HP full scale, e.g. 0.10 for ~10%).
    """
    if lp_max_fraction <= 0:
        raise ValueError("lp_max_fraction must be > 0, got %r" % (lp_max_fraction,))
    absolute_fraction_of_hp_max = db_to_relative_scale(dB, hp_reference_relative_scale, convention)
    return absolute_fraction_of_hp_max / lp_max_fraction


def shape_integration_factor(amp_profile):
    """Average/peak ratio of a synthesised shape envelope (generate_shape()
    always peak-normalises amp_profile to 1.0, so this is just the mean).
    This is the same quantity Bruker's shape tool reports as 'Integ.
    Factor' (a.k.a. Bp/B1): the ratio of a shape's AVERAGE B1 field to a
    rectangular pulse's B1 field of the same duration. It's needed because
    a shaped pulse's net rotation depends on the average amplitude over the
    pulse, not the peak -- e.g. a shape spending most of its duration at
    low amplitude needs a much higher PEAK amplitude than a rectangular
    pulse to deliver the same net rotation."""
    peak = np.max(np.abs(amp_profile))
    if peak <= 0:
        raise ValueError("shape_integration_factor: shape has zero peak amplitude")
    return float(np.mean(np.abs(amp_profile)) / peak)


def estimate_shape_relative_scale(hard_pulse_width, hard_pulse_relative_scale,
                                   shape_duration, target_rotation_deg,
                                   integration_factor, lp_max_fraction,
                                   power_adjust_dB=0.0):
    """Shaped-pulse power calculation, TRANSLATED from a working Bruker/
    TopSpin pulse programme you supplied (its cnst2/spw2 lines), adapted
    from Bruker's Watt/power-dB-referenced scale to X-Pulse's linear
    relative-AMPLITUDE scale (0.0..1.0, what Transmit2SetScale()/RFAsh0/
    RFAsh1 actually take).

    The Bruker macro computes (in dB, power convention):
        cnst2 = cnst11 - 20*log10((p1*totrot2)/(p2*90)) + 20*log10(integfac2) + cnst0
        spw2  = 10**(-cnst2/10)
    which reduces algebraically to a POWER ratio:
        spw2/plw1 = [(p1*totrot2)/(p2*90*integfac2)]**2 * 10**(-cnst0/10)
    Since power is proportional to amplitude^2, the equivalent AMPLITUDE
    ratio (which is what X-Pulse's relative scale actually is) is just the
    square root of that, i.e. no squaring/square-rooting needed at all if
    you work in amplitude terms directly -- which is what this function
    does:
        scale_needed = hard_pulse_relative_scale
                       * (hard_pulse_width * target_rotation_deg)
                       / (shape_duration * 90.0 * integration_factor)
                       * 10**(-power_adjust_dB / 20.0)

    hard_pulse_width/hard_pulse_relative_scale: your calibrated HARD 90
    pulse's width [us] and relative amplitude on the HIGH-power port (p1/
    plw1 in the Bruker macro).
    shape_duration: this shaped pulse's duration [us] (p2).
    target_rotation_deg: the rotation this shape is designed to produce
    (totrot2) -- 90 for E-BURP-1 (excitation), 180 for RE-BURP (refocusing).
    integration_factor: from shape_integration_factor() above -- computed
    directly from the ACTUAL synthesised shape here, not looked up from a
    table (integfac2 in the Bruker macro).
    power_adjust_dB: manual fine-tune (amplitude-convention dB), equivalent
    to the Bruker macro's cnst0 -- use this to correct any residual offset
    AFTER checking against a real nutation curve, not as a first guess.

    Result is referenced to the HIGH-power port (same port the hard pulse
    was calibrated on); divide by lp_max_fraction to re-express relative to
    the LOW-power port's own (smaller) full scale, which is what actually
    goes into RFAsh0/RFAsh1 -- done automatically below.
    """
    if shape_duration <= 0 or integration_factor <= 0 or lp_max_fraction <= 0:
        raise ValueError("estimate_shape_relative_scale: shape_duration, "
                          "integration_factor, and lp_max_fraction must all "
                          "be > 0 (got {0}, {1}, {2})".format(
                          shape_duration, integration_factor, lp_max_fraction))
    hp_equivalent_scale = (hard_pulse_relative_scale
                            * (hard_pulse_width * target_rotation_deg)
                            / (shape_duration * 90.0 * integration_factor)
                            * (10.0 ** (-power_adjust_dB / 20.0)))
    return hp_equivalent_scale / lp_max_fraction


def estimate_rf_duty_cycle(P, rf_on_time, comms):
    """Log estimated RF duty cycle for this scan and warn if it exceeds the
    (user-adjustable) MaxRFDuty parameter. No gradients in this sequence,
    so unlike the imaging versions there's no grad_on_time to report.

    IMPORTANT: MaxRFDuty is a conservative, editable placeholder (5%), NOT
    a vendor-confirmed rating -- confirm the real limit with Oxford
    Instruments.
    """
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR

    comms.log("Estimated RF duty cycle: {0:.3%} (limit {1:.1%})".format(rf_duty, P.MaxRFDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD, shorter P90sh/P180sh, or "
                  "lower RFAsh0/RFAsh1, or confirm with Oxford Instruments "
                  "that this is within the transmitter's rated duty cycle "
                  "before running unattended.".format(rf_duty, P.MaxRFDuty))


@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_Selective-Echo_X", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 15.01, ParameterTypes.Double, "X Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "X Freq Offset [Hz] -- ref. for pulse+acq; use PulseOffset to target the shaped pulses instead")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "X TX Freq Offset [ppm]")
    # THE clear, dedicated place to set the selective pulses' target
    # frequency -- see LWG_1D-Image-Echo-Selective_H.py's PulseOffset for
    # the full explanation. This is the parameter to SWEEP during
    # calibration to locate/tune the pulse's frequency response.
    PulseOffset = Parameter("PulseOffset", 0.0, ParameterTypes.Double, "Selective-Pulse Freq Offset [Hz] from SF+O1 (90/180 only, not acq) -- SWEEP to find target resonance")

    ReceiverPoints = Parameter("NP", 1024, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 34, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "100000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans
    NumScans = Parameter("NS", 4, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 0, ParameterTypes.Int32,"Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 100.0, ParameterTypes.Double, "X Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 500000, ParameterTypes.Int32, "Relaxation Delay [s]",
                              RD, min=100000, max=2000000000)

    # Shaped (chemical-shift-selective) pulses -- synthesised ON THE FLY
    # (see generate_shape()/shaped_pulse() above), no external shape file
    # needed. Default to the real literature pulses: EXCITATION = EBURP1
    # (Geen &amp; Freeman 1991, Table 2, nmax=8), REFOCUSING = REBURP (Table 8,
    # Np=256).
    ExcitationShape = Parameter("ExShape", "EBURP1", ParameterTypes.String, "Excitation Shape [EBURP1(default)/GAUSSIAN/SINC/BURP]")
    RefocusShape = Parameter("RefocusShape", "REBURP", ParameterTypes.String, "Refocusing Shape [REBURP(default)/GAUSSIAN/SINC/BURP]")
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00", ParameterTypes.String, "E-BURP-1 cosine coeffs A0..A8 (Geen&amp;Freeman'91 Tbl.2) -- used if ExShape=EBURP1/BURP")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01", ParameterTypes.String, "E-BURP-1 sine coeffs B0..B8 (Geen&amp;Freeman'91 Tbl.2) -- used if ExShape=EBURP1/BURP")
    RefocusBurpCoeffsA = Parameter("RefocusBurpCoeffsA", "0.49,-1.02,1.11,-1.57,0.83,-0.42,0.26,-0.16,0.10,-0.07,0.04,-0.03,0.01,-0.02,0.00,-0.01", ParameterTypes.String, "RE-BURP cosine coeffs A0..A15 (Geen&amp;Freeman'91 Tbl.8) -- used if RefocusShape=REBURP/BURP")
    RefocusBurpCoeffsB = Parameter("RefocusBurpCoeffsB", "", ParameterTypes.String, "RE-BURP sine coeffs (none published, purely real -- leave empty) -- used if RefocusShape=REBURP/BURP")
    P90sh = Parameter("P90sh", 5000.0, ParameterTypes.Double, "Shaped 90&#176; Width [&#956;s] (=shape resolution, 1pt/&#956;s) -- SWEEP to calibrate")
    P180sh = Parameter("P180sh", 10000.0, ParameterTypes.Double, "Shaped 180&#176; Width [&#956;s] (=shape resolution, 1pt/&#956;s) -- SWEEP to calibrate")
    # These are your calibration TARGETS -- start low and sweep upward
    # while watching the acquired signal (see design notes at top of file).
    TXAmplitude90 = Parameter("RFAsh0", 0.30, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0&#8230;1 of LP max] -- SWEEP to calibrate", RFA, min=0.0, max=1.0)
    TXAmplitude180 = Parameter("RFAsh1", 0.30, ParameterTypes.Double,
                               "Shaped 180&#176; TX Power [0&#8230;1 of LP max] -- SWEEP to calibrate", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", 0.10, ParameterTypes.Double, "LP Port Max as Fraction of HP [0&#8230;1] -- your measured value")

    # OPTIONAL power calculator -- OFF by default (PowerCalcMethod='manual'),
    # in which case RFAsh0/RFAsh1 above are used exactly as entered.
    #   'shape' (RECOMMENDED if you have a calibrated hard pulse): the
    #       physics-based calculation translated from your working Bruker
    #       pulse programme (see estimate_shape_relative_scale() above) --
    #       uses your hard-pulse calibration + each shape's OWN, actually-
    #       synthesised integration factor. No literature dB table needed.
    #   'db': legacy direct dB-attenuation route (see db_to_lp_relative_
    #       scale() above) -- use only if you already have a dB attenuation
    #       figure that is NOT derived from a hard-pulse ratio (e.g. quoted
    #       directly in a paper/datasheet).
    #   'manual' (default): RFAsh0/RFAsh1 above are used as entered.
    # Whichever method is used, VERIFY the result against a real nutation
    # curve (see design notes at top of file) before trusting it -- none of
    # these routes know your actual probe/coil/sample.
    PowerCalcMethod = Parameter("PowerCalcMethod", "manual", ParameterTypes.String, "RFAsh0/RFAsh1 source: manual(default)/shape(recommended)/db")
    RefAmplitude_HP = Parameter("RefAmplitude_HP", 0.40, ParameterTypes.Double, "Ref. hard-pulse rel. amplitude [0&#8230;1] on HP port, known flip angle -- for shape/db methods")
    HardPulseWidth = Parameter("P1Hard", 9.58, ParameterTypes.Double, "Ref. hard 90&#176; width [&#956;s] at RefAmplitude_HP on HP port -- for 'shape' method")
    ExcitationRotation = Parameter("ExRotation", 90.0, ParameterTypes.Double, "Target rotation of excitation shape [&#176;] (EBURP1=90) -- for 'shape' method")
    RefocusRotation = Parameter("RefocusRotation", 180.0, ParameterTypes.Double, "Target rotation of refocusing shape [&#176;] (REBURP=180) -- for 'shape' method")
    PowerAdjust_dB = Parameter("PowerAdjust_dB", 0.0, ParameterTypes.Double, "Manual fine-tune [dB] on top of 'shape' calc -- set AFTER nutation check")
    Excitation_dB = Parameter("Excitation_dB", 0.0, ParameterTypes.Double, "Excitation power as dB attenuation vs RefAmplitude_HP (+ve=less power) -- for 'db' method")
    Refocus_dB = Parameter("Refocus_dB", 0.0, ParameterTypes.Double, "Refocusing power as dB attenuation vs RefAmplitude_HP (+ve=less power) -- for 'db' method")
    DbConvention = Parameter("DbConvention", "amplitude", ParameterTypes.String, "dB convention: amplitude(B1,usual)/power -- check your source; for 'db' method")

    # Echo timing -- no gradients here, so Tau only needs to comfortably
    # clear the shaped pulses themselves (no PreGrad/GradSettle margin
    # needed, unlike the imaging version).
    Tau = Parameter("TAU", 20000, ParameterTypes.Int32, "Echo &#964; Delay [&#956;s] -- must exceed P90sh/2 + P180sh")

    # Mains-lock trigger -- OFF by default (see mains_lock_trigger()
    # docstring). Included here for consistency with the rest of the pp
    # family, even though this isn't an imaging sequence.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock Trigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rail (RF only -- no gradients in this sequence)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0&#8230;1]")

    # Phases (simple 2-step)
    PH1 = Parameter("PH1", "0,180", ParameterTypes.String, "Shaped 90&#176; Pulse Phase")
    PH2 = Parameter("PH2", "0", ParameterTypes.String, "Shaped 180&#176; Pulse Phase")
    PHRX = Parameter("PHRX", "0,180", ParameterTypes.String, "Acquisition Phase")

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

    # Total acquisition-window duration -- Receiver2()'s 'duration' argument
    # is the TOTAL window (points*dwell), NOT the per-point dwell.
    ReceiverTime = DW*(points+1)

    times = np.arange(0, points*DW, DW) / 1.0e6

    # Exact elapsed time of a shaped_pulse() call -- see
    # SHAPED_PULSE_FIXED_OVERHEAD above.
    ExcitationPulseWidth = P.P90sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD
    RefocusPulseWidth = P.P180sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD

    ExBurpA = parse_burp_coeffs(P.ExBurpCoeffsA)
    ExBurpB = parse_burp_coeffs(P.ExBurpCoeffsB)
    RefocusBurpA = parse_burp_coeffs(P.RefocusBurpCoeffsA)
    RefocusBurpB = parse_burp_coeffs(P.RefocusBurpCoeffsB)

    # ---- Optional power calculator override -----------------------------
    RFAsh0_effective = P.TXAmplitude90
    RFAsh1_effective = P.TXAmplitude180
    PowerMethod = str(P.PowerCalcMethod).strip().lower()

    if PowerMethod == 'shape':
        # Compute each shape's integration factor from the ACTUAL
        # synthesised envelope (same generate_shape() shaped_pulse() itself
        # calls at pulse time) -- no external table needed. This is a pure
        # numpy calculation, no hardware instructions, so calling it here
        # ahead of the scan loop is free and has no side effects.
        ExAmpProfile, _ = generate_shape(P.ExcitationShape, max(4, int(round(P.P90sh))), ExBurpA, ExBurpB)
        RefocusAmpProfile, _ = generate_shape(P.RefocusShape, max(4, int(round(P.P180sh))), RefocusBurpA, RefocusBurpB)
        ExIntegFactor = shape_integration_factor(ExAmpProfile)
        RefocusIntegFactor = shape_integration_factor(RefocusAmpProfile)
        RFAsh0_effective = estimate_shape_relative_scale(
            P.HardPulseWidth, P.RefAmplitude_HP, P.P90sh, P.ExcitationRotation,
            ExIntegFactor, P.LPMaxFraction, P.PowerAdjust_dB)
        RFAsh1_effective = estimate_shape_relative_scale(
            P.HardPulseWidth, P.RefAmplitude_HP, P.P180sh, P.RefocusRotation,
            RefocusIntegFactor, P.LPMaxFraction, P.PowerAdjust_dB)
        comms.log("PowerCalcMethod='shape': excitation integ.factor={0:.4f} "
                  "-> RFAsh0={1:.4f}; refocus integ.factor={2:.4f} -> "
                  "RFAsh1={3:.4f} (hard-pulse ref P1Hard={4}us @ "
                  "RefAmplitude_HP={5}, ExRotation={6}deg, RefocusRotation="
                  "{7}deg, LPMaxFraction={8:.2%}, PowerAdjust_dB={9:.2f}). "
                  "(Parameter-panel RFAsh0/RFAsh1 are ignored while "
                  "PowerCalcMethod='shape'.)"
                  .format(ExIntegFactor, RFAsh0_effective, RefocusIntegFactor,
                          RFAsh1_effective, P.HardPulseWidth, P.RefAmplitude_HP,
                          P.ExcitationRotation, P.RefocusRotation,
                          P.LPMaxFraction, P.PowerAdjust_dB))
        if RFAsh0_effective > 1.0 or RFAsh1_effective > 1.0:
            comms.log("WARNING: calculated relative scale exceeds 1.0 (LP "
                      "port cannot go higher than its own max) -- "
                      "RFAsh0_effective={0:.4f}, RFAsh1_effective={1:.4f}. "
                      "Check P1Hard/RefAmplitude_HP/ExRotation/RefocusRotation/"
                      "LPMaxFraction.".format(RFAsh0_effective, RFAsh1_effective))

    elif PowerMethod == 'db':
        RFAsh0_effective = db_to_lp_relative_scale(P.Excitation_dB, P.RefAmplitude_HP, P.LPMaxFraction, P.DbConvention)
        RFAsh1_effective = db_to_lp_relative_scale(P.Refocus_dB, P.RefAmplitude_HP, P.LPMaxFraction, P.DbConvention)
        comms.log("PowerCalcMethod='db': RFAsh0 overridden -- {0:.2f} dB "
                  "(convention='{1}') off RefAmplitude_HP={2} -> {3:.4f} "
                  "(LP-relative, LPMaxFraction={4:.2%}). RFAsh1 overridden "
                  "-- {5:.2f} dB -> {6:.4f}. (Parameter-panel RFAsh0/RFAsh1 "
                  "values are ignored while PowerCalcMethod='db'.)"
                  .format(P.Excitation_dB, P.DbConvention, P.RefAmplitude_HP,
                          RFAsh0_effective, P.LPMaxFraction,
                          P.Refocus_dB, RFAsh1_effective))
        if RFAsh0_effective > 1.0 or RFAsh1_effective > 1.0:
            comms.log("WARNING: dB-calculated relative scale exceeds 1.0 "
                      "(LP port cannot go higher than its own max) -- "
                      "RFAsh0_effective={0:.4f}, RFAsh1_effective={1:.4f}. "
                      "Check RefAmplitude_HP/Excitation_dB/Refocus_dB/"
                      "LPMaxFraction/DbConvention.".format(RFAsh0_effective, RFAsh1_effective))

    elif PowerMethod != 'manual':
        comms.log("WARNING: unrecognised PowerCalcMethod='{0}' -- falling "
                  "back to 'manual' (RFAsh0/RFAsh1 used as entered). Valid "
                  "values: 'manual', 'shape', 'db'.".format(P.PowerCalcMethod))

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

        # Switch to the LOW-POWER (LP) TX port for these shaped pulses --
        # per your measurement, LP tops out at ~LPMaxFraction (default 10%)
        # of the high-power port's full scale. This sequence's RF is 100%
        # shaped pulses, so LP is engaged once here, for the whole
        # experiment.
        Transmit2SelectPort(1,0)
        Transmit2LPEnable(1,1)
        Delay(10000) # relay settle -- see LWG_1D-Image-Echo-Selective_H.py
                      # for the same caution/value.

        Receiver2Preamp(128, P.ReceiverAttenuation)
        Receiver2Filter(200, ReceiverFilter)

        Phases.Reset()

    comms.log("LP port engaged for shaped pulses -- LPMaxFraction={0:.2%} "
              "of HP full scale. Pulse offset: PulseOffset={1} Hz "
              "(applied on top of SF+O1={2:.6f} MHz)."
              .format(P.LPMaxFraction, P.PulseOffset, Frequency))

    # ---- Duty-cycle warning (RF only, no gradients here) -----------------
    rf_on_time = ExcitationPulseWidth + RefocusPulseWidth
    estimate_rf_duty_cycle(P, rf_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver2FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            # Shaped 90 (excitation) pulse
            shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, RFAsh0_effective, P.TXEnableTime, Frequency, P.PulseOffset, ExBurpA, ExBurpB)
            # TAU
            safe_delay(P.Tau - ExcitationPulseWidth/2 - RefocusPulseWidth/2.0,
                       "first-TAU wait", comms)
            # Shaped 180 (refocusing) pulse
            shaped_pulse(P.P180sh, ph["PH2"], P.RefocusShape, RFAsh1_effective, P.TXEnableTime, Frequency, P.PulseOffset, RefocusBurpA, RefocusBurpB)
            # TAU -- ReceiverFilter.dead_time is reserved out of this wait
            # and paid back explicitly (as its own Delay, right before
            # Receiver2 below) rather than left out entirely.
            safe_delay(P.Tau - RefocusPulseWidth/2.0 - ReceiverFilter.dead_time,
                       "second-TAU wait", comms)
            # ACQU
            Channel2SetBasePhase(P.TXEnableTime,0)
            Receiver2Phase(P.TXEnableTime, ph["PHRX"])
            # Dead1: probe ring-down time.
            Delay(P.Dead1)
            # ReceiverFilter.dead_time: digital-filter settling time.
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
# 1. Claude - 04/08/26 - Initial version. Stripped-down, RF-only "core" of
#    LWG_1D-Image-Echo-Selective_H.py v4.0 -- same generate_shape()/
#    shaped_pulse() (on-the-fly EBURP1/REBURP/GAUSSIAN/SINC synthesis, LOW-
#    power port, dedicated PulseOffset Parameter, mains-lock-off default),
#    minus every gradient (no encode/read gradients, no GradientMatrix, no
#    axis selection) -- purely sequential (no 'with parallel:' needed, since
#    there's nothing to run concurrently with the RF). Built at your
#    request for a calibration tool: find the true P90sh/P180sh/RFAsh0/
#    RFAsh1/PulseOffset condition on a real spectrum via nutation before
#    trusting the imaging version's selective pulses.
# 2. Claude - 04/08/26 - dB POWER CALCULATOR: added db_to_relative_scale()/
#    db_to_lp_relative_scale() and an opt-in UseDbCalc/RefAmplitude_HP/
#    Excitation_dB/Refocus_dB/DbConvention Parameter group, per your
#    request to bridge Geen &amp; Freeman's dB-based power figures to the
#    X-Pulse's linear relative-power scale. I do NOT have the paper's
#    specific dB table (you gave me the Fourier coefficient tables, not a
#    power table), so this is general-purpose machinery for YOUR OWN dB
#    figure + YOUR OWN calibrated hard-pulse reference point -- it does not
#    assume or embed any particular literature number. OFF by default
#    (UseDbCalc=0): RFAsh0/RFAsh1 are used as entered, exactly as before.
#    See the accompanying markdown doc for the full calibration procedure
#    (this dB route AND the more fundamentally reliable nutation-curve
#    route, which doesn't need the paper's numbers at all).
# 3. Claude - 04/08/26 - SHAPE-INTEGRAL POWER CALCULATOR: added
#    shape_integration_factor()/estimate_shape_relative_scale(), TRANSLATED
#    from a working Bruker/TopSpin pulse programme you supplied (its
#    cnst11/cnst2/spw2 shaped-pulse power calculation), re-derived for
#    X-Pulse's linear relative-AMPLITUDE scale instead of Bruker's Watt/
#    power-dB scale (see estimate_shape_relative_scale()'s docstring for
#    the full algebra). Exposed as PowerCalcMethod='shape' (new preferred
#    option alongside the existing 'manual'/'db' methods -- UseDbCalc
#    renamed/generalised to PowerCalcMethod, a 3-way string selector, still
#    defaulting to 'manual' so nothing changes unless you opt in). Needs
#    only your existing hard-pulse calibration (P1Hard/RefAmplitude_HP) and
#    each shape's target rotation (ExRotation=90/RefocusRotation=180 by
#    default, matching E-BURP-1/RE-BURP) -- the shape's own integration
#    factor is computed directly from the SAME synthesised envelope
#    shaped_pulse() uses to drive the hardware, not looked up from a table,
#    so this needs none of Geen &amp; Freeman's literature dB figures. Kept the
#    legacy 'db' route (previously UseDbCalc=1) available for anyone who
#    already has a dB attenuation figure not derived from a hard-pulse
#    ratio. PowerAdjust_dB (equivalent to the Bruker macro's cnst0) is
#    available for post-hoc fine-tuning after nutation verification. See
#    the updated markdown doc for the full worked explanation.
#
# -----------------------------------------------------------------------------
