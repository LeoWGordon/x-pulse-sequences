#-------------------------------------------------------------------------------
# Name:        LWG_Slice-Selective-Echo_H.py
# Purpose:     TRUE spatial SLICE-selective spin echo (shaped E-BURP-1
#              excitation + shaped RE-BURP refocusing, gradient ON during
#              BOTH pulses), on the Oxford Instruments X-Pulse Broadband
#              Benchtop NMR Spectrometer (1H/19F channel). Axis-selectable
#              (x/y/z). Acquires a localised echo from one slab/slice -- no
#              frequency-encoded readout/imaging in this file.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     06/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.0
#
# Design notes -- READ THIS:
#
#  - CHEMICAL-SHIFT-SELECTIVE vs SLICE-SELECTIVE: as with
#    LWG_Slice-Selective-PulseAcquire_H.py, this is the genuinely spatial
#    counterpart of LWG_Selective-Echo_H.py (which leaves the gradient OFF
#    during its shaped pulses -- chemical-shift selectivity only). Here the
#    slice-select gradient (GradAxis/G1/FPX-Y-Z) is ON, at a constant
#    plateau, for the ENTIRE duration of BOTH shaped pulses.
#  - MATCHING THE EXCITATION AND REFOCUSING SLICES -- this is the part you
#    specifically flagged. E-BURP-1 (90) and RE-BURP (180) are DIFFERENT
#    pulse shapes with DIFFERENT bandwidth-time products (R = bandwidth x
#    duration; see R_EBURP1_90/R_REBURP_180 below). If you ran them both at
#    the SAME gradient strength AND the same duration, they would excite/
#    refocus DIFFERENTLY SIZED slices, and your echo would only reflect the
#    (narrower) OVERLAP region -- signal loss and a slice profile that is
#    not what either pulse's own width would suggest. The fix used here:
#    both pulses run on the SAME slice-select gradient strength (G1, same
#    axis), and the REFOCUSING pulse duration (P180sh) is derived from the
#    excitation duration (P90sh) via the ratio of the two R-values, so
#    bandwidth_ex/G1 == bandwidth_ref/G1, i.e. equal slice width, by
#    construction. This is automatic when AutoMatchSlice=1 (the default) --
#    P180sh is IGNORED and recomputed from P90sh every run (and logged, so
#    you can see the number). Set AutoMatchSlice=0 to take manual control of
#    P180sh (e.g. to deliberately probe a mismatched-slice condition) -- in
#    that case both slice widths are still calculated and logged, and you
#    get an explicit WARNING if they differ by more than a couple of
#    percent.
#  - WHY NO REPHASE LOBE ON THE 180: RE-BURP is a temporally SYMMETRIC pulse
#    (confirmed in the E-BURP-1/RE-BURP reconstruction sanity-check noted in
#    LWG_1D-Image-Echo-Selective_H.py's changelog item 4d -- "RE-BURP
#    reconstructs as expected for a refocusing pulse (symmetric about the
#    pulse centre)"). A symmetric 180 pulse is inherently self-refocusing
#    with respect to its OWN slice-select gradient (the dephasing picked up
#    in the first half is undone by the inverting second half) -- this is
#    precisely why RE-BURP-style pulses are the standard choice for
#    gradient-echo/slice-selective refocusing in MRI, and matches your own
#    working Bruker/TopSpin reference pulse programme (slice-sel_echo_TEST),
#    which likewise applies the SAME gradient (gpz1/gpnam1) under both the
#    Eburp2 excitation and the Reburp refocusing pulse with no separate
#    rephase gradient shown for either. The excitation pulse (E-BURP-1) IS
#    asymmetric, hence it alone gets an explicit rephase lobe (see
#    LWG_Slice-Selective-PulseAcquire_H.py's design notes) right after it,
#    BEFORE the TAU delay to the refocusing pulse.
#  - CRUSHER GRADIENTS: added around the refocusing pulse (same polarity,
#    equal area on each side -- the standard convention: this cancels any
#    magnetisation NOT correctly inverted by an imperfect real-world 180
#    (which would otherwise show up as an FID-like artefact riding on top of
#    the true echo), while leaving the correctly-refocused echo pathway
#    unaffected, since the 180 inverts the sign of dephasing it experiences,
#    so the two equal, same-sign crusher lobes cancel for that pathway only).
#    CrusherOn=1 by default. They run on the SAME axis as the slice-select
#    gradient for simplicity (this X-Pulse pp family has not demonstrated
#    simultaneous multi-axis gradients elsewhere) -- a small trade-off
#    (same-axis crushers can, in principle, interact with the slice profile
#    if set too strong; keep CrusherStrength modest, e.g. <=0.5, relative to
#    G1). Combined with an 8-step CYCLOPS+EXORCYCLE-style phase cycle (see
#    Phases below) for additional artefact suppression.
#  - EXCITED SLICE WIDTH: see LWG_Slice-Selective-PulseAcquire_H.py's design
#    notes for the full Bloch-simulation method behind R_EBURP1_90; the same
#    method (calibrated to an exact on-resonance 180, not 90, then FWHM of
#    the '-My recovered' refocusing profile) gives R_REBURP_180 = 5.96 for
#    RE-BURP, used as 5.9 below.
#  - GRADIENT CALIBRATION / PULSE CENTRE FREQUENCY: identical
#    GradCal_HzPerCm (uncalibrated placeholder) and SlicePPM (prominent,
#    near the top of Parameters) mechanism as
#    LWG_Slice-Selective-PulseAcquire_H.py -- see that file's design notes.
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
            self.comms.log("seqTime LWG_Slice-Selective-Echo_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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

# ---- Numerically-derived BURP bandwidth-time products --------------------
# See LWG_Slice-Selective-PulseAcquire_H.py's R_EBURP1_90 docstring for the
# full Bloch-simulation method (Rodrigues'-formula rotation, on-resonance
# flip angle calibrated EXACTLY via numerical search, then FWHM of the
# excitation/refocusing profile vs offset). R_REBURP_180 uses the same
# method but starting magnetisation M=(0,1,0), calibrated to an exact
# on-resonance 180 (i.e. -My fully recovered), and measuring the FWHM of the
# '-My recovered' profile vs offset: simulated R = 5.96, used as 5.9 below.
# Both are duration-independent properties of the pulse SHAPE alone.
R_EBURP1_90 = 4.9
R_REBURP_180 = 5.9

def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell

    t_scanTime = (P.P90sh + P.P180sh + P.Tau*2 + points*DW)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Spatially SLICE-selective spin echo (shaped E-BURP-1 "
                "excitation + shaped RE-BURP refocusing, BOTH with the "
                "slice-select gradient ON during the pulse), on the {H/F} "
                "channel. Excitation and refocusing slice widths are "
                "AUTOMATICALLY MATCHED (AutoMatchSlice=1 default) via the "
                "two shapes' bandwidth-time-product ratio. Crusher "
                "gradients flank the refocusing pulse. SlicePPM sets where "
                "the slice sits (ppm); ExcitedWidth for both pulses is "
                "calculated and logged/recorded every run.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,SlicePPM,TAU,P90sh,P180sh,AutoMatchSlice,RFAsh0,RFAsh1,G1,GradAxis"

    return basic

def get_gradient_functions(axisstr):
    """Map axis letter to gradient functions (same convention as the imaging
    pp family: x/y/z -> Gradient1/2/3)."""
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

def get_axis_fp(P, axisstr):
    """Per-axis calibration scaler (FPX/FPY/FPZ) matching GradAxis -- see
    LWG_Slice-Selective-PulseAcquire_H.py's identical helper."""
    if axisstr == 'x':
        return P.XGradNorm
    elif axisstr == 'y':
        return P.YGradNorm
    elif axisstr == 'z':
        return P.ZGradNorm
    else:
        exit()

def apply_gradient(axis, delta, ramp_time, gradient_value, settle_time, pre_grad):
    """Simple trapezoidal gradient lobe: ramp up, hold for delta, ramp down,
    settle. Used here for the crusher gradients (NOT for the slice-select
    gradients that run synchronised with the shaped pulses, which are
    written out explicitly in run() so their timing can be interleaved with
    shaped_pulse() inside a 'with parallel:' block) -- identical to the
    apply_gradient() already used in LWG_1D-Image-Echo_H.py /
    LWG_1D-Image-UTE_H.py."""
    Delay(pre_grad)
    gradient_slew_rate = abs(gradient_value / ramp_time)
    slew_rate_fn, gradient_fn = get_gradient_functions(axis)
    slew_rate_fn(1.0, gradient_slew_rate)
    gradient_fn(ramp_time, gradient_value)
    Delay(delta)
    gradient_fn(ramp_time, 0)
    gradient_fn(1.0, 0)
    Delay(settle_time)

def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical machinery to LWG_Selective-Echo_H.py /
    LWG_Slice-Selective-PulseAcquire_H.py; duplicated here so this file
    stays fully self-contained."""
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
                "supplied (BurpCoeffsA is empty). EBURP1/REBURP are "
                "defined in Geen, H. &amp; Freeman, R., J. Magn. Reson. 93, "
                "93-141 (1991) as a truncated Fourier series -- enter "
                "coefficients as ExBurpCoeffsA/B or RefocusBurpCoeffsA/B, or "
                "use 'GAUSSIAN' / 'SINC' instead, which need no "
                "coefficients.".format(shape_name))
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
# Exact fixed (non-parallel) instruction cost of shaped_pulse(), in us -- see
# LWG_1D-Image-Echo-Selective_H.py's identical constant for the full
# instruction-by-instruction derivation.

def shaped_pulse(duration, phase, shape, amplitude, txenabletime,
                  base_frequency, pulse_offset,
                  burp_coeffs_A=None, burp_coeffs_B=None):
    """Apply a shaped RF pulse, synthesised on the fly -- identical mechanism
    to LWG_Selective-Echo_H.py's shaped_pulse(). Deliberately does NOT touch
    any gradient -- the slice-select gradients are driven by SEPARATE
    parallel branches in run()."""
    n_steps = int(round(duration))
    if n_steps < 4:
        raise ValueError("shaped_pulse duration too short for step-wise "
                          "synthesis (need >=4 us so the shape has enough "
                          "points to be meaningful): got {0} us"
                          .format(duration))
    amp_profile, phase_profile = generate_shape(shape, n_steps, burp_coeffs_A, burp_coeffs_B)

    Channel1SetFrequency(1, base_frequency + pulse_offset)
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


def safe_delay(value, label, comms):
    """Delay() wrapper -- raises a clear, specific error instead of a
    silent/confusing failure if a computed delay would be negative."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU and/or PreGrad margins, or "
               "decrease P90sh/P180sh/RephaseFraction/CrusherTime, then "
               "retry.").format(label, value)
        comms.log(msg)
        raise ValueError(msg)
    Delay(value)


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger -- see the identical function in
    LWG_1D-Image-Echo_H.py for the full reasoning. OFF by default."""
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
    """See LWG_Selective-Echo_H.py's identical function for the full
    docstring."""
    if convention == 'amplitude':
        divisor = 20.0
    elif convention == 'power':
        divisor = 10.0
    else:
        raise ValueError("convention must be 'amplitude' or 'power', got %r" % (convention,))
    return reference_relative_scale * (10.0 ** (-dB / divisor))


def db_to_lp_relative_scale(dB, hp_reference_relative_scale, lp_max_fraction, convention='amplitude'):
    """See LWG_Selective-Echo_H.py's identical function for the full
    docstring."""
    if lp_max_fraction <= 0:
        raise ValueError("lp_max_fraction must be > 0, got %r" % (lp_max_fraction,))
    absolute_fraction_of_hp_max = db_to_relative_scale(dB, hp_reference_relative_scale, convention)
    return absolute_fraction_of_hp_max / lp_max_fraction


def shape_integration_factor(amp_profile):
    """Average/peak ratio of a synthesised shape envelope -- see
    LWG_Selective-Echo_H.py's identical function for the full docstring."""
    peak = np.max(np.abs(amp_profile))
    if peak <= 0:
        raise ValueError("shape_integration_factor: shape has zero peak amplitude")
    return float(np.mean(np.abs(amp_profile)) / peak)


def estimate_shape_relative_scale(hard_pulse_width, hard_pulse_relative_scale,
                                   shape_duration, target_rotation_deg,
                                   integration_factor, lp_max_fraction,
                                   power_adjust_dB=0.0):
    """See LWG_Selective-Echo_H.py's identical function for the full
    docstring (translated from your Bruker cnst2/spw2 macro)."""
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


def estimate_duty_cycles(P, rf_on_time, grad_on_time, comms):
    """Log estimated RF and gradient duty cycles -- see
    LWG_1D-Image-Echo_H.py's identical function for the full docstring.
    IMPORTANT: MaxRFDuty/MaxGradDuty are conservative placeholders, NOT a
    vendor-confirmed rating."""
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR
    grad_duty = grad_on_time / TR

    comms.log("Estimated duty cycle -- RF: {0:.3%} (limit {1:.1%}), "
              "Gradient [{2}-axis]: {3:.3%} (limit {4:.1%})".format(
              rf_duty, P.MaxRFDuty, P.Axis, grad_duty, P.MaxGradDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD or shorter P90sh/P180sh, "
                  "or confirm with Oxford Instruments that this is within "
                  "the transmitter's rated duty cycle before running "
                  "unattended.".format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} on the "
                  "{1}-axis coil exceeds MaxGradDuty ({2:.2%}). Consider a "
                  "longer RD, or confirm the rated duty cycle for the "
                  "X-Pulse gradient amplifier/coil with Oxford Instruments "
                  "-- this limit is a conservative placeholder, not a "
                  "vendor spec.".format(grad_duty, P.Axis, P.MaxGradDuty))


def compute_excited_width(P, comms, pulse_width_us, R_value, shape_label):
    """Compute the excited slice bandwidth [Hz] and width [cm] -- identical
    to LWG_Slice-Selective-PulseAcquire_H.py's compute_excited_width().
    Returns (bandwidth_hz, width_cm, grad_total)."""
    axis_fp = get_axis_fp(P, P.Axis)
    grad_total = P.G1 * axis_fp
    bandwidth_hz = R_value * 1.0e6 / pulse_width_us
    if abs(grad_total) < 1.0e-9 or abs(P.GradCal_HzPerCm) < 1.0e-9:
        comms.log("WARNING: cannot compute excited width for the {0} pulse "
                  "-- G1*FP{1} or GradCal_HzPerCm is ~0 (grad_total={2:.4f}, "
                  "GradCal_HzPerCm={3}). Excited bandwidth={4:.1f} Hz, but "
                  "width cannot be converted to a physical length."
                  .format(shape_label, P.Axis.upper(), grad_total,
                          P.GradCal_HzPerCm, bandwidth_hz))
        return bandwidth_hz, None, grad_total
    width_cm = bandwidth_hz / (P.GradCal_HzPerCm * abs(grad_total))
    comms.log("*** Excited {0} slice: bandwidth = {1:.1f} Hz (R={2:.2f}/"
              "{3:.1f}us), width = {4:.3f} cm ({5:.2f} mm) on the {6}-axis "
              "at G1*FP={7:.4f} -- GradCal_HzPerCm={8} (VERIFY this "
              "calibration on a real phantom before trusting the width "
              "number). ***"
              .format(shape_label, bandwidth_hz, R_value, pulse_width_us,
                      width_cm, width_cm*10.0, P.Axis.upper(), grad_total,
                      P.GradCal_HzPerCm))
    return bandwidth_hz, width_cm, grad_total


@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_Slice-Selective-Echo_H", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 59.7, ParameterTypes.Double, "H/F Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "H/F Freq Offset [Hz] -- sequence-wide ref.; use SlicePPM to position the slice")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "H/F TX Freq Offset [ppm]")

    # *** THE clear, obvious place to set the slice/pulse centre frequency ***
    SlicePPM = Parameter("SlicePPM", 0.0, ParameterTypes.Double, "*** Slice/Pulse Centre Freq [ppm] rel. SF+O1 -- SET to position the excited slice (both shaped pulses) ***")
    PulseOffsetTrim = Parameter("PulseOffsetTrim", 0.0, ParameterTypes.Double, "Fine-trim [Hz] on top of SlicePPM x SF -- sub-ppm adjustment only, normally 0")

    ReceiverPoints = Parameter("NP", 1024, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 34, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "100000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans
    NumScans = Parameter("NS", 8, ParameterTypes.Int32,"Scans")
    DS = Parameter("DS", 0, ParameterTypes.Int32,"Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 100.0, ParameterTypes.Double, "H/F Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 500000, ParameterTypes.Int32, "Relaxation Delay [s]",
                              RD, min=100000, max=2000000000)

    # Shaped (spatially slice-selective) pulses -- synthesised ON THE FLY.
    # Defaults = literature pulses: EXCITATION = E-BURP-1 (Geen &amp; Freeman
    # 1991, Table 2, nmax=8), REFOCUSING = RE-BURP (Table 8, Np=256).
    ExcitationShape = Parameter("ExShape", "EBURP1", ParameterTypes.String, "Excitation (90) Shape [EBURP1(default)/GAUSSIAN/SINC/BURP] -- see ExBurpCoeffsA/B")
    RefocusShape = Parameter("RefocusShape", "REBURP", ParameterTypes.String, "Refocusing (180) Shape [REBURP(default)/GAUSSIAN/SINC/BURP] -- see RefocusBurpCoeffsA/B")
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00", ParameterTypes.String, "E-BURP-1 cosine coeffs A0..A8 (Geen &amp; Freeman'91 Tbl.2)")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01", ParameterTypes.String, "E-BURP-1 sine coeffs B0..B8 (Geen &amp; Freeman'91 Tbl.2)")
    RefocusBurpCoeffsA = Parameter("RefocusBurpCoeffsA", "0.49,-1.02,1.11,-1.57,0.83,-0.42,0.26,-0.16,0.10,-0.07,0.04,-0.03,0.01,-0.02,0.00,-0.01", ParameterTypes.String, "RE-BURP cosine coeffs A0..A15 (Geen &amp; Freeman'91 Tbl.8)")
    RefocusBurpCoeffsB = Parameter("RefocusBurpCoeffsB", "", ParameterTypes.String, "RE-BURP sine coeffs -- none published, leave empty")

    P90sh = Parameter("P90sh", 5000.0, ParameterTypes.Double, "Shaped 90&#176; (excitation) Width [&#956;s] -- also sets shape res. + excitation slice width")
    # P180sh's displayed default already matches P90sh's default via the
    # R_REBURP_180/R_EBURP1_90 ratio (5.9/4.9), so the two slices are
    # consistent even before AutoMatchSlice recomputes it at run time.
    P180sh = Parameter("P180sh", 6020.0, ParameterTypes.Double, "Shaped 180&#176; (refocusing) Width [&#956;s] -- IGNORED/recomputed from P90sh if AutoMatchSlice=1")
    AutoMatchSlice = Parameter("AutoMatchSlice", 1, ParameterTypes.Int32, "*** Auto-derive P180sh from P90sh to match slices [1=On(default,RECOMMENDED), 0=manual P180sh] ***")

    TXAmplitude90 = Parameter("RFAsh0", 0.30, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0&#8230;1.0 of LP max] -- CALIBRATE via nutation, see LWG_Selective-Echo_H.py", RFA, min=0.0, max=1.0)
    TXAmplitude180 = Parameter("RFAsh1", 0.30, ParameterTypes.Double,
                               "Shaped 180&#176; TX Power [0&#8230;1.0 of LP max] -- CALIBRATE via nutation, see LWG_Selective-Echo_H.py", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", 0.10, ParameterTypes.Double, "LP Port Max as Fraction of HP [0&#8230;1.0] -- your measured value")

    # OPTIONAL power calculator -- OFF by default. See LWG_Selective-Echo_H.py
    # / LWG_Selective-Pulse-Power-Calibration.md for the full explanation.
    PowerCalcMethod = Parameter("PowerCalcMethod", "manual", ParameterTypes.String, "RFAsh0/RFAsh1 source: manual(default)/shape(recommended)/db")
    RefAmplitude_HP = Parameter("RefAmplitude_HP", 0.40, ParameterTypes.Double, "Ref. hard-pulse rel. amplitude [0&#8230;1] on HP port, known flip angle -- for shape/db methods")
    HardPulseWidth = Parameter("P1Hard", 9.58, ParameterTypes.Double, "Ref. hard 90&#176; width [&#956;s] at RefAmplitude_HP on HP port -- for 'shape' method")
    ExcitationRotation = Parameter("ExRotation", 90.0, ParameterTypes.Double, "Target rotation of excitation shape [&#176;] (EBURP1=90) -- for 'shape' method")
    RefocusRotation = Parameter("RefocusRotation", 180.0, ParameterTypes.Double, "Target rotation of refocusing shape [&#176;] (REBURP=180) -- for 'shape' method")
    PowerAdjust_dB = Parameter("PowerAdjust_dB", 0.0, ParameterTypes.Double, "Manual fine-tune [dB] on top of 'shape' calc -- for PowerCalcMethod='shape'")
    Excitation_dB = Parameter("Excitation_dB", 0.0, ParameterTypes.Double, "Excitation power as dB attenuation vs RefAmplitude_HP -- for 'db' method")
    Refocus_dB = Parameter("Refocus_dB", 0.0, ParameterTypes.Double, "Refocusing power as dB attenuation vs RefAmplitude_HP -- for 'db' method")
    DbConvention = Parameter("DbConvention", "amplitude", ParameterTypes.String, "dB convention: amplitude(B1,usual)/power -- for 'db' method")

    # Slice-select gradient -- SAME gradient (G1, same axis) is used for
    # BOTH pulses, by design (see 'MATCHING THE EXCITATION AND REFOCUSING
    # SLICES' above) -- there is deliberately no separate G-for-180 knob.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Slice-Select Gradient Strength [-1.0&#8230;1.0] -- used for BOTH pulses")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Slice-Select Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z", ParameterTypes.String, "Gradient Axes")

    RampTime = Parameter("D70", 200.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 200.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s]")
    PreGrad = Parameter("D75", 200.0, ParameterTypes.Double, "Pre-Gradient Time (before each ramp-up) [&#956;s]")
    RephaseFraction = Parameter("RephaseFrac", 0.5, ParameterTypes.Double, "Excitation-only slice rephase lobe area, fraction of excitation gradient area -- 0.5 default, fine-tune empirically")

    # Echo timing
    Tau = Parameter("TAU", 20000, ParameterTypes.Int32, "Echo &#964; Delay [&#956;s], pulse-centre to pulse-centre -- must exceed half pulse widths + crusher/gradient margins")

    # Crusher gradients around the refocusing pulse (same axis, same
    # polarity, equal area both sides -- see design notes)
    CrusherOn = Parameter("CrusherOn", 1, ParameterTypes.Int32, "Crusher gradients around the 180 [1=On(default), 0=Off]")
    CrusherStrength = Parameter("CrusherG", 0.30, ParameterTypes.Double, "Crusher gradient strength [-1.0&#8230;1.0] -- keep modest relative to G1 (shares slice-select axis)")
    CrusherTime = Parameter("CrusherT", 1000.0, ParameterTypes.Double, "Crusher gradient hold time [&#956;s] (each side)")

    # *** Gradient calibration -- UNCALIBRATED PLACEHOLDER, see design notes ***
    GradCal_HzPerCm = Parameter("GradCal", 1000.0, ParameterTypes.Double, "*** PLACEHOLDER, UNCALIBRATED *** Gradient calibration [Hz/cm] at G1*FP=1.0 -- determine empirically before trusting ExcitedWidth")

    # Mains-lock trigger -- OFF by default.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rails (conservative placeholders, not vendor specs)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0&#8230;1.0]")

    # Phases -- combined 2-step (90) x EXORCYCLE 4-step (180) = 8-step cycle.
    # ph1 cancels simple 90/DC artefacts; ph2 (EXORCYCLE: 0,90,180,270 on the
    # refocusing pulse) cancels signal from magnetisation NOT correctly
    # inverted by an imperfect real 180 (that pathway picks up 1x ph2 per
    # step, the true echo picks up 2x ph2 per step -- the receiver is
    # cycled to match the TRUE echo's 2x-per-step pattern, so only the
    # artefact pathway is cancelled by the sum over the 4 EXORCYCLE steps).
    PH1 = Parameter("PH1", "0,0,0,0,180,180,180,180", ParameterTypes.String, "Shaped 90&#176; Pulse Phase")
    PH2 = Parameter("PH2", "0,90,180,270,0,90,180,270", ParameterTypes.String, "Shaped 180&#176; Pulse Phase (EXORCYCLE)")
    PHRX = Parameter("PHRX", "0,180,0,180,180,0,180,0", ParameterTypes.String, "Acquisition Phase")

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

    # SlicePPM (ppm relative to SF+O1) -> Hz, plus fine-trim -- THE pulse
    # centre frequency, applied to BOTH shaped pulses.
    PulseOffsetHz = (P.SlicePPM * P.FrequencyBase) + P.PulseOffsetTrim

    ReceiverFilter = Filter(P.Filter)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell

    ReceiverTime = DW*(points+1)
    times = np.arange(0, points*DW, DW) / 1.0e6

    # ---- Slice matching: derive P180sh from P90sh unless overridden -------
    RValueRatio = R_REBURP_180 / R_EBURP1_90
    if int(P.AutoMatchSlice) != 0:
        P180sh_eff = P.P90sh * RValueRatio
        comms.log("AutoMatchSlice=1: P180sh overridden to {0:.2f} us (from "
                  "P90sh={1}us x R_REBURP_180/R_EBURP1_90={2:.4f}) so the "
                  "excitation and refocusing pulses excite/refocus the SAME "
                  "slice width at the same G1. The panel's P180sh value is "
                  "ignored while AutoMatchSlice=1."
                  .format(P180sh_eff, P.P90sh, RValueRatio))
    else:
        P180sh_eff = P.P180sh
        comms.log("AutoMatchSlice=0: using manually-entered P180sh={0}us -- "
                  "slice widths will be checked and a WARNING issued below "
                  "if they don't match.".format(P180sh_eff))

    ExcitationPulseWidth = P.P90sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD
    RefocusPulseWidth = P180sh_eff + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD

    ExBurpA = parse_burp_coeffs(P.ExBurpCoeffsA)
    ExBurpB = parse_burp_coeffs(P.ExBurpCoeffsB)
    RefocusBurpA = parse_burp_coeffs(P.RefocusBurpCoeffsA)
    RefocusBurpB = parse_burp_coeffs(P.RefocusBurpCoeffsB)

    # ---- Optional power calculator override -------------------------------
    RFAsh0_effective = P.TXAmplitude90
    RFAsh1_effective = P.TXAmplitude180
    PowerMethod = str(P.PowerCalcMethod).strip().lower()

    if PowerMethod == 'shape':
        ExAmpProfile, _ = generate_shape(P.ExcitationShape, max(4, int(round(P.P90sh))), ExBurpA, ExBurpB)
        RefocusAmpProfile, _ = generate_shape(P.RefocusShape, max(4, int(round(P180sh_eff))), RefocusBurpA, RefocusBurpB)
        ExIntegFactor = shape_integration_factor(ExAmpProfile)
        RefocusIntegFactor = shape_integration_factor(RefocusAmpProfile)
        RFAsh0_effective = estimate_shape_relative_scale(
            P.HardPulseWidth, P.RefAmplitude_HP, P.P90sh, P.ExcitationRotation,
            ExIntegFactor, P.LPMaxFraction, P.PowerAdjust_dB)
        RFAsh1_effective = estimate_shape_relative_scale(
            P.HardPulseWidth, P.RefAmplitude_HP, P180sh_eff, P.RefocusRotation,
            RefocusIntegFactor, P.LPMaxFraction, P.PowerAdjust_dB)
        comms.log("PowerCalcMethod='shape': excitation integ.factor={0:.4f} "
                  "-> RFAsh0={1:.4f}; refocus integ.factor={2:.4f} -> "
                  "RFAsh1={3:.4f}.".format(ExIntegFactor, RFAsh0_effective,
                          RefocusIntegFactor, RFAsh1_effective))
        if RFAsh0_effective > 1.0 or RFAsh1_effective > 1.0:
            comms.log("WARNING: calculated relative scale exceeds 1.0 -- "
                      "RFAsh0_effective={0:.4f}, RFAsh1_effective={1:.4f}."
                      .format(RFAsh0_effective, RFAsh1_effective))
    elif PowerMethod == 'db':
        RFAsh0_effective = db_to_lp_relative_scale(P.Excitation_dB, P.RefAmplitude_HP, P.LPMaxFraction, P.DbConvention)
        RFAsh1_effective = db_to_lp_relative_scale(P.Refocus_dB, P.RefAmplitude_HP, P.LPMaxFraction, P.DbConvention)
        comms.log("PowerCalcMethod='db': RFAsh0={0:.4f}, RFAsh1={1:.4f}."
                  .format(RFAsh0_effective, RFAsh1_effective))
        if RFAsh0_effective > 1.0 or RFAsh1_effective > 1.0:
            comms.log("WARNING: dB-calculated relative scale exceeds 1.0 -- "
                      "RFAsh0_effective={0:.4f}, RFAsh1_effective={1:.4f}."
                      .format(RFAsh0_effective, RFAsh1_effective))
    elif PowerMethod != 'manual':
        comms.log("WARNING: unrecognised PowerCalcMethod='{0}' -- falling "
                  "back to 'manual'.".format(P.PowerCalcMethod))

    # ---- Excited slice width: CALCULATED AND RECORDED every run, for BOTH
    # pulses ------------------------------------------------------------
    ExBW_Hz, ExWidth_cm, GradTotal = compute_excited_width(
        P, comms, P.P90sh, R_EBURP1_90, "EXCITATION (90)")
    RefocusBW_Hz, RefocusWidth_cm, _ = compute_excited_width(
        P, comms, P180sh_eff, R_REBURP_180, "REFOCUSING (180)")
    if ExWidth_cm is not None and RefocusWidth_cm is not None:
        mismatch = abs(ExWidth_cm - RefocusWidth_cm) / ExWidth_cm
        if mismatch > 0.02:
            comms.log("WARNING: excitation slice ({0:.3f} cm) and "
                      "refocusing slice ({1:.3f} cm) differ by {2:.1%} -- "
                      "the echo will only reflect their OVERLAP, narrower "
                      "than either pulse's own width. Set AutoMatchSlice=1, "
                      "or manually set P180sh = P90sh x {3:.4f}."
                      .format(ExWidth_cm, RefocusWidth_cm, mismatch, RValueRatio))
        else:
            comms.log("Excitation and refocusing slice widths match to "
                      "within {0:.2%} ({1:.3f} cm vs {2:.3f} cm) -- good."
                      .format(mismatch, ExWidth_cm, RefocusWidth_cm))

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
                "jc_slicePPM={0}".format(P.SlicePPM),
                "jc_excitedBW_ex_Hz={0:.2f}".format(ExBW_Hz),
                "jc_excitedWidth_ex_cm={0}".format("{0:.4f}".format(ExWidth_cm) if ExWidth_cm is not None else "uncalibrated"),
                "jc_excitedBW_ref_Hz={0:.2f}".format(RefocusBW_Hz),
                "jc_excitedWidth_ref_cm={0}".format("{0:.4f}".format(RefocusWidth_cm) if RefocusWidth_cm is not None else "uncalibrated"),
                "jc_P180sh_effective_us={0:.2f}".format(P180sh_eff),
                    ]
        return jc_data

    comms.log("Sequence Time: {0}".format(time_calculation(Parameters)))
    seqAcqu = CallBack1D(Parameters, comms, ReceiverFilter, jcamp_meta())

    def RecvCallback(acqData, scan):
        seqAcqu.process_data(scan, acqData)

    with sequential:

        mains_lock_trigger(P, comms)

        Transmit1SelectPort(1,0)
        Transmit1LPEnable(1,1)
        Delay(10000) # relay settle

        Receiver1Preamp(128, P.ReceiverAttenuation)
        Receiver1Filter(200, ReceiverFilter)

        GradientMatrix(200,Matrix.T)

        Phases.Reset()

    comms.log("LP port engaged for the shaped pulses -- LPMaxFraction="
              "{0:.2%} of HP full scale. Slice/pulse offset: SlicePPM={1} "
              "ppm x SF={2}MHz + PulseOffsetTrim={3}Hz = {4:.2f} Hz "
              "(applied on top of SF+O1={5:.6f} MHz)."
              .format(P.LPMaxFraction, P.SlicePPM, P.FrequencyBase,
                      P.PulseOffsetTrim, PulseOffsetHz, Frequency))

    # ---- Exact timing accounting (all in us) -------------------------------
    # NOTE ON DESIGN: each slice-select gradient is synchronised with its
    # shaped pulse by simple SEQUENTIAL composition (ramp the gradient up,
    # call shaped_pulse() completely unmodified, ramp it down/rephase) --
    # NOT by running them in separate 'with parallel:' branches. A gradient
    # command's amplitude persists until the next one changes it, so this
    # leaves the gradient flat for exactly the shaped pulse's duration
    # without nesting a 'with parallel:' inside shaped_pulse()'s own
    # internal 3-way parallel block (Transmit1/amplitude/phase) -- a
    # combination not demonstrated anywhere in the Pulse Sequence
    # Programming manual's examples (which only ever show flat, single-
    # level 'with parallel:' blocks). This also means there is only ONE
    # timeline per pulse (no branches to keep in sync), so the 'TAU' maths
    # below is a straightforward sum, not a branch-padding exercise.

    # Excitation event: ramp-up, shaped pulse dwell, ramp-down, rephase lobe
    # (ramp/hold/ramp), settle.
    ExGradHoldRephase = P.RephaseFraction * ExcitationPulseWidth
    ExGradEventTime = (P.PreGrad + P.RampTime + ExcitationPulseWidth + P.RampTime
                        + P.RampTime + ExGradHoldRephase + P.RampTime
                        + P.GradSettle)
    # Time from the excitation pulse's CENTRE to the end of its gradient
    # event (ramp-down + rephase + settle, plus the second half of the
    # shaped pulse itself) -- how much still needs to be 'paid for' out of
    # the first TAU.
    ExCentreToEventEnd = (ExcitationPulseWidth/2.0 + P.RampTime + P.RampTime
                           + ExGradHoldRephase + P.RampTime + P.GradSettle)

    # Refocusing event: ramp-up, shaped pulse dwell, ramp-down, settle --
    # NO rephase lobe (self-refocusing, see design notes).
    RefocusGradEventTime = (P.PreGrad + P.RampTime + RefocusPulseWidth + P.RampTime
                         + P.GradSettle)
    RefocusLeadToCentre = P.PreGrad + P.RampTime + RefocusPulseWidth/2.0
    RefocusCentreToEventEnd = RefocusPulseWidth/2.0 + P.RampTime + P.GradSettle

    CrusherBlockTime = (P.PreGrad + P.RampTime + P.CrusherTime + P.RampTime + P.GradSettle) if int(P.CrusherOn) != 0 else 0.0

    rf_on_time = ExcitationPulseWidth + RefocusPulseWidth
    grad_on_time = ExGradEventTime + RefocusGradEventTime + 2*CrusherBlockTime
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver1FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            slew_rate_fn, gradient_fn = get_gradient_functions(P.Axis)

            # ---- Excitation: E-BURP-1 + slice-select gradient + rephase --
            # (sequential composition -- see GradEventTime design note above)
            Delay(P.PreGrad)
            slew_rate_fn(1.0, abs(P.G1/P.RampTime))
            gradient_fn(P.RampTime, P.G1)
            shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, RFAsh0_effective, P.TXEnableTime, Frequency, PulseOffsetHz, ExBurpA, ExBurpB)
            gradient_fn(P.RampTime, 0)
            slew_rate_fn(1.0, abs(P.G1/P.RampTime))
            gradient_fn(P.RampTime, -P.G1)
            safe_delay(ExGradHoldRephase, "excitation rephase-lobe hold (check RephaseFraction)", comms)
            gradient_fn(P.RampTime, 0)
            gradient_fn(1.0, 0)
            Delay(P.GradSettle)

            # ---- Crusher (pre-180) ----------------------------------------
            if int(P.CrusherOn) != 0:
                apply_gradient(P.Axis, P.CrusherTime, P.RampTime, P.CrusherStrength, P.GradSettle, P.PreGrad)

            # ---- First TAU: excitation-centre to refocusing-centre --------
            safe_delay(P.Tau - ExCentreToEventEnd - CrusherBlockTime - RefocusLeadToCentre,
                       "first TAU wait", comms)

            # ---- Refocusing: RE-BURP + slice-select gradient (matched) ---
            # NO rephase lobe -- RE-BURP is self-refocusing by symmetry.
            Delay(P.PreGrad)
            slew_rate_fn(1.0, abs(P.G1/P.RampTime))
            gradient_fn(P.RampTime, P.G1)
            shaped_pulse(P180sh_eff, ph["PH2"], P.RefocusShape, RFAsh1_effective, P.TXEnableTime, Frequency, PulseOffsetHz, RefocusBurpA, RefocusBurpB)
            gradient_fn(P.RampTime, 0)
            gradient_fn(1.0, 0)
            Delay(P.GradSettle)

            # ---- Crusher (post-180) ----------------------------------------
            if int(P.CrusherOn) != 0:
                apply_gradient(P.Axis, P.CrusherTime, P.RampTime, P.CrusherStrength, P.GradSettle, P.PreGrad)

            # ---- Second TAU: refocusing-centre to acquisition -------------
            # ReceiverFilter.dead_time is reserved out of this wait and paid
            # back explicitly (its own Delay, right before Receiver1 below).
            safe_delay(P.Tau - RefocusCentreToEventEnd - CrusherBlockTime - ReceiverFilter.dead_time,
                       "second TAU wait", comms)

            # ---- ACQU -------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 06/08/26 - Initial version. Built as the genuinely SPATIALLY
#    slice-selective counterpart to LWG_Selective-Echo_H.py (which is
#    frequency/chemical-shift-selective only -- confirmed by you). Adds a
#    slice-select gradient synchronised with each shaped pulse by simple
#    SEQUENTIAL composition (ramp the gradient up, call shaped_pulse()
#    completely unmodified, ramp it down) rather than a 'with parallel:'
#    branch, so the gradient sits flat for exactly the shaped pulse's
#    duration without nesting a 'with parallel:' inside shaped_pulse()'s own
#    internal 3-way parallel block -- a combination not demonstrated
#    anywhere in the Pulse Sequence Programming manual's examples. Runs on
#    the SAME gradient (G1, same axis) for both the E-BURP-1 excitation and
#    RE-BURP refocusing pulses. Adds AutoMatchSlice (default on), which derives
#    P180sh from P90sh via the ratio of the two shapes' numerically-derived
#    bandwidth-time products (R_EBURP1_90=4.9, R_REBURP_180=5.9 -- see their
#    docstrings for the Bloch-simulation method), so the excitation and
#    refocusing slices are the SAME width by construction -- this is the
#    "equally sized slices" requirement, addressed directly. Adds an
#    excitation-only rephase gradient lobe (RephaseFraction, default 0.5;
#    RE-BURP needs none, being self-refocusing by symmetry -- see design
#    notes). Adds crusher gradients around the refocusing pulse (CrusherOn,
#    default on) plus an 8-step combined CYCLOPS+EXORCYCLE phase cycle, for
#    artefact suppression against imperfect real-world 180 pulses. Adds
#    SlicePPM (prominent, ppm) as the pulse-centre-frequency control,
#    matching LWG_Slice-Selective-PulseAcquire_H.py. Adds
#    compute_excited_width(), calculating and comms.log()ing (and recording
#    into JCAMP metadata) both pulses' excited bandwidth/width every run,
#    with an explicit WARNING if AutoMatchSlice=0 and the two widths
#    mismatch by more than 2%. Everything else (CallBack1D/Parameters
#    scaffolding, generate_shape()/shaped_pulse(), the three-way
#    PowerCalcMethod calculator, mains_lock_trigger(), estimate_duty_cycles(),
#    safe_delay() guard, CRLF line endings) is carried over from
#    LWG_Selective-Echo_H.py / LWG_1D-Image-Echo_H.py.
#
# -----------------------------------------------------------------------------
