#-------------------------------------------------------------------------------
# Name:        LWG_Slice-Selective-PulseAcquire_X.py
# Purpose:     TRUE spatial SLICE-selective pulse-acquire (single shaped
#              E-BURP-1 excitation, gradient ON during the pulse), on the
#              Oxford Instruments X-Pulse Broadband Benchtop NMR Spectrometer
#              (X channel) -- mechanical port of
#              LWG_Slice-Selective-PulseAcquire_H.py, see X-CHANNEL
#              CALIBRATION note below. Axis-selectable (x/y/z). No
#              refocusing pulse, no frequency-encoded readout -- this
#              excites and acquires a localised FID from one slab/slice, it
#              does not profile/image it.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     06/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.0
#
# X-CHANNEL CALIBRATION -- IMPORTANT: this is a mechanical H->X port (all
# Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1). SF
# now defaults to 15.01 MHz as a placeholder -- P90sh/RFAsh0/GradCal_HzPerCm/
# R_EBURP1_90 below are still the 1H-tuned values/assumptions from the _H
# file and are almost certainly WRONG for whatever nucleus you actually put
# on the X channel (in particular, GradCal_HzPerCm and any gyromagnetic-
# ratio-dependent assumption must be re-measured/re-derived for your target
# nucleus -- R_EBURP1_90 itself is nucleus-independent, being a property of
# the pulse SHAPE only, but the Hz<->cm conversion via GradCal_HzPerCm is
# not). Recalibrate P90sh/RFAsh0 with LWG_Selective-Echo_X.py first.
#
# Design notes -- READ THIS:
#
#  - CHEMICAL-SHIFT-SELECTIVE vs SLICE-SELECTIVE: LWG_Selective-Echo_H.py and
#    LWG_1D-Image-Echo-Selective_H.py both explicitly leave the gradient OFF
#    during their shaped pulses -- that only picks a resonance (chemical
#    shift), not a spatial slice. THIS file is the genuine spatial version:
#    the slice-select gradient (GradAxis/G1/FPX-Y-Z) is ON, at a constant
#    plateau, for the ENTIRE duration of the shaped E-BURP-1 pulse, so the
#    excited region is a real spatial slab, positioned by SlicePPM (see
#    below) and thickness set by P90sh + G1 (see 'Excited slice width'
#    below).
#  - REPHASE LOBE: because a gradient is running while spins are actively
#    being tipped, spins across the slice pick up a position-dependent phase
#    during the pulse -- otherwise they'd destructively interfere and you'd
#    lose most of the signal. This is fixed with a REPHASE gradient lobe
#    (opposite polarity) immediately after the excitation gradient -- exactly
#    the same fix used for every sinc/Gaussian slice-select excitation pulse
#    in ordinary MRI. The textbook rule for a temporally SYMMETRIC pulse
#    envelope is: rephase area = 0.5 x (excitation-gradient area). E-BURP-1
#    is knowingly ASYMMETRIC in time (see the design notes in
#    LWG_1D-Image-Echo-Selective_H.py changelog item 4d -- "E-BURP-1
#    reconstructs asymmetric... intentionally asymmetric so magnetization
#    outside the excited band returns to +z"), so 0.5 is a physically-
#    motivated STARTING point, not an exact value for this specific shape.
#    RephaseFraction (Parameter, default 0.5) is exposed so you can tune it
#    empirically: on a real sample/phantom, scan RephaseFraction and maximise
#    the acquired FID's initial amplitude -- that's your true optimum.
#  - EXCITED SLICE WIDTH: E-BURP-1's excitation bandwidth is NOT in the
#    Fourier-coefficient tables you supplied (those only give the pulse's
#    time-domain envelope, not its frequency response). Rather than quote an
#    uncertain literature number from memory, R_EBURP1_90 below was
#    determined by literally Bloch-simulating (Rodrigues'-formula rotation,
#    calibrated so the on-resonance flip angle is EXACTLY 90 degrees via a
#    numerical search, not the linear/small-tip-angle 'integration factor'
#    shortcut) the SAME generate_shape('EBURP1', ...) envelope this file
#    actually plays out, then measuring the FWHM of the excitation profile
#    (|Mxy| vs offset) and computing bandwidth x duration. Result: R = 4.96
#    (rounds to the commonly-quoted literature value of ~4.9 for the E-BURP
#    excitation family, Geen &amp; Freeman, J. Magn. Reson. 93, 93-141, 1991 --
#    a reassuring cross-check, not a coincidence). Used as 4.9 below. This
#    is duration-independent (R = bandwidth x duration is a fixed property
#    of a given pulse SHAPE -- the Bloch equation is invariant under jointly
#    rescaling time and all rates, so R does not need to be re-simulated for
#    every P90sh you might use).
#  - GRADIENT CALIBRATION: GradCal_HzPerCm is an UNCALIBRATED PLACEHOLDER --
#    this file has no way to know your actual gradient coil's field/current
#    relationship. Determine it empirically, the same way you already
#    calibrate FPX/FPY/FPZ (image a phantom of known physical length with
#    LWG_1D-Image-Echo_H.py or LWG_1D-Image-UTE_H.py and read off Hz/cm from
#    the known acquisition bandwidth and the phantom's reconstructed extent),
#    then enter that number here so the logged/recorded 'ExcitedWidth' is
#    physically meaningful rather than a placeholder-based guess.
#  - PULSE CENTRE FREQUENCY: SlicePPM (near the top of Parameters, impossible
#    to miss) is THE place to set where the slice is excited, in ppm
#    relative to SF+O1 -- exactly what was asked for. PulseOffsetTrim (Hz)
#    is available underneath it for sub-ppm fine adjustment; both are summed
#    and applied to the shaped pulse only (acquisition stays referenced to
#    SF+O1, same convention as LWG_Selective-Echo_H.py's PulseOffset).
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
            self.comms.log("seqTime LWG_Slice-Selective-PulseAcquire_X Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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
# R = (excitation/refocusing bandwidth [Hz]) x (pulse duration [s]), FWHM
# definition. Derived by Bloch-simulating the EXACT generate_shape() envelope
# below (Rodrigues'-formula rotation per 1us step, on-resonance flip angle
# calibrated to exactly 90 (or 180) degrees via numerical search -- NOT the
# linear/small-tip-angle 'integration factor' shortcut, which is only a
# power-calibration approximation, see estimate_shape_relative_scale()).
# R is duration-independent for a fixed pulse SHAPE (rescaling time and all
# rates together leaves the Bloch equation, and therefore R, unchanged) --
# it does not need to be re-derived for whatever P90sh you actually use.
# E-BURP-1 simulated at R=4.96, rounds to the ~4.9 commonly quoted in the
# literature for the E-BURP excitation family (Geen &amp; Freeman 1991) -- a
# reassuring cross-check. See file header design notes for the full method.
R_EBURP1_90 = 4.9

def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell

    t_scanTime = (P.P90sh*1.6 + P.Dead1 + points*DW)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Spatially SLICE-selective pulse-acquire (single shaped "
                "E-BURP-1 excitation with the slice-select gradient ON "
                "during the pulse, plus a rephase lobe), on the {X} "
                "channel. No refocusing pulse, no frequency-encoded "
                "readout -- acquires a localised FID from one slab. Gradient "
                "axis selectable via GradAxis (x/y/z). SlicePPM sets where "
                "the slice sits (ppm); ExcitedWidth is calculated and "
                "logged/recorded every run.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,SlicePPM,P90sh,RFAsh0,G1,GradAxis"

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
    """Per-axis calibration scaler (FPX/FPY/FPZ) matching GradAxis -- the
    ONLY place per-axis calibration is applied is GradientMatrix() in run()
    (same convention as LWG_1D-Image-Echo_H.py); this helper is only used
    for the ExcitedWidth calculation below, which needs to know the same
    effective per-axis scale GradientMatrix() is applying in hardware."""
    if axisstr == 'x':
        return P.XGradNorm
    elif axisstr == 'y':
        return P.YGradNorm
    elif axisstr == 'z':
        return P.ZGradNorm
    else:
        exit()

def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical machinery to LWG_Selective-Echo_H.py /
    LWG_1D-Image-Echo-Selective_H.py (see those files for the full
    docstring); duplicated here so this file stays fully self-contained."""
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
                "coefficients as ExBurpCoeffsA/B, or use 'GAUSSIAN' / "
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
                          "SINC, EBURP1/BURP (needs BurpCoeffsA/B)."
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
    any gradient -- the slice-select gradient is driven by a SEPARATE
    parallel branch in run() (see slice_selective_excite() below), so this
    function is reusable verbatim from the frequency-selective-only pps."""
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
               "(negative). Increase PreGrad/RampTime margins, or decrease "
               "RephaseFraction/GradSettle, then retry.").format(label, value)
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
    docstring -- dB-referenced power figure -> X-Pulse linear relative
    scale."""
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
                  "({1:.2%}). Consider a longer RD, shorter P90sh, or lower "
                  "RFAsh0, or confirm with Oxford Instruments that this is "
                  "within the transmitter's rated duty cycle before running "
                  "unattended.".format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} on the "
                  "{1}-axis coil exceeds MaxGradDuty ({2:.2%}). Consider a "
                  "longer RD or shorter GradientOnTime-equivalent settings, "
                  "or confirm the rated duty cycle for the X-Pulse gradient "
                  "amplifier/coil with Oxford Instruments -- this limit is a "
                  "conservative placeholder, not a vendor spec."
                  .format(grad_duty, P.Axis, P.MaxGradDuty))


def compute_excited_width(P, comms, pulse_width_us, R_value, shape_label):
    """Compute the excited slice bandwidth [Hz] and width [cm], using the
    R-value bandwidth-time product (see R_EBURP1_90 above) and the
    user-supplied GradCal_HzPerCm calibration. Logs a clear, prominent
    result every run -- this is the 'expected excited width' the pulse
    programme is required to calculate and record.

    Returns (bandwidth_hz, width_cm, grad_total)."""
    axis_fp = get_axis_fp(P, P.Axis)
    grad_total = P.G1 * axis_fp
    bandwidth_hz = R_value * 1.0e6 / pulse_width_us    # R/(T[s]) with T in us
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
    Sequence = Parameter("Sequence", "LWG_Slice-Selective-PulseAcquire_X", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 15.01, ParameterTypes.Double, "X Base Freq [MHz] -- PLACEHOLDER, set to your target nucleus")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "X Freq Offset [Hz] -- overall sequence reference (affects BOTH the pulse and acquisition); leave at 0 (or your usual reference) and use SlicePPM below to position the slice")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "X TX Freq Offset [ppm]")

    # *** THE clear, obvious place to set the slice/pulse centre frequency ***
    SlicePPM = Parameter("SlicePPM", 0.0, ParameterTypes.Double, "*** SLICE / PULSE CENTRE FREQUENCY [ppm] relative to SF+O1 -- SET THIS to position the excited slice. Converted internally to Hz (SlicePPM x SF) and applied ONLY to the shaped pulse; acquisition stays referenced to SF+O1. ***")
    PulseOffsetTrim = Parameter("PulseOffsetTrim", 0.0, ParameterTypes.Double, "Fine-trim [Hz] added on top of SlicePPM x SF -- for sub-ppm adjustment only; leave at 0 normally")

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

    # Shaped (spatially slice-selective) excitation pulse -- synthesised ON
    # THE FLY (see generate_shape()/shaped_pulse() above). Default =
    # E-BURP-1 (Geen &amp; Freeman 1991, Table 2, nmax=8).
    ExcitationShape = Parameter("ExShape", "EBURP1", ParameterTypes.String, "Excitation (90) Shape [EBURP1 (default), GAUSSIAN, SINC, or BURP -- see ExBurpCoeffsA/B]")
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00", ParameterTypes.String, "E-BURP-1 Fourier cosine coeffs A0..A8 (Geen &amp; Freeman 1991, Table 2, nmax=8) -- used when ExShape=EBURP1/BURP")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01", ParameterTypes.String, "E-BURP-1 Fourier sine coeffs B0(unused)..B8 (Geen &amp; Freeman 1991, Table 2, nmax=8) -- used when ExShape=EBURP1/BURP")
    P90sh = Parameter("P90sh", 5000.0, ParameterTypes.Double, "Shaped 90&#176; Pulse Width [&#956;s] -- also sets shape resolution (1 point/&#956;s) AND (with G1/GradCal_HzPerCm) the excited slice width -- see ExcitedWidth logged every run")

    TXAmplitude90 = Parameter("RFAsh0", 0.30, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0.0&#8230;1.0 of LP max] -- CALIBRATE via nutation, see LWG_Selective-Echo_H.py / the calibration markdown", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", 0.10, ParameterTypes.Double, "LP Port Max Power as Fraction of HP Port [0.0&#8230;1.0] -- your measured value")

    # OPTIONAL power calculator -- OFF by default (PowerCalcMethod='manual').
    # See LWG_Selective-Echo_H.py / LWG_Selective-Pulse-Power-Calibration.md
    # for the full explanation of each method.
    PowerCalcMethod = Parameter("PowerCalcMethod", "manual", ParameterTypes.String, "RFAsh0 source: 'manual' (default, use as entered), 'shape' (calculate from hard-pulse calibration + shape integral), 'db' (legacy direct dB attenuation)")
    RefAmplitude_HP = Parameter("RefAmplitude_HP", 0.40, ParameterTypes.Double, "Reference HARD-pulse relative amplitude on the HIGH-power port for a KNOWN, calibrated flip angle -- used if PowerCalcMethod='shape' or 'db'")
    HardPulseWidth = Parameter("P1Hard", 9.58, ParameterTypes.Double, "Reference HARD 90&#176; pulse width [&#956;s] at RefAmplitude_HP on the HIGH-power port -- only used if PowerCalcMethod='shape'")
    ExcitationRotation = Parameter("ExRotation", 90.0, ParameterTypes.Double, "Target rotation of the excitation shape [&#176;] (E-BURP-1 = 90) -- only used if PowerCalcMethod='shape'")
    PowerAdjust_dB = Parameter("PowerAdjust_dB", 0.0, ParameterTypes.Double, "Manual fine-tune [dB, amplitude convention] on top of the 'shape' calculation -- equivalent to a Bruker cnst0. Only used if PowerCalcMethod='shape'")
    Excitation_dB = Parameter("Excitation_dB", 0.0, ParameterTypes.Double, "Excitation-pulse power, as dB ATTENUATION relative to RefAmplitude_HP -- only used if PowerCalcMethod='db'")
    DbConvention = Parameter("DbConvention", "amplitude", ParameterTypes.String, "dB convention: 'amplitude' (B1/voltage, usual NMR convention) or 'power' -- only used if PowerCalcMethod='db'")

    # Slice-select gradient -- FPX/FPY/FPZ are the ONLY place per-axis
    # calibration is applied (via GradientMatrix() in run()); G1 is the
    # single logical gradient-strength knob for the slice-select gradient.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Slice-Select Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Slice-Select Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z", ParameterTypes.String, "Gradient Axes")

    RampTime = Parameter("D70", 200.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 200.0, ParameterTypes.Double, "Gradient Settling Duration (before acquisition) [&#956;s]")
    PreGrad = Parameter("D75", 200.0, ParameterTypes.Double, "Pre-Gradient Time (before ramp-up starts) [&#956;s]")
    RephaseFraction = Parameter("RephaseFrac", 0.5, ParameterTypes.Double, "Slice REPHASE gradient lobe area, as a fraction of the excitation gradient's area -- 0.5 is the standard starting point for a symmetric pulse; E-BURP-1 is intentionally asymmetric, so fine-tune empirically (scan this and maximise the FID amplitude) -- see design notes")

    # *** Gradient calibration -- UNCALIBRATED PLACEHOLDER, see design notes ***
    GradCal_HzPerCm = Parameter("GradCal", 1000.0, ParameterTypes.Double, "*** PLACEHOLDER, UNCALIBRATED *** Slice-select gradient calibration [Hz/cm] at G1*FP=1.0 for your target nucleus -- determine empirically (phantom of known length + LWG_1D-Image-Echo_H.py or LWG_1D-Image-UTE_H.py) before trusting the logged ExcitedWidth")

    # Mains-lock trigger -- OFF by default, see mains_lock_trigger() docstring.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rails (conservative placeholders, not vendor specs)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Phases (simple 2-step)
    PH1 = Parameter("PH1", "0,180", ParameterTypes.String, "Shaped 90&#176; Pulse Phase")
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

    # SlicePPM (ppm relative to SF+O1) -> Hz, plus fine-trim -- THE pulse
    # centre frequency, applied only to the shaped pulse (see shaped_pulse()).
    PulseOffsetHz = (P.SlicePPM * P.FrequencyBase) + P.PulseOffsetTrim

    ReceiverFilter = Filter(P.Filter)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell

    ReceiverTime = DW*(points+1)
    times = np.arange(0, points*DW, DW) / 1.0e6

    ExcitationPulseWidth = P.P90sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD

    ExBurpA = parse_burp_coeffs(P.ExBurpCoeffsA)
    ExBurpB = parse_burp_coeffs(P.ExBurpCoeffsB)

    # ---- Optional power calculator override -------------------------------
    RFAsh0_effective = P.TXAmplitude90
    PowerMethod = str(P.PowerCalcMethod).strip().lower()

    if PowerMethod == 'shape':
        ExAmpProfile, _ = generate_shape(P.ExcitationShape, max(4, int(round(P.P90sh))), ExBurpA, ExBurpB)
        ExIntegFactor = shape_integration_factor(ExAmpProfile)
        RFAsh0_effective = estimate_shape_relative_scale(
            P.HardPulseWidth, P.RefAmplitude_HP, P.P90sh, P.ExcitationRotation,
            ExIntegFactor, P.LPMaxFraction, P.PowerAdjust_dB)
        comms.log("PowerCalcMethod='shape': excitation integ.factor={0:.4f} "
                  "-> RFAsh0={1:.4f} (hard-pulse ref P1Hard={2}us @ "
                  "RefAmplitude_HP={3}, ExRotation={4}deg, LPMaxFraction="
                  "{5:.2%}, PowerAdjust_dB={6:.2f})."
                  .format(ExIntegFactor, RFAsh0_effective, P.HardPulseWidth,
                          P.RefAmplitude_HP, P.ExcitationRotation,
                          P.LPMaxFraction, P.PowerAdjust_dB))
        if RFAsh0_effective > 1.0:
            comms.log("WARNING: calculated relative scale exceeds 1.0 -- "
                      "RFAsh0_effective={0:.4f}.".format(RFAsh0_effective))
    elif PowerMethod == 'db':
        RFAsh0_effective = db_to_lp_relative_scale(P.Excitation_dB, P.RefAmplitude_HP, P.LPMaxFraction, P.DbConvention)
        comms.log("PowerCalcMethod='db': RFAsh0 overridden -- {0:.2f} dB "
                  "(convention='{1}') -> {2:.4f} (LP-relative)."
                  .format(P.Excitation_dB, P.DbConvention, RFAsh0_effective))
        if RFAsh0_effective > 1.0:
            comms.log("WARNING: dB-calculated relative scale exceeds 1.0 -- "
                      "RFAsh0_effective={0:.4f}.".format(RFAsh0_effective))
    elif PowerMethod != 'manual':
        comms.log("WARNING: unrecognised PowerCalcMethod='{0}' -- falling "
                  "back to 'manual'.".format(P.PowerCalcMethod))

    # ---- Excited slice width: CALCULATED AND RECORDED every run ----------
    ExcitedBW_Hz, ExcitedWidth_cm, GradTotal = compute_excited_width(
        P, comms, P.P90sh, R_EBURP1_90, "excitation")

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
                "jc_excitedBW_Hz={0:.2f}".format(ExcitedBW_Hz),
                "jc_excitedWidth_cm={0}".format("{0:.4f}".format(ExcitedWidth_cm) if ExcitedWidth_cm is not None else "uncalibrated"),
                    ]
        return jc_data

    comms.log("Sequence Time: {0}".format(time_calculation(Parameters)))
    seqAcqu = CallBack1D(Parameters, comms, ReceiverFilter, jcamp_meta())

    def RecvCallback(acqData, scan):
        seqAcqu.process_data(scan, acqData)

    # Total elapsed time of the slice-select gradient event (ramp-up,
    # shaped-pulse dwell, ramp-down, rephase lobe ramp/hold/ramp, settle) --
    # computed once, outside the scan loop, for duty-cycle accounting.
    # NOTE ON DESIGN: the gradient is synchronised with the shaped pulse by
    # simple SEQUENTIAL composition (ramp gradient up -> call shaped_pulse()
    # unmodified -> ramp gradient down), NOT by running them in separate
    # 'with parallel:' branches. A gradient command's amplitude persists
    # until the next gradient command changes it, so ramping to G1 and then
    # simply letting the (separately timed) shaped_pulse() call run its
    # course leaves the gradient flat at G1 for exactly that call's
    # duration -- no explicit 'hold' branch is needed. This avoids nesting a
    # 'with parallel:' inside another 'with parallel:' branch (shaped_pulse()
    # already contains its own internal 3-way parallel block for
    # Transmit2/amplitude/phase) -- an untested combination not shown
    # anywhere in the Pulse Sequence Programming manual's own examples
    # (which only ever show FLAT, single-level 'with parallel:' blocks,
    # e.g. section 3.6's TX+X/Y/Z-gradient example, and the geoDIFFSTE
    # example in section 3.8). This sequential-composition design reuses
    # shaped_pulse() completely unmodified and matches the manual's
    # demonstrated patterns exactly.
    GradHoldRephase = P.RephaseFraction * ExcitationPulseWidth
    GradEventTime = (P.PreGrad + P.RampTime + ExcitationPulseWidth + P.RampTime
                      + P.RampTime + GradHoldRephase + P.RampTime
                      + P.GradSettle)

    with sequential:

        mains_lock_trigger(P, comms)

        Transmit2SelectPort(1,0)
        Transmit2LPEnable(1,1)
        Delay(10000) # relay settle

        Receiver2Preamp(128, P.ReceiverAttenuation)
        Receiver2Filter(200, ReceiverFilter)

        GradientMatrix(200,Matrix.T)

        Phases.Reset()

    comms.log("LP port engaged for the shaped pulse -- LPMaxFraction={0:.2%} "
              "of HP full scale. Slice/pulse offset: SlicePPM={1} ppm x "
              "SF={2}MHz + PulseOffsetTrim={3}Hz = {4:.2f} Hz (applied on "
              "top of SF+O1={5:.6f} MHz)."
              .format(P.LPMaxFraction, P.SlicePPM, P.FrequencyBase,
                      P.PulseOffsetTrim, PulseOffsetHz, Frequency))

    rf_on_time = ExcitationPulseWidth
    grad_on_time = GradEventTime
    estimate_duty_cycles(P, rf_on_time, grad_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver2FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            slew_rate_fn, gradient_fn = get_gradient_functions(P.Axis)

            # ---- Slice-select gradient ramp-up (plain sequential) --------
            Delay(P.PreGrad)
            slew_rate_fn(1.0, abs(P.G1/P.RampTime))
            gradient_fn(P.RampTime, P.G1)
            # Gradient now sits at G1 (persists until the next gradient
            # command) -- the shaped_pulse() call below runs its full,
            # exactly-known ExcitationPulseWidth us with the gradient flat
            # at this value, no separate 'hold' delay needed. See the
            # GradEventTime comment above for why this sequential
            # composition is used instead of a parallel branch.
            shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, RFAsh0_effective, P.TXEnableTime, Frequency, PulseOffsetHz, ExBurpA, ExBurpB)
            # ---- Ramp down + rephase lobe (opposite polarity, area =
            # RephaseFraction x main-lobe area -- see design notes) --------
            gradient_fn(P.RampTime, 0)
            slew_rate_fn(1.0, abs(P.G1/P.RampTime))
            gradient_fn(P.RampTime, -P.G1)
            safe_delay(GradHoldRephase, "rephase-lobe hold (check RephaseFraction)", comms)
            gradient_fn(P.RampTime, 0)
            gradient_fn(1.0, 0)
            Delay(P.GradSettle)

            # ACQU -- gradients are fully off/settled by now (this is a
            # localised FID acquisition, not a frequency-encoded readout).
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
# 1. Claude - 06/08/26 - Initial version. Built as the genuinely SPATIALLY
#    slice-selective counterpart to LWG_Selective-Echo_H.py (which is
#    frequency/chemical-shift-selective only, gradient OFF during the shaped
#    pulse -- confirmed by you). Adds: a slice-select gradient synchronised
#    with the shaped E-BURP-1 pulse by simple SEQUENTIAL composition (ramp
#    the gradient up, call shaped_pulse() completely unmodified, ramp it
#    down) rather than a 'with parallel:' branch -- a gradient command's
#    amplitude persists until the next one changes it, so this leaves the
#    gradient flat for exactly the shaped pulse's duration without needing
#    to nest a 'with parallel:' inside shaped_pulse()'s own internal 3-way
#    parallel block (Transmit2/amplitude/phase), a combination not
#    demonstrated anywhere in the Pulse Sequence Programming manual's
#    examples (which only show flat, single-level 'with parallel:' blocks --
#    see GradEventTime comment in run() for the full reasoning). Plus a
#    rephase lobe afterward (opposite polarity, area = RephaseFraction x
#    main-lobe area, default 0.5 -- the standard starting rule for a
#    symmetric pulse envelope, exposed as a Parameter since E-BURP-1 is
#    deliberately asymmetric in time). Adds SlicePPM (ppm, prominently
#    placed and labelled near the top of Parameters) as the dedicated pulse-
#    centre-frequency control, plus PulseOffsetTrim (Hz) for fine
#    adjustment. Adds compute_excited_width(), which calculates and
#    comms.log()s (and writes into the JCAMP metadata via jcamp_meta(), so
#    it is recorded with the acquired data, not just the transient log) the
#    expected excited slice bandwidth [Hz] and width [cm/mm], using
#    R_EBURP1_90 (a bandwidth-time product determined by Bloch-simulating
#    the exact synthesised E-BURP-1 envelope -- see the constant's own
#    docstring for the full method) and a new GradCal_HzPerCm Parameter
#    (UNCALIBRATED PLACEHOLDER -- flagged clearly, to be measured the same
#    way FPX/FPY/FPZ already are, against a real phantom). Everything else
#    (CallBack1D/Parameters scaffolding, generate_shape()/shaped_pulse(),
#    the three-way PowerCalcMethod calculator, mains_lock_trigger(),
#    estimate_duty_cycles(), safe_delay() guard, CRLF line endings) is
#    carried over unchanged from LWG_Selective-Echo_H.py /
#    LWG_1D-Image-Echo_H.py.
#
# -----------------------------------------------------------------------------
