#-------------------------------------------------------------------------------
# Name:        LWG_1D-Image-Echo-Selective_H.py
# Purpose:     One-dimensional MRI / image-echo (frequency-encoded profile)
#              sequence, EXCITED AND REFOCUSED WITH SHAPED, CHEMICAL-SHIFT-
#              SELECTIVE PULSES, so the resulting profile reflects only the
#              spatial distribution of ONE chosen resonance. X-Pulse
#              Broadband Benchtop NMR Spectrometer (1H/19F channel).
#              Axis-selectable (x/y/z).
#
# Author:      Claude, for L. Gordon (DTU), built on LWG_1D-Image-Echo_H.py
#
# Created:     04/08/2026
# Revised:     04/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     4.0
#
# IMPORTANT (v4.0): if your v3.0 profile was just 'burning a hole' at the
# target frequency instead of forming a clean selective image, the most
# likely cause was uncalibrated RF power on the wrong (high-power) TX port
# -- see changelog item 5 below, and calibrate P90sh/P180sh/RFAsh0/RFAsh1
# with the new LWG_Selective-Echo_H.py (no imaging, just a frequency-
# selective spin echo) before re-trying this sequence.
#
# Design notes:
#  - This is LWG_1D-Image-Echo_H.py (v3.0, itself derived from your proven
#    LWG_diffprof_2_H.py) with the hard 90 pulse and composite 90-270
#    refocusing pulse replaced by shaped excitation/refocusing pulses, so it
#    inherits the same branch-timing structure, dead-time handling,
#    safe_delay() guards, and duty-cycle warnings that are already working
#    for you.
#  - The read/phase gradients are OFF during both shaped pulses -- this is
#    CHEMICAL-SHIFT selectivity (pick one resonance out of a spectrum), not
#    spatial SLICE selectivity (which would require a gradient running
#    during the shaped pulse). If you actually want slice-selective
#    excitation (a gradient on during the shaped pulse, e.g. for true 2D/3D
#    slice-select MRI), say so -- that needs an extra gradient event
#    synchronised with the shaped pulse and isn't in this version.
#  - SHAPES ARE SYNTHESISED ON THE FLY IN PYTHON, not loaded from an
#    external .csv shape file. generate_shape()/shaped_pulse() below drive
#    amplitude+phase step-by-step (1 us/step) alongside a single continuous
#    Transmit1() call, exactly the technique already proven for the
#    WURST-20 pulse in your 'sequence elements.py' template -- just
#    modulating amplitude+phase instead of amplitude+frequency. This needs
#    nothing pre-loaded on the hardware and is fully reproducible from this
#    file alone.
#  - DEFAULT SHAPES ARE NOW THE REAL LITERATURE PULSES: excitation =
#    E-BURP-1, refocusing = RE-BURP (Geen, H. & Freeman, R., J. Magn.
#    Reson. 93, 93-141 (1991)), reconstructed from their published
#    truncated-Fourier-series coefficients (E-BURP-1: Table 2, nmax=8;
#    RE-BURP: Table 8, Np=256), transcribed from the two coefficient files
#    you provided. 'GAUSSIAN'/'SINC' are still available (no coefficients
#    needed, exact closed-form) as a simpler fallback if you ever want one.
#  - Selective pulses are long (ms, not us) to get useful selectivity on a
#    60 MHz benchtop -- default P90sh/P180sh/TAU below are a starting point,
#    not a calibrated value for any particular resonance separation; you
#    will need to tune pulse width (roughly bandwidth ~ 1/duration,
#    shape-dependent) against your actual spectrum.
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
            self.comms.log("seqTime LWG_1D-Image-Echo-Selective_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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
    # message).
    t_scanTime = (P.P90sh + P.P180sh + P.Tau*2 + points*DW)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("1D MRI / image-echo (frequency-encoded profile) sequence, "
                "with SHAPED, chemical-shift-selective excitation and "
                "refocusing pulses, on the {H/F} channel.\nSpin-echo "
                "(shaped-90-TAU-shaped-180-TAU-echo) with a dephase gradient "
                "before the refocusing pulse and a matched readout gradient "
                "spanning the acquisition window. Gradients are OFF during "
                "both shaped pulses (chemical-shift selectivity, not slice "
                "selectivity). Gradient axis selectable via GradAxis (x/y/z).")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,TAU,D70,D73,D71,D75,G1,GradAxis,P90sh,P180sh,PulseOffset"

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
    in run() below, NOT here."""
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
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps
    ON THE FLY -- no external .csv shape file needed. Driven into the pulse
    step-by-step (1 us/step) via Transmit1SetScale/Channel1SetBasePhase,
    exactly the technique already proven for the WURST-20 pulse in your
    'sequence elements.py' / FID_1D_WURST_HX.py (there it modulates
    amplitude+FREQUENCY; here it modulates amplitude+PHASE instead, which is
    the right pair of controls for an amplitude/phase-modulated selective
    pulse like a Gaussian, sinc, or BURP-family shape).

    shape_name options:
      'GAUSSIAN' -- truncated Gaussian (+/-2.5 sigma, <1% truncation),
                    constant phase. Exact closed-form, no external data
                    needed. Good general-purpose selective excitation --
                    similar behaviour to the vendor's 'GAUSS201'.
      'SINC'     -- truncated sinc x Hamming window (+/-5 main-lobe widths
                    each side). Sharper/more rectangular excitation profile
                    than a Gaussian. Exact closed-form.
      'EBURP1' / 'REBURP' / 'BURP' -- reconstructs a Geen & Freeman
                    (J. Magn. Reson. 93, 93-141, 1991) BURP-family pulse
                    from its truncated Fourier series (burp_coeffs_A /
                    burp_coeffs_B). Defaults are pre-populated from the
                    paper's Table 2 (E-BURP-1, nmax=8) and Table 8 (RE-BURP,
                    Np=256), transcribed from the two files you added
                    (ExBurpCoeffsA/B, RefBurpCoeffsA/B Parameters below) --
                    RE-BURP is purely amplitude/0-180-phase modulated in the
                    paper (no B_n terms), which is why RefBurpCoeffsB is
                    empty by default; that's expected, not an error.
    """
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
        # Only the cosine (A_n) coefficients are required -- a BURP pulse
        # with no B_n terms (like RE-BURP, see above) is a perfectly valid,
        # purely-real Fourier series, not a missing-data error.
        if not burp_coeffs_A:
            raise ValueError(
                "shape '{0}' requested but no Fourier A-coefficients were "
                "supplied (BurpCoeffsA is empty). EBURP1/REBURP are "
                "defined in Geen, H. & Freeman, R., J. Magn. Reson. 93, "
                "93-141 (1991) as a truncated Fourier series -- enter "
                "coefficients as ExBurpCoeffsA/B or RefBurpCoeffsA/B "
                "(comma-separated A0,A1,A2,... and B0,B1,B2,...), or use "
                "'GAUSSIAN' / 'SINC' instead, which need no coefficients."
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
# Exact fixed (non-parallel) instruction cost of shaped_pulse(), in us:
#   Channel1SetFrequency(1,...)  [pulse-offset ON]   = 1
#   Transmit1SetScale(5,...)                         = 5
#   Channel1SetBasePhase(5,...)                       = 5
#   Transmit1BlankingOn(1)                             = 1
#   Transmit1BlankingOff(1)                            = 1
#   Channel1SetFrequency(1,...)  [pulse-offset OFF]  = 1
#                                                      -----
#                                                        14
# (plus txenabletime and n_steps, added separately wherever this constant
# is used). Every hardware-instruction argument here is itself a real
# elapsed-time cost in this framework (confirmed against LWG_diffprof_2_H.py
# / ninety270()'s '+4' overhead, which is the same kind of sum) -- so this
# is an exact count, not a guess. Used consistently below (and in run())
# wherever shaped_pulse()'s true elapsed time needs to be known, so the
# echo stays centred on TAU as accurately as this framework allows.

def shaped_pulse(duration, phase, shape, amplitude, txenabletime,
                  base_frequency, pulse_offset,
                  burp_coeffs_A=None, burp_coeffs_B=None):
    """Apply a chemical-shift-selective shaped RF pulse, SYNTHESISED ON THE
    FLY in Python -- no external .csv shape file needed (see
    generate_shape() above). Mechanism: a single continuous Transmit1
    (duration) call runs in one 'with parallel:' branch, while two other
    branches step Transmit1SetScale (amplitude) and Channel1SetBasePhase
    (phase) once per microsecond -- identical in structure to the WURST-20
    pulse already proven in your 'sequence elements.py' template, just
    modulating amplitude+phase instead of amplitude+frequency.

    base_frequency / pulse_offset: the RF carrier is set to
    (base_frequency + pulse_offset) for the duration of THIS shaped pulse
    ONLY, then restored to base_frequency immediately afterward. This keeps
    acquisition referenced to a stable base_frequency (SF+O1) while the
    selective pulse itself targets an independently-tunable PulseOffset --
    exactly what you want for calibration: sweep PulseOffset to find/tune
    the pulse's frequency response without the acquired spectrum's
    frequency axis shifting every time you retune it.

    duration is in microseconds and is rounded to the nearest integer
    number of 1 us steps (n_steps); the shape is resolved at n_steps points,
    so duration should be at least a few tens of microseconds for a
    reasonable profile (your P90sh/P180sh defaults of 5000/10000 us give
    5000/10000 points, far more than enough)."""
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
    """Delay() wrapper for every timing-formula-derived (as opposed to
    fixed) delay in the sequence. If TAU/GradientOnTime/RampTime/PreGrad/
    P90sh/P180sh don't leave enough room, the raw formula would go negative
    -- which can silently desynchronise the read gradient from the
    acquisition window ('leaking' signal/artefacts outside the intended
    echo) instead of failing loudly. This turns that failure mode into a
    specific, immediate error instead."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase TAU and/or PreGrad, or decrease "
               "P90sh/P180sh/GradientOnTime/RampTime/GradSettle, then retry."
               ).format(label, value)
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

    OFF (UseMainsLock=0) by default for every sequence in this imaging
    family, per your instruction -- mains triggering adds a non-
    deterministic pre-scan delay that works against reproducible
    gradient-echo timing. Set UseMainsLock=1 (and confirm MainsLockChannel)
    only if you determine your X-Pulse needs mains-synchronous triggering.
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


def parse_burp_coeffs(s):
    """Parse a comma-separated 'A0,A1,A2,...' string Parameter into a list
    of floats. Returns [] for an empty/whitespace-only string (the safe
    default -- generate_shape() then raises a clear error if EBURP1/REBURP
    is selected without real coefficients, instead of silently using
    zeros)."""
    s = (s or "").strip()
    if not s:
        return []
    return [float(x) for x in s.split(",") if x.strip() != ""]


def estimate_duty_cycles(P, rf_on_time, grad_on_time, comms):
    """Log estimated RF and gradient duty cycles for this scan and warn if
    they exceed the (user-adjustable) MaxRFDuty / MaxGradDuty parameters.

    IMPORTANT: MaxRFDuty and MaxGradDuty are conservative, editable
    placeholders (5% / 10%), NOT a vendor-confirmed rating -- no public duty
    cycle specification for the X-Pulse 60 MHz RF/gradient amplifiers was
    available when this was written. Confirm the real limits with Oxford
    Instruments and adjust accordingly. Shaped pulses here are typically
    MUCH longer (ms) than hard pulses (us), so RF duty cycle is worth
    watching more closely on this sequence than on the hard-pulse version.
    """
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR
    grad_duty = grad_on_time / TR

    comms.log("Estimated duty cycle -- RF: {0:.3%} (limit {1:.1%}), "
              "Gradient [{2}-axis]: {3:.3%} (limit {4:.1%})".format(
              rf_duty, P.MaxRFDuty, P.Axis, grad_duty, P.MaxGradDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Shaped pulses here can be a large fraction of "
                  "RD -- consider a longer RD, shorter P90sh/P180sh, or lower "
                  "RFAsh0/RFAsh1, or confirm with Oxford Instruments that "
                  "this is within the transmitter's rated duty cycle before "
                  "running unattended.".format(rf_duty, P.MaxRFDuty))
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
    Sequence = Parameter("Sequence", "LWG_1D-Image-Echo-Selective_H", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 59.7, ParameterTypes.Double, "H/F Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "H/F Freq Offset [Hz] -- overall sequence reference (affects BOTH the pulses and acquisition); leave at 0 (or your usual reference) and use PulseOffset below to target the selective pulses at a specific resonance")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "H/F TX Freq Offset [ppm]")
    # THE clear, dedicated place to set the selective pulses' target
    # frequency: applied ONLY during the shaped pulses (see shaped_pulse()),
    # NOT to acquisition, so the acquired spectrum's frequency axis stays
    # fixed at SF+O1 while you tune/calibrate this independently.
    PulseOffset = Parameter("PulseOffset", 0.0, ParameterTypes.Double, "Selective-Pulse Frequency Offset [Hz] -- targets ONLY the shaped 90/180 pulses at (SF+O1+PulseOffset); acquisition stays referenced to SF+O1. THIS is the parameter to sweep/set when calibrating or picking which resonance to image.")

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

    # Shaped (chemical-shift-selective) pulses -- replace the hard 90 /
    # composite 90-270 refocusing pulse used in LWG_1D-Image-Echo_H.py.
    # Shapes are SYNTHESISED ON THE FLY (see generate_shape()/shaped_pulse()
    # above) -- no external shape file needed. Now default to the real
    # literature pulses: EXCITATION = EBURP1 (Geen & Freeman 1991, Table 2,
    # nmax=8), REFOCUSING = REBURP (Table 8, Np=256) -- coefficients
    # transcribed from the two files you provided. 'GAUSSIAN'/'SINC' remain
    # available (no coefficients needed) if you want a simpler fallback.
    ExcitationShape = Parameter("ExShape", "EBURP1", ParameterTypes.String, "Excitation (90) Shape [EBURP1 (default), GAUSSIAN, SINC, or BURP -- see ExBurpCoeffsA/B]")
    RefocusShape = Parameter("RefShape", "REBURP", ParameterTypes.String, "Refocusing (180) Shape [REBURP (default), GAUSSIAN, SINC, or BURP -- see RefBurpCoeffsA/B]")
    # E-BURP-1 (excitation), Geen & Freeman 1991, Table 2, nmax=8 column.
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00", ParameterTypes.String, "E-BURP-1 Fourier cosine coeffs A0..A8 (Geen & Freeman 1991, Table 2, nmax=8) -- used when ExShape=EBURP1/BURP")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01", ParameterTypes.String, "E-BURP-1 Fourier sine coeffs B0(unused)..B8 (Geen & Freeman 1991, Table 2, nmax=8) -- used when ExShape=EBURP1/BURP")
    # RE-BURP (refocusing), Geen & Freeman 1991, Table 8, Np=256 column --
    # purely real/cosine (no B_n terms; RE-BURP is amplitude/0-180-phase
    # modulated only), hence RefBurpCoeffsB is empty by default.
    RefBurpCoeffsA = Parameter("RefBurpCoeffsA", "0.49,-1.02,1.11,-1.57,0.83,-0.42,0.26,-0.16,0.10,-0.07,0.04,-0.03,0.01,-0.02,0.00,-0.01", ParameterTypes.String, "RE-BURP Fourier cosine coeffs A0..A15 (Geen & Freeman 1991, Table 8, Np=256) -- used when RefShape=REBURP/BURP")
    RefBurpCoeffsB = Parameter("RefBurpCoeffsB", "", ParameterTypes.String, "RE-BURP Fourier sine coeffs (none published -- RE-BURP is purely real; leave empty) -- used when RefShape=REBURP/BURP")
    P90sh = Parameter("P90sh", 5000.0, ParameterTypes.Double, "Shaped 90&#176; Pulse Width [&#956;s] -- also sets the shape resolution (1 point/&#956;s)")
    P180sh = Parameter("P180sh", 10000.0, ParameterTypes.Double, "Shaped 180&#176; Pulse Width [&#956;s] -- also sets the shape resolution (1 point/&#956;s)")
    # UNCALIBRATED starting values -- these now scale power on the LOW-POWER
    # (LP) port (see run(), 'Switch to LOW-POWER TX port'), not the
    # high-power port these numbers were previously (also uncalibrated)
    # scaling. LP's own full scale (1.0) is only ~LPMaxFraction of the HP
    # port's full scale, per your measurement -- so 0.30 here is a modest,
    # cautious starting guess, NOT a calibrated flip angle. Use
    # LWG_Selective-Echo_H.py (the calibration-only pp, no imaging/
    # gradients) to find the real 90/180 condition via a nutation sweep
    # before trusting quantitative results from this sequence. Using the
    # wrong power here was almost certainly why the selective image just
    # 'burned a hole' rather than forming a clean echo -- with a badly
    # wrong flip angle the pulses don't coherently refocus, so instead of a
    # clean spin echo you mostly destroy/saturate signal at that frequency.
    TXAmplitude90 = Parameter("RFAsh0", 0.30, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0.0&#8230;1.0 of LP max] -- UNCALIBRATED, see LWG_Selective-Echo_H.py", RFA, min=0.0, max=1.0)
    TXAmplitude180 = Parameter("RFAsh1", 0.30, ParameterTypes.Double,
                               "Shaped 180&#176; TX Power [0.0&#8230;1.0 of LP max] -- UNCALIBRATED, see LWG_Selective-Echo_H.py", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", 0.10, ParameterTypes.Double, "LP Port Max Power as Fraction of HP Port [0.0&#8230;1.0] -- your measured value, for reference/logging only; re-measure and update if it drifts")

    # Gradients -- FPX/FPY/FPZ are the ONLY place per-axis calibration is
    # applied (via GradientMatrix() in run()); G1 is the single logical
    # gradient-strength knob used for both encode and readout lobes.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Gradient Strength [-1.0&#8230;1.0]")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")

    # Sequence specific
    GradientOnTime = Parameter("D71", 1000.0, ParameterTypes.Double, "Gradient Duration [&#956;s]")
    RampTime = Parameter("D70", 100.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSettle = Parameter("D73", 100.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s]")
    Tau = Parameter("TAU", 20000, ParameterTypes.Int32, "Echo Time [&#956;s] -- must comfortably exceed P90sh/2 + P180sh (shaped pulses are long)")
    PreGrad = Parameter("D75", 100.0, ParameterTypes.Double, "Pre-Gradient Time [&#956;s]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")

    # Mains-lock trigger -- OFF by default for imaging (see
    # mains_lock_trigger() docstring for why). Channel is unconfirmed for
    # X-Pulse specifically.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock ExternalTrigger Channel [1-3, unconfirmed for X-Pulse -- see manual 3.7.19]")

    # Duty-cycle guard rails (see estimate_duty_cycles() docstring: these are
    # conservative, user-adjustable placeholders, not vendor-confirmed specs)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Phases (simple 2-step, since there's no composite-pulse phase cycle
    # to account for any more)
    PH1 = Parameter("PH1", "0,180", ParameterTypes.String, "Shaped 90&#176; Pulse Phase")
    PH2 = Parameter("PH2", "0", ParameterTypes.String, "Shaped 180&#176; Pulse Phase")
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
    Channel1SetFrequency(10,Frequency)
    Channel1RestartSynth(10)

    ReceiverFilter = Filter(P.Filter)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell

    # Total acquisition-window duration -- Receiver1()'s 'duration' argument
    # is the TOTAL window (points*dwell), NOT the per-point dwell.
    ReceiverTime = DW*(points+1)

    times = np.arange(0, points*DW, DW) / 1.0e6

    # Exact elapsed time of a shaped_pulse() call for the excitation/
    # refocusing pulses respectively -- see SHAPED_PULSE_FIXED_OVERHEAD
    # above for the derivation (this is an exact instruction-cost sum, not
    # an estimate). Used consistently in every surrounding wait so the echo
    # stays centred on TAU.
    ExcitationPulseWidth = P.P90sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD
    RefocusPulseWidth = P.P180sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD

    # Parsed once, outside the scan loop -- see generate_shape()/
    # parse_burp_coeffs() above. Excitation and refocusing pulses use
    # DIFFERENT published coefficient sets (E-BURP-1 vs RE-BURP), so they
    # are kept as separate pairs and each is only used if the matching
    # ExShape/RefShape is actually set to EBURP1/REBURP/BURP.
    ExBurpA = parse_burp_coeffs(P.ExBurpCoeffsA)
    ExBurpB = parse_burp_coeffs(P.ExBurpCoeffsB)
    RefBurpA = parse_burp_coeffs(P.RefBurpCoeffsA)
    RefBurpB = parse_burp_coeffs(P.RefBurpCoeffsB)

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

        # Switch to the LOW-POWER (LP) TX port for these shaped pulses --
        # per your measurement, LP tops out at ~LPMaxFraction (default 10%)
        # of the high-power port's full scale, and the vendor's own long
        # shaped-pulse sequence (WET-FID_H.py) uses LP for exactly this
        # reason. This sequence's RF is 100% shaped pulses (no hard pulses
        # to switch back for), so LP is engaged once here, for the whole
        # experiment. Using the wrong (high-power) port at an uncalibrated
        # low relative-scale is the most likely reason the selective image
        # was 'burning a hole' instead of forming a clean echo -- see the
        # RFAsh0/RFAsh1 Parameter comments and design notes for more.
        Transmit1SelectPort(1,0)
        Transmit1LPEnable(1,1)
        Delay(10000) # relay settle -- same caution/value as the HP-path
                      # settle used elsewhere in this pp family; not
                      # vendor-confirmed as the true minimum for this
                      # specific transition, but this is a one-off setup
                      # cost, so kept conservative.

        Receiver1Preamp(128, P.ReceiverAttenuation)
        Receiver1Filter(200, ReceiverFilter)

        # Gradient Setup -- the ONLY place FPX/FPY/FPZ scale the gradients
        GradientMatrix(200,Matrix.T)

        Phases.Reset()

    comms.log("LP port engaged for shaped pulses -- LPMaxFraction={0:.2%} "
              "of HP full scale (your measured value; RFAsh0/RFAsh1 are "
              "relative to LP's OWN max, i.e. LPMaxFraction*HP). Pulse "
              "offset for the selective pulses: PulseOffset={1} Hz "
              "(applied on top of SF+O1={2:.6f} MHz)."
              .format(P.LPMaxFraction, P.PulseOffset, Frequency))

    # ---- Duty-cycle warnings (computed once, outside the scan loop) ------
    rf_on_time = ExcitationPulseWidth + RefocusPulseWidth
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
                    # 90 delay -- mirrors the exact elapsed time of the
                    # shaped excitation pulse block below
                    Delay(ExcitationPulseWidth)
                    # Encode Gradient
                    apply_gradient(P.Axis, P.GradientOnTime, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                    # TAU - gradient time
                    safe_delay(P.Tau - (P.GradientOnTime + 2*P.RampTime + P.GradSettle + 2) - RefocusPulseWidth/2.0,
                               "gradient first-TAU wait", comms) # first tau delay
                    # 180 delay -- park here for exactly as long as the
                    # shaped refocusing pulse block takes
                    safe_delay(RefocusPulseWidth,
                               "gradient refocus-width wait", comms)
                    # TAU
                    safe_delay(P.Tau-P.PreGrad, "gradient second-TAU wait", comms) # second tau delay
                    # Read Gradient
                    apply_gradient(P.Axis, P.GradientOnTime*2+P.RampTime*2, P.RampTime, P.G1, P.GradSettle, P.PreGrad)
                with sequential:
                    # Shaped 90 (excitation) pulse -- gradients are off here,
                    # so this selects a resonance (chemical shift), not a
                    # spatial slice. PulseOffset targets this pulse (and the
                    # refocusing pulse below) at the desired resonance,
                    # independent of the sequence-wide O1.
                    shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, P.TXAmplitude90, P.TXEnableTime, Frequency, P.PulseOffset, ExBurpA, ExBurpB)
                    # TAU
                    safe_delay(P.Tau - ExcitationPulseWidth/2 - RefocusPulseWidth/2.0,
                               "RF first-TAU wait", comms) # first tau delay
                    # Shaped 180 (refocusing) pulse
                    shaped_pulse(P.P180sh, ph["PH2"], P.RefocusShape, P.TXAmplitude180, P.TXEnableTime, Frequency, P.PulseOffset, RefBurpA, RefBurpB)
                    # TAU -- ReceiverFilter.dead_time is reserved out of this
                    # wait and paid back explicitly (as its own Delay, right
                    # before Receiver1 below) rather than left out entirely.
                    safe_delay(P.Tau - RefocusPulseWidth/2.0 + 1.*P.PreGrad - ReceiverFilter.dead_time,
                               "RF second-TAU wait", comms)
                    # ACQU
                    Channel1SetBasePhase(P.TXEnableTime,0)
                    Receiver1Phase(P.TXEnableTime, ph["PHRX"])
                    # Dead1: probe ring-down time.
                    Delay(P.Dead1)
                    # ReceiverFilter.dead_time: digital-filter settling time
                    # (varies with the Filter bandwidth you choose) -- see
                    # LWG_1D-Image-Echo_H.py changelog item (h) for why this
                    # matters.
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
# 1. Claude - 04/08/26 (v1.0, SUPERSEDED) - Initial version, adapted from
#    LWG_1D-Image-Echo_H.py (v3.0) by replacing the hard 90 pulse and
#    composite 90-270 refocusing pulse with independent shaped
#    (Transmit1ShapedPulse) excitation and refocusing pulses, using the
#    shaped-pulse pattern proven in the vendor's WET-FID_H.py (Default
#    pps). This version required an external .csv shape file ('GAUSS201')
#    to already exist on the X-Pulse hardware computer.
# 2. Claude - 04/08/26 (v2.0, THIS VERSION) - ON-THE-FLY SHAPE SYNTHESIS:
#    you confirmed you don't have eBURP1/REBURP shape files, so
#    Transmit1ShapedPulse() + external .csv is replaced with
#    generate_shape()/shaped_pulse(), which synthesise the RF envelope in
#    Python and drive it into the pulse step-by-step (1 us/step) via
#    Transmit1SetScale (amplitude) + Channel1SetBasePhase (phase) run in
#    parallel with a single continuous Transmit1(duration) call -- the same
#    technique already proven for the WURST-20 pulse in your 'sequence
#    elements.py' template (there: amplitude+frequency; here:
#    amplitude+phase). 'GAUSSIAN' (new default for both ExShape/RefShape,
#    replacing 'GAUSS201') and 'SINC' are exact, closed-form, ready to use.
#    'EBURP1'/'REBURP'/'BURP' reconstruct a genuine Geen & Freeman (1991)
#    Fourier-series pulse IF you supply BurpCoeffsA/BurpCoeffsB (new string
#    Parameters, empty by default) -- I deliberately did NOT hard-code
#    published BURP coefficients from memory, since I don't have a verified
#    copy of the original table and didn't want to risk silently shipping
#    wrong numbers for something that drives real RF hardware;
#    generate_shape() raises a clear, specific error instead if you select
#    EBURP1/REBURP without filling these in. Also added mains-lock support
#    (see item 3) for consistency with the rest of the imaging family.
#    Everything else (gradient encode/decode structure, apply_gradient(),
#    safe_delay() guards, dead-time handling, duty-cycle warnings, CRLF
#    line endings) is carried over unchanged. See the design notes at the
#    top of this file for the chemical-shift-vs-slice-selectivity
#    distinction.
# 3. Claude - 04/08/26 - MAINS LOCK: added UseMainsLock/MainsLockChannel
#    parameters and mains_lock_trigger(), called at the very start of the
#    first 'with sequential:' block. Defaults to OFF (UseMainsLock=0), per
#    your instruction -- see mains_lock_trigger() docstring for the
#    reasoning. MainsLockChannel is a Parameter, not hard-coded, since the
#    manual doesn't list an X-Pulse row for the mains-lock trigger channel.
# 4. Claude - 04/08/26 (v3.0, THIS VERSION) - REAL EBURP1/REBURP
#    COEFFICIENTS: you provided the actual Geen & Freeman (1991)
#    coefficient tables (Table 2, E-BURP-1; Table 8, RE-BURP), so:
#      a) Fixed a validation bug in generate_shape(): it previously
#         required BOTH BurpCoeffsA and BurpCoeffsB to be non-empty before
#         reconstructing an EBURP1/REBURP/BURP shape. RE-BURP is purely
#         real (amplitude/0-180-phase only, no B_n terms in the published
#         table) -- requiring a non-empty B array made it impossible to
#         ever select REBURP correctly. Now only the A (cosine)
#         coefficients are required; missing/empty B is correctly treated
#         as all-zero.
#      b) Split the single BurpCoeffsA/BurpCoeffsB pair into two pairs --
#         ExBurpCoeffsA/ExBurpCoeffsB (used when ExcitationShape is
#         EBURP1/BURP) and RefBurpCoeffsA/RefBurpCoeffsB (used when
#         RefocusShape is REBURP/BURP) -- since E-BURP-1 and RE-BURP are
#         different pulses with different coefficients; the previous
#         single shared pair could not represent both at once.
#      c) Pre-populated real defaults: ExBurpCoeffsA/B from Table 2's
#         nmax=8 column (the most complete fit given in the paper),
#         RefBurpCoeffsA from Table 8's Np=256 column (RefBurpCoeffsB left
#         empty -- correctly so, see (a)). ExcitationShape/RefocusShape now
#         default to 'EBURP1'/'REBURP' (previously 'GAUSSIAN'/'GAUSSIAN',
#         used as a stand-in while coefficients were unavailable).
#      d) Sanity-checked the reconstructed envelopes numerically before
#         shipping: RE-BURP reconstructs as expected for a refocusing
#         pulse (symmetric about the pulse centre); E-BURP-1 reconstructs
#         asymmetric, matching its known literature behaviour (E-BURP-1 is
#         intentionally asymmetric so magnetization outside the excited
#         band returns to +z rather than to some other net phase).
#         'GAUSSIAN'/'SINC' remain available as a simpler, coefficient-free
#         fallback if you ever want one.
# 5. Claude - 04/08/26 (v4.0, THIS VERSION) - You reported the selective
#    image just 'burned a hole' at the target frequency rather than forming
#    a clean image. Diagnosis/fixes:
#      a) LOW-POWER PORT: v1.0-v3.0 left the transmitter on the HIGH-power
#         port (set once, at the top, for the whole sequence) and applied
#         RFAsh0/RFAsh1=0.10 there -- i.e. an UNCALIBRATED 10% of full HIGH
#         power for a multi-millisecond shaped pulse. You confirmed the
#         selective pulses should instead use the LOW-power port (which
#         tops out at ~10% of HP), matching the vendor's own long-shaped-
#         pulse sequence (WET-FID_H.py). Added the LP switch in run()'s
#         setup block (new LPMaxFraction Parameter documents the ~10%
#         relationship). This is the most likely root cause of the
#         hole-burning: at the wrong power on the wrong port, the pulses
#         almost certainly delivered a badly wrong flip angle, so instead
#         of a coherent 90-180 spin echo you got something closer to
#         partial saturation/dephasing at that frequency -- which looks
#         exactly like a 'hole' rather than a clean, refocused signal.
#      b) NEW PulseOffset PARAMETER: added a dedicated Hz offset, applied
#         ONLY during the two shaped pulses (via Channel1SetFrequency,
#         set before each shaped pulse and restored immediately after --
#         see SHAPED_PULSE_FIXED_OVERHEAD/shaped_pulse()), separate from
#         the sequence-wide O1. Acquisition stays referenced to SF+O1;
#         PulseOffset is the one clear, dedicated place to target/calibrate
#         which resonance the selective pulses act on, per your request
#         that all frequency-selective pps have this. Added to
#         sequence_basic() so it's prominent in the parameter list.
#      c) EXACT PULSE-WIDTH ACCOUNTING: replaced the previous approximate
#         '+6'/'+5' fixed-overhead guesses (used when estimating a
#         shaped_pulse() call's elapsed time, for centring the echo on TAU)
#         with an exact instruction-cost sum, SHAPED_PULSE_FIXED_OVERHEAD=14
#         (verified against the same accounting method used successfully in
#         LWG_1D-Image-Echo_H.py's ninety270()/RefocusPulseWidth). Also
#         confirmed numerically that 'with parallel:' branches do NOT need
#         to be exactly equal-duration for this framework to work (the
#         block's real time is the max of its branches) -- so this was a
#         precision improvement, not a fix for a compile-breaking bug.
#      d) RFAsh0/RFAsh1 bumped from 0.10 to 0.30 as a more reasonable (but
#         still explicitly UNCALIBRATED) LP-relative starting point --
#         calibrate for real with LWG_Selective-Echo_H.py before trusting
#         quantitative flip angles.
#
# -----------------------------------------------------------------------------
