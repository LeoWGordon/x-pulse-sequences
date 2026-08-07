#-------------------------------------------------------------------------------
# Name:        LWG_Selective-PulseAcquire_H.py
# Purpose:     Frequency-selective PULSE-ACQUIRE (shaped-90 then acquire, NO
#              refocusing pulse, NO echo, NO gradients). The most basic
#              calibration step for a selective pulse: find P90sh/RFAsh0/
#              PulseOffset in isolation, before adding a refocusing pulse
#              (LWG_Selective-Echo_H.py) or imaging gradients
#              (LWG_1D-Image-Echo-Selective_H.py). X-Pulse Broadband
#              Benchtop NMR Spectrometer (1H/19F channel).
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     06/08/2026
# Revised:     07/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.2
#
# Design notes:
#  - This is LWG_Selective-Echo_H.py with the refocusing pulse, TAU-TAU
#    echo timing, and every Ref*/P180sh/RFAsh1 Parameter removed -- shaped-
#    90 pulse, short fixed dead time, acquire. Simplest possible way to
#    find the EXCITATION condition (P90sh/RFAsh0/PulseOffset) on a real
#    resonance: sweep RFAsh0 (or P90sh) and watch the acquired signal
#    maximise. Do this FIRST, then move to LWG_Selective-Echo_H.py to also
#    calibrate the refocusing pulse using the same PulseOffset/P90sh/RFAsh0.
#  - Same generate_shape()/shaped_pulse() on-the-fly EBURP/GAUSSIAN/SINC
#    synthesis, LOW-power port, dedicated PulseOffset Parameter as the rest
#    of this pulse-programme family -- carry calibrated values straight
#    across between all of them.
#  - Parameter descriptions are kept SHORT here (units/critical info first)
#    since SpinFlow's parameter panel has limited display width -- see the
#    longer prose in LWG_Selective-Echo_H.py's comments/docstrings if you
#    need the full rationale for any of these.
#  - POWER: RFAsh0's default is a ONE-TIME approximate calculation from a
#    calibrated hard 90 (see _HARD_P90_WIDTH_US/_HARD_P90_AMPLITUDE/
#    _DEFAULT_RFASH0 above the Parameters block) -- not a live/reselectable
#    calculation. Checked against Oxford's own WET-FID_H.py/WET-FID-AUTO_H.py
#    (Default pps folder): those use a directly-calibrated shaped-pulse
#    power Parameter with NO hard-pulse-derived formula at all (their "AUTO"
#    variant only auto-detects the shaped pulse's FREQUENCY via a pilot
#    scan, never its power) -- so there's no vendor formalism to defer to
#    here; RFAsh0 is meant to be swept/recalibrated by hand regardless of
#    what its starting default is.
#  - SHAPE: default is E-BURP-1 pending real E-BURP-2 Fourier coefficients
#    (or a Bruker shape-file export) from Leo -- generate_shape() already
#    accepts 'EBURP2' as a shape name (identical Fourier-series code path
#    to EBURP1), so switching is a one-line default change once the real
#    numbers are in hand (see _DEFAULT_SHAPE_NAME above).
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
            self.comms.log("seqTime LWG_Selective-PulseAcquire_H Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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

    t_scanTime = (P.P90sh + points*DW)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Frequency-selective PULSE-ACQUIRE (shaped-90 then acquire), "
                "NO refocusing pulse/echo/gradients, on the {H/F} channel. "
                "Most basic calibration tool in this family: use PulseOffset "
                "to target a resonance and sweep P90sh/RFAsh0 to find the "
                "excitation 90 condition before adding a refocusing pulse "
                "(LWG_Selective-Echo_H.py) or imaging gradients.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,P90sh,RFAsh0,PulseOffset"

    return basic

def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical to LWG_Selective-Echo_H.py's version; duplicated
    here so this file stays fully self-contained.

    shape_name: 'GAUSSIAN', 'SINC', or 'EBURP1'/'EBURP2'/'BURP' (needs
    burp_coeffs_A -- see ExBurpCoeffsA/B Parameters)."""
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

    elif shape_name in ('EBURP1', 'EBURP2', 'BURP'):
        if not burp_coeffs_A:
            raise ValueError(
                "shape '{0}' requested but no Fourier A-coefficients were "
                "supplied. EBURP1/EBURP2 are defined in Geen, H. &amp; "
                "Freeman, R., J. Magn. Reson. 93, 93-141 (1991) as a "
                "truncated Fourier series -- enter coefficients as "
                "ExBurpCoeffsA/B, or use 'GAUSSIAN' / 'SINC' instead, which "
                "need no coefficients."
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
                          "SINC, EBURP1/EBURP2/BURP (needs BurpCoeffsA/B)."
                          .format(shape_name))

    peak = np.max(np.abs(amp))
    if peak > 0:
        amp = amp/peak
    return amp, phase


SHAPED_PULSE_FIXED_OVERHEAD = 14
# See LWG_Selective-Echo_H.py's shaped_pulse() for the exact
# instruction-cost derivation of this constant (identical here).

def shaped_pulse(duration, phase, shape, amplitude, txenabletime,
                  base_frequency, pulse_offset,
                  burp_coeffs_A=None, burp_coeffs_B=None):
    """Apply a chemical-shift-selective shaped RF pulse, synthesised on the
    fly -- identical mechanism to LWG_Selective-Echo_H.py's shaped_pulse()."""
    n_steps = int(round(duration))
    if n_steps < 4:
        raise ValueError("shaped_pulse duration too short for step-wise "
                          "synthesis (need >=4 us so the shape has enough "
                          "points to be meaningful): got {0} us"
                          .format(duration))
    amp_profile, phase_profile = generate_shape(shape, n_steps, burp_coeffs_A, burp_coeffs_B)

    # pulse_offset is in Hz (see PulseOffset Parameter); base_frequency
    # is in MHz -- *1.0e-6 converts Hz->MHz before adding. Omitting this
    # conversion sends the synth a frequency off by ~1e6x whenever
    # PulseOffset != 0, which fails hardware init (a real bug fixed
    # 06/08/2026 -- it was invisible with PulseOffset=0, the default).
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


def mains_lock_trigger(P, comms):
    """Optionally emit a mains-line trigger before the first pulse event.
    OFF by default -- see LWG_Selective-Echo_H.py for the full rationale."""
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


def shape_integration_factor(amp_profile, phase_profile_deg):
    """COHERENT (signed) average/peak ratio of a synthesised shape envelope
    (a.k.a. Bp/B1, or Bruker's 'Integ. Factor'): a shape spending most of
    its duration at low amplitude needs a higher PEAK amplitude than a
    rectangular pulse to deliver the same net rotation, since on-resonance
    rotation depends on the AVERAGE B1 over the pulse, not the peak.

    MUST use the SIGNED profile (amp*cos(radians(phase)), reconstructing
    the real-valued Fourier-series b1(t) that generate_shape() folds into
    amp=|b1|/phase=0-or-180), NOT mean(|amp_profile|) -- BURP-family shapes
    (and SINC) have genuine negative lobes that partially CANCEL in the
    on-resonance rotation (all sub-pulses share the same fixed axis when
    phase is 0/180 and offset=0, so the rotations commute and simply add
    algebraically). Using mean(|amp|) instead double-counts those lobes as
    if they always added constructively, overestimating the effective
    average B1 by a large factor for shapes with substantial negative
    lobes (EBURP1: ~4.9x too high -- confirmed 07/08/26 by Bloch-simulating
    the resulting 'default' RFAsh0 in leonmr/bloch_sim.py and finding it
    excited barely 18 degrees on-resonance instead of the intended 90).

    Used once, at import time below, to turn a calibrated hard-90 into an
    approximate starting RFAsh0 default -- see the comment above
    TXAmplitude90's Parameter() call. VERIFY any shape/power combination
    with leonmr/bloch_sim.py (or a real nutation curve) before trusting it
    quantitatively -- this is still only an approximation.
    """
    peak = np.max(np.abs(amp_profile))
    if peak <= 0:
        raise ValueError("shape_integration_factor: shape has zero peak amplitude")
    signed_profile = np.asarray(amp_profile) * np.cos(np.deg2rad(phase_profile_deg))
    return float(np.mean(signed_profile) / peak)


def estimate_rf_duty_cycle(P, rf_on_time, comms):
    """Log estimated RF duty cycle and warn if it exceeds MaxRFDuty
    (conservative, editable placeholder -- confirm with Oxford Instruments)."""
    TR = float(P.RecycleDelay)
    rf_duty = rf_on_time / TR

    comms.log("Estimated RF duty cycle: {0:.3%} (limit {1:.1%})".format(rf_duty, P.MaxRFDuty))

    if rf_duty > P.MaxRFDuty:
        comms.log("WARNING: estimated RF duty cycle {0:.2%} exceeds MaxRFDuty "
                  "({1:.2%}). Consider a longer RD, shorter P90sh, or lower "
                  "RFAsh0, or confirm with Oxford Instruments that this is "
                  "within the transmitter's rated duty cycle before running "
                  "unattended.".format(rf_duty, P.MaxRFDuty))


# ---- Approximate default RFAsh0, calculated ONCE (at import time, not per
# run) from a calibrated hard 90 pulse on the HP port -- replaces the
# PowerCalcMethod machinery this file used to have (manual/shape/db
# selector, RefAmplitude_HP/HardPulseWidth/ExRotation/PowerAdjust_dB/
# Excitation_dB/DbConvention Parameters), simplified out 07/08/26 as more
# complexity than a calibration tool needs: RFAsh0 is meant to be swept by
# hand anyway (that's this file's whole purpose), so a single reasonable
# STARTING VALUE is enough -- it doesn't need a live, re-selectable
# calculation. ALWAYS verify with a real nutation curve before trusting
# this quantitatively; it doesn't know your actual coil/sample.
#
#   RFAsh0 (LP-relative) ~= HardAmplitude * HardWidth
#                            / (ShapeWidth * IntegrationFactor) / LPMaxFraction
#
# (target rotation is 90 deg for both the hard reference and this
# excitation shape, so the two 90s in the full Bruker-style formula --
# see LWG_Selective-Echo_H.py's changelog for the original derivation --
# cancel out.)
_HARD_P90_WIDTH_US = 9.58          # this repo's usual calibrated hard-90 width (P90/P1Hard)
_HARD_P90_AMPLITUDE = 0.4          # ...at this HP-relative amplitude (RFA0)
_LP_MAX_FRACTION = 0.10            # measured LP-vs-HP max-power ratio (LPMaxFraction default, below)
_DEFAULT_SHAPE_NAME = 'EBURP1'     # TODO: switch to 'EBURP2' once real coefficients are supplied
_DEFAULT_SHAPE_WIDTH_US = 5000.0   # matches P90sh's own default, below
_DEFAULT_EXBURP_A = "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00"
_DEFAULT_EXBURP_B = "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01"

_default_amp_profile, _default_phase_profile = generate_shape(
    _DEFAULT_SHAPE_NAME, int(_DEFAULT_SHAPE_WIDTH_US),
    parse_burp_coeffs(_DEFAULT_EXBURP_A), parse_burp_coeffs(_DEFAULT_EXBURP_B))
_DEFAULT_RFASH0 = round(
    _HARD_P90_AMPLITUDE * _HARD_P90_WIDTH_US
    / (_DEFAULT_SHAPE_WIDTH_US * shape_integration_factor(_default_amp_profile, _default_phase_profile))
    / _LP_MAX_FRACTION, 4)


@ParameterBlock
class Parameters:

	# Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_Selective-PulseAcquire_H", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBase = Parameter("SF", 59.7, ParameterTypes.Double, "H/F Base Freq [MHz]")
    FrequencyOffset = Parameter("O1", 0.0, ParameterTypes.Double, "H/F Freq Offset [Hz] -- ref. for pulse+acq; use PulseOffset to target the shaped pulse instead")
    TxPPM = Parameter("TxPPM", 0.0, ParameterTypes.Double, "H/F TX Freq Offset [ppm]")
    PulseOffset = Parameter("PulseOffset", 0.0, ParameterTypes.Double, "Shaped-Pulse Freq Offset [Hz] from SF+O1 -- SWEEP to find your target resonance")

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

    # Shaped (chemical-shift-selective) excitation pulse -- ON THE FLY, no
    # shape file needed. Default shape is set by _DEFAULT_SHAPE_NAME above
    # (currently EBURP1 -- switch to EBURP2 there once real coefficients
    # are available; ExBurpCoeffsA/B below will need updating to match).
    ExcitationShape = Parameter("ExShape", _DEFAULT_SHAPE_NAME, ParameterTypes.String, "Excitation Shape [EBURP1/EBURP2(default)/GAUSSIAN/SINC/BURP]")
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", _DEFAULT_EXBURP_A, ParameterTypes.String, "E-BURP cosine coeffs A0..A8 (Geen&amp;Freeman'91) -- used if ExShape=EBURP1/EBURP2/BURP")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", _DEFAULT_EXBURP_B, ParameterTypes.String, "E-BURP sine coeffs B0..B8 (Geen&amp;Freeman'91) -- used if ExShape=EBURP1/EBURP2/BURP")
    P90sh = Parameter("P90sh", _DEFAULT_SHAPE_WIDTH_US, ParameterTypes.Double, "Shaped 90&#176; Width [&#956;s] (=shape resolution, 1pt/&#956;s) -- SWEEP to calibrate")
    TXAmplitude90 = Parameter("RFAsh0", _DEFAULT_RFASH0, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0&#8230;1 of LP max] -- default approx. from hard 90, SWEEP to calibrate", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", _LP_MAX_FRACTION, ParameterTypes.Double, "LP Port Max as Fraction of HP [0&#8230;1] -- your measured value")

    # Mains-lock trigger -- OFF by default.
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock Trigger Channel [1-3, unconfirmed for X-Pulse]")

    # Duty-cycle guard rail (RF only -- no gradients in this sequence)
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0&#8230;1]")

    # Phases (simple 2-step)
    PH1 = Parameter("PH1", "0,180", ParameterTypes.String, "Shaped 90&#176; Pulse Phase")
    PHRX = Parameter("PHRX", "0,180", ParameterTypes.String, "Acquisition Phase")

def run(comms):

    P = Parameters

    Phases = PhasesManager(P)
    Phases.Reset()

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

    ExcitationPulseWidth = P.P90sh + P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD

    ExBurpA = parse_burp_coeffs(P.ExBurpCoeffsA)
    ExBurpB = parse_burp_coeffs(P.ExBurpCoeffsB)

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

        # Switch to the LOW-POWER (LP) TX port for this shaped pulse.
        Transmit1SelectPort(1,0)
        Transmit1LPEnable(1,1)
        Delay(10000) # relay settle -- see LWG_Selective-Echo_H.py

        Receiver1Preamp(128, P.ReceiverAttenuation)
        Receiver1Filter(200, ReceiverFilter)

        Phases.Reset()

    comms.log("LP port engaged for shaped pulse -- LPMaxFraction={0:.2%} of "
              "HP full scale. Pulse offset: PulseOffset={1} Hz (applied on "
              "top of SF+O1={2:.6f} MHz)."
              .format(P.LPMaxFraction, P.PulseOffset, Frequency))

    # ---- Duty-cycle warning (RF only, no gradients here) -----------------
    rf_on_time = ExcitationPulseWidth
    estimate_rf_duty_cycle(P, rf_on_time, comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans,P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Receiver1FilterFlush(200, ReceiverFilter)
            # RD
            Delay(P.RecycleDelay-9.0e4)

            # Shaped 90 (excitation) pulse -- straight to acquisition, no
            # refocusing pulse and no TAU/echo delay.
            shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, P.TXAmplitude90, P.TXEnableTime, Frequency, P.PulseOffset, ExBurpA, ExBurpB)
            # ACQU
            Channel1SetBasePhase(P.TXEnableTime,0)
            Receiver1Phase(P.TXEnableTime, ph["PHRX"])
            # Dead1: probe ring-down time.
            Delay(P.Dead1)
            # ReceiverFilter.dead_time: digital-filter settling time.
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
# 1. Claude - 06/08/26 - Initial version. Stripped-down, excitation-only
#    "core" of LWG_Selective-Echo_H.py: same generate_shape()/shaped_pulse()
#    machinery (on-the-fly EBURP1/GAUSSIAN/SINC synthesis, LOW-power port,
#    dedicated PulseOffset Parameter, mains-lock-off default), minus the
#    refocusing pulse, TAU-TAU echo timing, and every Ref*/P180sh/RFAsh1
#    Parameter -- shaped-90 straight into acquisition. Built at your
#    request as an even more basic calibration tool than
#    LWG_Selective-Echo_H.py, for finding P90sh/RFAsh0/PulseOffset in
#    isolation before layering on refocusing-pulse calibration.
#    Also written with SHORT Parameter descriptions throughout (units/
#    critical info first, verbose rationale left to comments) since
#    SpinFlow's parameter panel has limited display width -- per your
#    request to make sure the useful parameters are actually readable
#    there, not just present.
# 2. Claude - 07/08/26 - POWER-CALC SIMPLIFIED OUT: removed the
#    PowerCalcMethod machinery entirely (manual/shape/db selector, plus the
#    RefAmplitude_HP/HardPulseWidth/ExcitationRotation/PowerAdjust_dB/
#    Excitation_dB/DbConvention Parameters and db_to_relative_scale()/
#    db_to_lp_relative_scale()/estimate_shape_relative_scale() functions),
#    per your feedback that it was over-engineered for what's fundamentally
#    a manual-sweep calibration tool. Checked Oxford's own WET-FID_H.py /
#    WET-FID-AUTO_H.py (Default pps folder) for a vendor formalism to use
#    instead -- they use a directly-calibrated Parameter with NO hard-
#    pulse-derived formula, so implemented your own suggestion instead:
#    RFAsh0's Parameter DEFAULT is now a single approximate value,
#    calculated ONCE at import time (not re-selectable at runtime) from a
#    calibrated hard 90 on the HP port and this shape's own integration
#    factor (see _HARD_P90_WIDTH_US/_HARD_P90_AMPLITUDE/_DEFAULT_RFASH0
#    above the Parameters block) -- kept shape_integration_factor() (the
#    only piece of the old machinery that's genuinely needed) since it's
#    what turns the shape's actual synthesised envelope into that estimate.
#    generate_shape() also now accepts 'EBURP2' (identical Fourier-series
#    code path to EBURP1) as prep for switching the default shape once real
#    E-BURP-2 coefficients/a Bruker shape file are supplied -- ExShape
#    still defaults to EBURP1 for now (see _DEFAULT_SHAPE_NAME).
# 3. Claude - 07/08/26 - FIXED shape_integration_factor(): was using
#    mean(|amp_profile|), which OVERESTIMATES the effective average B1 for
#    any shape with genuine negative lobes (BURP-family, SINC) by treating
#    every lobe as if it added constructively. Found by cross-checking
#    item 2's new RFAsh0 default against leonmr/bloch_sim.py's Bloch
#    simulation: the resulting pulse only excited ~18 degrees on-resonance
#    instead of the intended 90. Fixed to use the COHERENT (signed) mean
#    (amp*cos(radians(phase)), reconstructing the real-valued Fourier-
#    series b1(t)) -- correct because all sub-pulses share a fixed axis
#    on-resonance when phase is 0/180, so their rotations simply add
#    algebraically rather than by magnitude. For EBURP1 this raised
#    _DEFAULT_RFASH0 from 0.0231 to 0.1135 (~4.9x), confirmed to give a
#    clean 90 degree on-resonance excitation (Mz=0, |Mxy|=1) in
#    bloch_sim.py. GAUSSIAN (never negative) is unaffected by this fix;
#    SINC and any future BURP-family shape are.
#
# -----------------------------------------------------------------------------
