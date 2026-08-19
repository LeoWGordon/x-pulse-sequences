#-------------------------------------------------------------------------------
# Name:        LWG_Slice-Selective-Echo_H.py
# Purpose:     TRUE spatial SLICE-selective spin echo (shaped E-BURP-2
#              excitation + shaped RE-BURP refocusing, gradient ON during
#              BOTH pulses), on the Oxford Instruments X-Pulse Broadband
#              Benchtop NMR Spectrometer (1H/19F channel). Axis-selectable
#              (x/y/z). Acquires a localised echo from one slab/slice -- no
#              frequency-encoded readout/imaging in this file.
#
# Author:      Claude, for L. Gordon (DTU)
#
# Created:     06/08/2026
# Revised:     07/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.3
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
#    specifically flagged. E-BURP-2 (90) and RE-BURP (180) are DIFFERENT
#    pulse shapes with DIFFERENT bandwidth-time products (R = bandwidth x
#    duration; see R_EBURP2_90/R_REBURP_180 below). If you ran them both at
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
#    rephase gradient shown for either. The excitation pulse (E-BURP-2) IS
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
#  - EXCITED SLICE WIDTH: R_EBURP2_90 and R_REBURP_180 (the bandwidth-time
#    products) are now taken DIRECTLY from the real Bruker TopSpin shape
#    files' own ##$SHAPE_BWFAC metadata (EBurp2: 4.952000E00, ReBurp:
#    5.814000E00 -- see _EBURP2_TABLE_B1/_REBURP_TABLE_B1's own comments for
#    the source) -- no re-derivation needed, since these are the vendor's
#    own numbers for the exact shapes now driving hardware. (An earlier
#    version of this file, when the pulses were hand-reconstructed E-BURP-1/
#    RE-BURP, Bloch-simulated its own R values -- got 4.96/5.96,
#    reassuringly close to the current 4.952/5.814, but that was for
#    different pulses and has been superseded.) See
#    LWG_Slice-Selective-PulseAcquire_H.py's design notes for the original
#    Bloch-simulation method, still valid as an independent cross-check.
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

# ---- BURP bandwidth-time products -----------------------------------------
# R = (excitation/refocusing bandwidth [Hz]) x (pulse duration [s]), FWHM
# definition. Taken DIRECTLY from the real Bruker TopSpin EBurp2/ReBurp
# shape files' own ##$SHAPE_BWFAC metadata (see _EBURP2_TABLE_B1/
# _REBURP_TABLE_B1's own comments for the source and cross-validation) --
# no re-derivation needed, since these are the vendor's own numbers for the
# exact shapes now driving hardware. See LWG_Slice-Selective-
# PulseAcquire_H.py's design notes for the original Bloch-simulation method
# (Rodrigues'-formula rotation, on-resonance flip angle calibrated exactly
# via numerical search, then FWHM of the excitation/refocusing profile vs
# offset -- got 4.96/5.96, reassuringly close), still valid as an
# independent cross-check. Both are duration-independent properties of the
# pulse SHAPE alone (the Bloch equation is invariant under jointly
# rescaling time and all rates, so R does not need to be re-derived for
# whatever P90sh/P180sh you actually use).
R_EBURP2_90 = 4.952
R_REBURP_180 = 5.814

def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)
    points = P.ReceiverPoints
    DW = ReceiverFilter.dwell
    ReceiverTime = DW*(points+1)

    # Matches run()'s actual per-scan timing exactly (purely sequential --
    # no 'with parallel:' -- so this is every Delay()/duration argument in
    # the per-scan block summed directly, then algebraically simplified).
    # RD was previously OMITTED entirely -- by far the dominant term --
    # and the whole gradient/rephase/crusher event was collapsed into a
    # bare P90sh+P180sh+Tau*2 guess. The exact sum has striking
    # cancellations, all confirmed by symbolic expansion: RefocusPulseWidth
    # cancels out completely (RefocusLeadToCentre/RefocusCentreToEventEnd
    # each subtract half of it out of the TAUs, exactly matching the
    # refocusing event's own contribution); ExGradHoldRephase (and hence
    # RephaseFraction) cancels out completely too (subtracted out of
    # ExCentreToEventEnd exactly as it was added); and CrusherBlockTime
    # cancels out completely AND IS CRUSHERON-INDEPENDENT (each crusher's
    # own elapsed time is exactly compensated by being subtracted out of
    # the following TAU, whether or not a crusher actually runs). What
    # survives: RD + 2*Tau + half the excitation pulse width + one
    # PreGrad/RampTime pair (from ExCentreToEventEnd/RefocusLeadToCentre
    # not being perfectly symmetric) + Dead1 + 2.5*TXEnableTime (net,
    # after cancellation against the TAU-embedded half-TXEnableTime terms)
    # + half SHAPED_PULSE_FIXED_OVERHEAD + acquisition + a 205us fixed
    # small-instruction remainder.
    t_scanTime = (P.RecycleDelay + 2*P.Tau + P.P90sh/2.0 + P.PreGrad
                  + P.RampTime + P.Dead1 + ReceiverTime
                  + 2.5*P.TXEnableTime + SHAPED_PULSE_FIXED_OVERHEAD/2.0
                  + 205.0)
    t_acqTime = t_scanTime * (P.NumScans + P.DS)
    return t_acqTime

def sequence_description():

    seq_desc = ("Spatially SLICE-selective spin echo (shaped E-BURP-2 "
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

    basic = "NS,RD,NP,Filter,SlicePPM,TAU,P90sh,P180sh,AutoMatchSlice,RFAsh0,RFAsh1,G1,GradAxis,Probe,Dead1"

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

# E-BURP-2 (90 deg excitation) and RE-BURP (180 deg refocusing), as
# 1000-point SIGNED b1(tau) tables (peak-normalised to 1.0, tau uniformly
# spaced over [0,1)) -- exported DIRECTLY from real Bruker TopSpin shape
# files (JCAMP-DX 'Shape Data', ##$SHAPE_PARAMETERS= Type: EBurp2/ReBurp),
# supplied by L. Gordon 07/08/2026. Used instead of hand-entered Fourier-
# coefficient reconstructions (as EBURP1/BURP still are, via ExBurpCoeffsA/
# B or RefocusBurpCoeffsA/B) because these are numerically EXACT to
# Bruker's own reference shapes, with no risk of a transcribed-coefficient
# error driving real hardware. amp/phase are recovered from each table's
# sign (amp=|value|, phase=0 or 180) by _resample_burp_table() below -- see
# shape_integration_factor()'s docstring for why the SIGNED (not absolute-
# value) profile matters.
#
# CROSS-VALIDATED 07/08/2026: recomputing shape_integration_factor() on
# each table reproduces that shape file's own quoted SHAPE_INTEGFAC to 6
# significant figures (EBURP2: 6.102960E-02, REBURP: 7.981016E-02) --
# independent confirmation both that these tables were transcribed
# correctly AND that shape_integration_factor()'s coherent/signed-mean
# definition is the same one Bruker itself uses (see LWG_Selective-
# PulseAcquire_H.py's changelog for the full validation story).
_EBURP2_TABLE_B1 = (
    0.054881, 0.053111, 0.051408, 0.049769, 0.048195, 0.046685, 0.045236, 0.043849, 0.042521, 0.041253,
    0.040042, 0.038888, 0.037789, 0.036744, 0.035752, 0.034812, 0.033922, 0.033082, 0.032289, 0.031543,
    0.030842, 0.030184, 0.029569, 0.028995, 0.028461, 0.027965, 0.027505, 0.027081, 0.026690, 0.026332,
    0.026005, 0.025707, 0.025437, 0.025193, 0.024974, 0.024778, 0.024604, 0.024450, 0.024316, 0.024198,
    0.024097, 0.024010, 0.023935, 0.023873, 0.023820, 0.023776, 0.023739, 0.023709, 0.023682, 0.023659,
    0.023638, 0.023618, 0.023597, 0.023574, 0.023549, 0.023519, 0.023484, 0.023443, 0.023395, 0.023338,
    0.023272, 0.023196, 0.023109, 0.023010, 0.022899, 0.022774, 0.022635, 0.022482, 0.022313, 0.022128,
    0.021927, 0.021710, 0.021475, 0.021222, 0.020952, 0.020663, 0.020356, 0.020030, 0.019685, 0.019322,
    0.018939, 0.018537, 0.018117, 0.017677, 0.017219, 0.016741, 0.016245, 0.015731, 0.015199, 0.014648,
    0.014080, 0.013494, 0.012892, 0.012273, 0.011637, 0.010986, 0.010319, 0.009637, 0.008940, 0.008230,
    0.007505, 0.006768, 0.006018, 0.005256, 0.004482, 0.003697, 0.002901, 0.002095, 0.001280, 0.000456,
    -0.000377, -0.001218, -0.002066, -0.002921, -0.003783, -0.004651, -0.005525, -0.006403, -0.007286, -0.008174,
    -0.009065, -0.009960, -0.010858, -0.011758, -0.012661, -0.013565, -0.014471, -0.015379, -0.016287, -0.017197,
    -0.018106, -0.019016, -0.019926, -0.020836, -0.021745, -0.022654, -0.023562, -0.024469, -0.025375, -0.026280,
    -0.027184, -0.028086, -0.028987, -0.029886, -0.030784, -0.031680, -0.032574, -0.033466, -0.034357, -0.035245,
    -0.036132, -0.037017, -0.037900, -0.038781, -0.039660, -0.040537, -0.041412, -0.042285, -0.043155, -0.044024,
    -0.044890, -0.045754, -0.046615, -0.047474, -0.048331, -0.049185, -0.050037, -0.050885, -0.051731, -0.052574,
    -0.053413, -0.054249, -0.055082, -0.055911, -0.056736, -0.057558, -0.058375, -0.059187, -0.059995, -0.060797,
    -0.061595, -0.062387, -0.063173, -0.063953, -0.064727, -0.065493, -0.066253, -0.067005, -0.067750, -0.068486,
    -0.069214, -0.069933, -0.070643, -0.071343, -0.072034, -0.072713, -0.073382, -0.074040, -0.074686, -0.075320,
    -0.075942, -0.076551, -0.077147, -0.077729, -0.078297, -0.078851, -0.079390, -0.079914, -0.080423, -0.080916,
    -0.081392, -0.081852, -0.082295, -0.082721, -0.083130, -0.083521, -0.083894, -0.084248, -0.084584, -0.084901,
    -0.085199, -0.085478, -0.085737, -0.085977, -0.086196, -0.086396, -0.086576, -0.086735, -0.086874, -0.086993,
    -0.087090, -0.087168, -0.087224, -0.087260, -0.087274, -0.087268, -0.087241, -0.087193, -0.087124, -0.087034,
    -0.086923, -0.086792, -0.086639, -0.086466, -0.086271, -0.086056, -0.085820, -0.085564, -0.085286, -0.084988,
    -0.084670, -0.084330, -0.083971, -0.083590, -0.083189, -0.082768, -0.082326, -0.081864, -0.081381, -0.080878,
    -0.080354, -0.079810, -0.079245, -0.078659, -0.078053, -0.077426, -0.076778, -0.076109, -0.075419, -0.074707,
    -0.073974, -0.073220, -0.072444, -0.071646, -0.070826, -0.069983, -0.069118, -0.068230, -0.067319, -0.066385,
    -0.065427, -0.064445, -0.063439, -0.062408, -0.061353, -0.060272, -0.059166, -0.058034, -0.056876, -0.055692,
    -0.054480, -0.053242, -0.051976, -0.050682, -0.049359, -0.048008, -0.046628, -0.045219, -0.043781, -0.042312,
    -0.040813, -0.039284, -0.037724, -0.036132, -0.034510, -0.032856, -0.031170, -0.029452, -0.027702, -0.025920,
    -0.024105, -0.022258, -0.020378, -0.018466, -0.016521, -0.014543, -0.012532, -0.010489, -0.008414, -0.006306,
    -0.004166, -0.001995, 0.000209, 0.002444, 0.004710, 0.007007, 0.009334, 0.011692, 0.014079, 0.016495,
    0.018939, 0.021412, 0.023912, 0.026439, 0.028992, 0.031570, 0.034174, 0.036801, 0.039452, 0.042124,
    0.044819, 0.047534, 0.050269, 0.053023, 0.055795, 0.058584, 0.061389, 0.064209, 0.067043, 0.069889,
    0.072748, 0.075618, 0.078497, 0.081385, 0.084280, 0.087182, 0.090089, 0.093001, 0.095915, 0.098831,
    0.101749, 0.104666, 0.107582, 0.110495, 0.113405, 0.116310, 0.119209, 0.122102, 0.124986, 0.127862,
    0.130728, 0.133584, 0.136427, 0.139258, 0.142074, 0.144877, 0.147663, 0.150434, 0.153186, 0.155921,
    0.158637, 0.161334, 0.164009, 0.166664, 0.169297, 0.171907, 0.174494, 0.177058, 0.179597, 0.182110,
    0.184599, 0.187061, 0.189497, 0.191906, 0.194287, 0.196640, 0.198966, 0.201262, 0.203529, 0.205766,
    0.207974, 0.210151, 0.212298, 0.214414, 0.216499, 0.218552, 0.220574, 0.222563, 0.224521, 0.226445,
    0.228337, 0.230195, 0.232020, 0.233812, 0.235569, 0.237292, 0.238980, 0.240634, 0.242252, 0.243834,
    0.245381, 0.246891, 0.248365, 0.249801, 0.251200, 0.252561, 0.253884, 0.255167, 0.256412, 0.257616,
    0.258780, 0.259904, 0.260985, 0.262025, 0.263021, 0.263975, 0.264884, 0.265748, 0.266567, 0.267340,
    0.268065, 0.268743, 0.269372, 0.269952, 0.270481, 0.270959, 0.271386, 0.271759, 0.272078, 0.272343,
    0.272552, 0.272704, 0.272798, 0.272834, 0.272810, 0.272725, 0.272579, 0.272370, 0.272098, 0.271760,
    0.271358, 0.270888, 0.270351, 0.269745, 0.269070, 0.268324, 0.267507, 0.266617, 0.265654, 0.264616,
    0.263504, 0.262315, 0.261050, 0.259707, 0.258285, 0.256785, 0.255204, 0.253543, 0.251800, 0.249976,
    0.248069, 0.246079, 0.244005, 0.241847, 0.239605, 0.237277, 0.234864, 0.232365, 0.229780, 0.227109,
    0.224351, 0.221506, 0.218574, 0.215555, 0.212449, 0.209255, 0.205974, 0.202606, 0.199150, 0.195608,
    0.191978, 0.188261, 0.184458, 0.180568, 0.176592, 0.172530, 0.168382, 0.164149, 0.159830, 0.155427,
    0.150940, 0.146368, 0.141713, 0.136975, 0.132154, 0.127251, 0.122267, 0.117201, 0.112054, 0.106827,
    0.101520, 0.096134, 0.090670, 0.085127, 0.079506, 0.073808, 0.068034, 0.062183, 0.056256, 0.050255,
    0.044178, 0.038028, 0.031804, 0.025506, 0.019136, 0.012694, 0.006180, -0.000406, -0.007062, -0.013789,
    -0.020586, -0.027452, -0.034387, -0.041390, -0.048462, -0.055602, -0.062808, -0.070081, -0.077421, -0.084826,
    -0.092296, -0.099831, -0.107430, -0.115093, -0.122819, -0.130607, -0.138457, -0.146368, -0.154340, -0.162371,
    -0.170461, -0.178608, -0.186814, -0.195075, -0.203392, -0.211763, -0.220188, -0.228664, -0.237191, -0.245768,
    -0.254393, -0.263065, -0.271782, -0.280544, -0.289347, -0.298190, -0.307072, -0.315991, -0.324944, -0.333930,
    -0.342946, -0.351990, -0.361059, -0.370151, -0.379264, -0.388394, -0.397540, -0.406697, -0.415863, -0.425034,
    -0.434208, -0.443381, -0.452550, -0.461711, -0.470860, -0.479993, -0.489108, -0.498198, -0.507262, -0.516294,
    -0.525289, -0.534245, -0.543156, -0.552017, -0.560825, -0.569574, -0.578259, -0.586876, -0.595421, -0.603887,
    -0.612270, -0.620564, -0.628766, -0.636869, -0.644868, -0.652759, -0.660535, -0.668192, -0.675725, -0.683128,
    -0.690396, -0.697524, -0.704506, -0.711338, -0.718015, -0.724530, -0.730880, -0.737059, -0.743063, -0.748886,
    -0.754523, -0.759971, -0.765225, -0.770279, -0.775130, -0.779774, -0.784206, -0.788422, -0.792419, -0.796192,
    -0.799739, -0.803055, -0.806138, -0.808984, -0.811591, -0.813955, -0.816074, -0.817945, -0.819567, -0.820938,
    -0.822054, -0.822916, -0.823520, -0.823867, -0.823955, -0.823783, -0.823350, -0.822656, -0.821700, -0.820484,
    -0.819006, -0.817267, -0.815267, -0.813008, -0.810490, -0.807715, -0.804683, -0.801396, -0.797856, -0.794065,
    -0.790025, -0.785737, -0.781205, -0.776431, -0.771418, -0.766167, -0.760684, -0.754969, -0.749028, -0.742863,
    -0.736477, -0.729875, -0.723060, -0.716035, -0.708806, -0.701375, -0.693748, -0.685927, -0.677917, -0.669723,
    -0.661348, -0.652798, -0.644076, -0.635187, -0.626135, -0.616924, -0.607560, -0.598045, -0.588386, -0.578585,
    -0.568648, -0.558578, -0.548381, -0.538059, -0.527617, -0.517059, -0.506390, -0.495612, -0.484730, -0.473748,
    -0.462668, -0.451496, -0.440234, -0.428885, -0.417452, -0.405940, -0.394351, -0.382687, -0.370952, -0.359149,
    -0.347279, -0.335346, -0.323351, -0.311298, -0.299187, -0.287022, -0.274804, -0.262534, -0.250215, -0.237848,
    -0.225435, -0.212976, -0.200474, -0.187930, -0.175344, -0.162717, -0.150052, -0.137348, -0.124606, -0.111828,
    -0.099013, -0.086164, -0.073279, -0.060360, -0.047408, -0.034423, -0.021406, -0.008357, 0.004723, 0.017834,
    0.030975, 0.044145, 0.057345, 0.070572, 0.083826, 0.097107, 0.110414, 0.123746, 0.137101, 0.150479,
    0.163879, 0.177299, 0.190738, 0.204195, 0.217667, 0.231154, 0.244654, 0.258164, 0.271682, 0.285207,
    0.298736, 0.312267, 0.325797, 0.339323, 0.352842, 0.366352, 0.379849, 0.393330, 0.406792, 0.420230,
    0.433641, 0.447022, 0.460367, 0.473674, 0.486936, 0.500151, 0.513313, 0.526418, 0.539460, 0.552436,
    0.565338, 0.578163, 0.590905, 0.603558, 0.616118, 0.628577, 0.640931, 0.653174, 0.665299, 0.677301,
    0.689174, 0.700912, 0.712509, 0.723958, 0.735253, 0.746389, 0.757359, 0.768157, 0.778777, 0.789213,
    0.799458, 0.809508, 0.819355, 0.828994, 0.838420, 0.847625, 0.856606, 0.865356, 0.873870, 0.882143,
    0.890169, 0.897944, 0.905462, 0.912720, 0.919712, 0.926435, 0.932883, 0.939054, 0.944943, 0.950547,
    0.955862, 0.960886, 0.965615, 0.970048, 0.974180, 0.978011, 0.981538, 0.984760, 0.987675, 0.990283,
    0.992581, 0.994570, 0.996249, 0.997618, 0.998677, 0.999426, 0.999867, 1.000000, 0.999826, 0.999347,
    0.998564, 0.997480, 0.996096, 0.994416, 0.992441, 0.990175, 0.987621, 0.984783, 0.981663, 0.978267,
    0.974598, 0.970660, 0.966458, 0.961996, 0.957280, 0.952315, 0.947106, 0.941658, 0.935976, 0.930067,
    0.923937, 0.917592, 0.911037, 0.904279, 0.897324, 0.890178, 0.882850, 0.875344, 0.867667, 0.859827,
    0.851830, 0.843682, 0.835392, 0.826966, 0.818410, 0.809731, 0.800937, 0.792034, 0.783030, 0.773930,
    0.764743, 0.755474, 0.746130, 0.736719, 0.727246, 0.717718, 0.708141, 0.698523, 0.688868, 0.679184,
    0.669476, 0.659750, 0.650012, 0.640268, 0.630523, 0.620783, 0.611053, 0.601337, 0.591642, 0.581972,
    0.572331, 0.562725, 0.553158, 0.543633, 0.534155, 0.524729, 0.515357, 0.506044, 0.496793, 0.487607,
    0.478490, 0.469444, 0.460473, 0.451580, 0.442766, 0.434034, 0.425387, 0.416827, 0.408355, 0.399974,
    0.391686, 0.383491, 0.375393, 0.367391, 0.359487, 0.351682, 0.343978, 0.336374, 0.328873, 0.321475,
    0.314179, 0.306988, 0.299900, 0.292917, 0.286039, 0.279265, 0.272597, 0.266033, 0.259574, 0.253220,
    0.246971, 0.240826, 0.234785, 0.228847, 0.223013, 0.217282, 0.211653, 0.206126, 0.200700, 0.195376,
    0.190151, 0.185026, 0.180000, 0.175072, 0.170242, 0.165508, 0.160871, 0.156329, 0.151882, 0.147529,
    0.143269, 0.139101, 0.135025, 0.131040, 0.127145, 0.123340, 0.119623, 0.115993, 0.112451, 0.108995,
    0.105624, 0.102338, 0.099136, 0.096017, 0.092979, 0.090024, 0.087148, 0.084353, 0.081636, 0.078998,
    0.076437, 0.073952, 0.071543, 0.069208, 0.066947, 0.064759, 0.062643, 0.060598, 0.058624, 0.056718,
)

_REBURP_TABLE_B1 = (
    -0.074916, -0.074916, -0.074855, -0.074733, -0.074552, -0.074310, -0.074011, -0.073655, -0.073244, -0.072780,
    -0.072266, -0.071702, -0.071093, -0.070442, -0.069749, -0.069020, -0.068258, -0.067465, -0.066645, -0.065801,
    -0.064937, -0.064057, -0.063164, -0.062260, -0.061351, -0.060438, -0.059525, -0.058615, -0.057711, -0.056815,
    -0.055930, -0.055058, -0.054201, -0.053361, -0.052539, -0.051737, -0.050955, -0.050195, -0.049456, -0.048740,
    -0.048045, -0.047372, -0.046720, -0.046087, -0.045474, -0.044878, -0.044298, -0.043732, -0.043178, -0.042635,
    -0.042099, -0.041568, -0.041040, -0.040512, -0.039982, -0.039446, -0.038902, -0.038346, -0.037777, -0.037191,
    -0.036586, -0.035959, -0.035307, -0.034629, -0.033921, -0.033182, -0.032410, -0.031602, -0.030759, -0.029877,
    -0.028956, -0.027995, -0.026993, -0.025950, -0.024866, -0.023741, -0.022574, -0.021367, -0.020120, -0.018833,
    -0.017509, -0.016149, -0.014754, -0.013325, -0.011865, -0.010376, -0.008859, -0.007318, -0.005755, -0.004171,
    -0.002570, -0.000954, 0.000673, 0.002311, 0.003955, 0.005604, 0.007254, 0.008903, 0.010549, 0.012190,
    0.013822, 0.015445, 0.017056, 0.018652, 0.020234, 0.021798, 0.023343, 0.024869, 0.026374, 0.027858,
    0.029319, 0.030757, 0.032171, 0.033562, 0.034930, 0.036274, 0.037595, 0.038894, 0.040170, 0.041424,
    0.042658, 0.043872, 0.045067, 0.046244, 0.047404, 0.048548, 0.049678, 0.050794, 0.051898, 0.052990,
    0.054073, 0.055146, 0.056212, 0.057270, 0.058322, 0.059369, 0.060411, 0.061449, 0.062483, 0.063514,
    0.064542, 0.065567, 0.066589, 0.067608, 0.068623, 0.069634, 0.070641, 0.071643, 0.072638, 0.073627,
    0.074608, 0.075579, 0.076541, 0.077490, 0.078426, 0.079348, 0.080254, 0.081141, 0.082010, 0.082857,
    0.083681, 0.084480, 0.085254, 0.085999, 0.086715, 0.087400, 0.088052, 0.088669, 0.089251, 0.089795,
    0.090301, 0.090767, 0.091193, 0.091576, 0.091917, 0.092214, 0.092467, 0.092674, 0.092836, 0.092952,
    0.093022, 0.093046, 0.093023, 0.092953, 0.092838, 0.092676, 0.092468, 0.092214, 0.091915, 0.091572,
    0.091184, 0.090753, 0.090279, 0.089763, 0.089205, 0.088606, 0.087968, 0.087289, 0.086573, 0.085818,
    0.085026, 0.084198, 0.083333, 0.082434, 0.081499, 0.080531, 0.079529, 0.078494, 0.077426, 0.076326,
    0.075193, 0.074029, 0.072834, 0.071607, 0.070348, 0.069059, 0.067739, 0.066387, 0.065005, 0.063591,
    0.062146, 0.060670, 0.059163, 0.057624, 0.056054, 0.054452, 0.052818, 0.051153, 0.049456, 0.047727,
    0.045966, 0.044174, 0.042349, 0.040493, 0.038606, 0.036687, 0.034737, 0.032757, 0.030746, 0.028705,
    0.026634, 0.024535, 0.022406, 0.020250, 0.018067, 0.015856, 0.013620, 0.011359, 0.009073, 0.006764,
    0.004432, 0.002078, -0.000297, -0.002691, -0.005105, -0.007537, -0.009986, -0.012451, -0.014931, -0.017426,
    -0.019933, -0.022454, -0.024986, -0.027529, -0.030081, -0.032643, -0.035212, -0.037790, -0.040373, -0.042963,
    -0.045558, -0.048158, -0.050761, -0.053367, -0.055977, -0.058588, -0.061200, -0.063814, -0.066427, -0.069041,
    -0.071653, -0.074264, -0.076872, -0.079478, -0.082080, -0.084678, -0.087272, -0.089860, -0.092441, -0.095015,
    -0.097581, -0.100138, -0.102685, -0.105221, -0.107744, -0.110255, -0.112750, -0.115230, -0.117692, -0.120136,
    -0.122560, -0.124961, -0.127340, -0.129693, -0.132019, -0.134317, -0.136585, -0.138820, -0.141022, -0.143187,
    -0.145314, -0.147402, -0.149448, -0.151449, -0.153405, -0.155313, -0.157171, -0.158977, -0.160730, -0.162427,
    -0.164067, -0.165647, -0.167167, -0.168624, -0.170016, -0.171343, -0.172602, -0.173792, -0.174912, -0.175961,
    -0.176936, -0.177838, -0.178664, -0.179414, -0.180087, -0.180682, -0.181197, -0.181633, -0.181989, -0.182263,
    -0.182456, -0.182566, -0.182593, -0.182536, -0.182395, -0.182169, -0.181858, -0.181461, -0.180978, -0.180407,
    -0.179750, -0.179004, -0.178169, -0.177245, -0.176231, -0.175125, -0.173928, -0.172638, -0.171255, -0.169777,
    -0.168203, -0.166532, -0.164763, -0.162895, -0.160926, -0.158855, -0.156681, -0.154401, -0.152015, -0.149521,
    -0.146916, -0.144201, -0.141372, -0.138428, -0.135367, -0.132188, -0.128889, -0.125468, -0.121923, -0.118253,
    -0.114456, -0.110530, -0.106474, -0.102286, -0.097964, -0.093508, -0.088916, -0.084187, -0.079319, -0.074312,
    -0.069164, -0.063874, -0.058442, -0.052867, -0.047148, -0.041285, -0.035276, -0.029122, -0.022823, -0.016377,
    -0.009786, -0.003048, 0.003836, 0.010867, 0.018044, 0.025367, 0.032837, 0.040452, 0.048213, 0.056120,
    0.064172, 0.072370, 0.080713, 0.089201, 0.097833, 0.106611, 0.115533, 0.124600, 0.133811, 0.143166,
    0.152666, 0.162311, 0.172100, 0.182033, 0.192111, 0.202333, 0.212700, 0.223211, 0.233866, 0.244665,
    0.255608, 0.266694, 0.277924, 0.289295, 0.300808, 0.312462, 0.324255, 0.336185, 0.348252, 0.360453,
    0.372786, 0.385249, 0.397837, 0.410549, 0.423379, 0.436325, 0.449380, 0.462541, 0.475801, 0.489154,
    0.502594, 0.516112, 0.529702, 0.543354, 0.557061, 0.570811, 0.584595, 0.598402, 0.612221, 0.626040,
    0.639845, 0.653625, 0.667364, 0.681050, 0.694667, 0.708200, 0.721633, 0.734951, 0.748136, 0.761172,
    0.774043, 0.786730, 0.799216, 0.811483, 0.823513, 0.835289, 0.846793, 0.858006, 0.868910, 0.879489,
    0.889725, 0.899600, 0.909099, 0.918204, 0.926900, 0.935171, 0.943002, 0.950380, 0.957291, 0.963722,
    0.969661, 0.975097, 0.980019, 0.984419, 0.988288, 0.991618, 0.994403, 0.996637, 0.998317, 0.999439,
    1.000000, 1.000000, 0.999439, 0.998317, 0.996637, 0.994403, 0.991618, 0.988288, 0.984419, 0.980019,
    0.975097, 0.969661, 0.963722, 0.957291, 0.950380, 0.943002, 0.935171, 0.926900, 0.918204, 0.909099,
    0.899600, 0.889725, 0.879489, 0.868910, 0.858006, 0.846793, 0.835289, 0.823513, 0.811483, 0.799216,
    0.786730, 0.774043, 0.761172, 0.748136, 0.734951, 0.721633, 0.708200, 0.694667, 0.681050, 0.667364,
    0.653625, 0.639845, 0.626040, 0.612221, 0.598402, 0.584595, 0.570811, 0.557061, 0.543354, 0.529702,
    0.516112, 0.502594, 0.489154, 0.475801, 0.462541, 0.449380, 0.436325, 0.423379, 0.410549, 0.397837,
    0.385249, 0.372786, 0.360453, 0.348252, 0.336185, 0.324255, 0.312462, 0.300808, 0.289295, 0.277924,
    0.266694, 0.255608, 0.244665, 0.233866, 0.223211, 0.212700, 0.202333, 0.192111, 0.182033, 0.172100,
    0.162311, 0.152666, 0.143166, 0.133811, 0.124600, 0.115533, 0.106611, 0.097833, 0.089201, 0.080713,
    0.072370, 0.064172, 0.056120, 0.048213, 0.040452, 0.032837, 0.025367, 0.018044, 0.010867, 0.003836,
    -0.003048, -0.009786, -0.016377, -0.022823, -0.029122, -0.035276, -0.041285, -0.047148, -0.052867, -0.058442,
    -0.063874, -0.069164, -0.074312, -0.079319, -0.084187, -0.088916, -0.093508, -0.097964, -0.102286, -0.106474,
    -0.110530, -0.114456, -0.118253, -0.121923, -0.125468, -0.128889, -0.132188, -0.135367, -0.138428, -0.141372,
    -0.144201, -0.146916, -0.149521, -0.152015, -0.154401, -0.156681, -0.158855, -0.160926, -0.162895, -0.164763,
    -0.166532, -0.168203, -0.169777, -0.171255, -0.172638, -0.173928, -0.175125, -0.176231, -0.177245, -0.178169,
    -0.179004, -0.179750, -0.180407, -0.180978, -0.181461, -0.181858, -0.182169, -0.182395, -0.182536, -0.182593,
    -0.182566, -0.182456, -0.182263, -0.181989, -0.181633, -0.181197, -0.180682, -0.180087, -0.179414, -0.178664,
    -0.177838, -0.176936, -0.175961, -0.174912, -0.173792, -0.172602, -0.171343, -0.170016, -0.168624, -0.167167,
    -0.165647, -0.164067, -0.162427, -0.160730, -0.158977, -0.157171, -0.155313, -0.153405, -0.151449, -0.149448,
    -0.147402, -0.145314, -0.143187, -0.141022, -0.138820, -0.136585, -0.134317, -0.132019, -0.129693, -0.127340,
    -0.124961, -0.122560, -0.120136, -0.117692, -0.115230, -0.112750, -0.110255, -0.107744, -0.105221, -0.102685,
    -0.100138, -0.097581, -0.095015, -0.092441, -0.089860, -0.087272, -0.084678, -0.082080, -0.079478, -0.076872,
    -0.074264, -0.071653, -0.069041, -0.066427, -0.063814, -0.061200, -0.058588, -0.055977, -0.053367, -0.050761,
    -0.048158, -0.045558, -0.042963, -0.040373, -0.037790, -0.035212, -0.032643, -0.030081, -0.027529, -0.024986,
    -0.022454, -0.019933, -0.017426, -0.014931, -0.012451, -0.009986, -0.007537, -0.005105, -0.002691, -0.000297,
    0.002078, 0.004432, 0.006764, 0.009073, 0.011359, 0.013620, 0.015856, 0.018067, 0.020250, 0.022406,
    0.024535, 0.026634, 0.028705, 0.030746, 0.032757, 0.034737, 0.036687, 0.038606, 0.040493, 0.042349,
    0.044174, 0.045966, 0.047727, 0.049456, 0.051153, 0.052818, 0.054452, 0.056054, 0.057624, 0.059163,
    0.060670, 0.062146, 0.063591, 0.065005, 0.066387, 0.067739, 0.069059, 0.070348, 0.071607, 0.072834,
    0.074029, 0.075193, 0.076326, 0.077426, 0.078494, 0.079529, 0.080531, 0.081499, 0.082434, 0.083333,
    0.084198, 0.085026, 0.085818, 0.086573, 0.087289, 0.087968, 0.088606, 0.089205, 0.089763, 0.090279,
    0.090753, 0.091184, 0.091572, 0.091915, 0.092214, 0.092468, 0.092676, 0.092838, 0.092953, 0.093023,
    0.093046, 0.093022, 0.092952, 0.092836, 0.092674, 0.092467, 0.092214, 0.091917, 0.091576, 0.091193,
    0.090767, 0.090301, 0.089795, 0.089251, 0.088669, 0.088052, 0.087400, 0.086715, 0.085999, 0.085254,
    0.084480, 0.083681, 0.082857, 0.082010, 0.081141, 0.080254, 0.079348, 0.078426, 0.077490, 0.076541,
    0.075579, 0.074608, 0.073627, 0.072638, 0.071643, 0.070641, 0.069634, 0.068623, 0.067608, 0.066589,
    0.065567, 0.064542, 0.063514, 0.062483, 0.061449, 0.060411, 0.059369, 0.058322, 0.057270, 0.056212,
    0.055146, 0.054073, 0.052990, 0.051898, 0.050794, 0.049678, 0.048548, 0.047404, 0.046244, 0.045067,
    0.043872, 0.042658, 0.041424, 0.040170, 0.038894, 0.037595, 0.036274, 0.034930, 0.033562, 0.032171,
    0.030757, 0.029319, 0.027858, 0.026374, 0.024869, 0.023343, 0.021798, 0.020234, 0.018652, 0.017056,
    0.015445, 0.013822, 0.012190, 0.010549, 0.008903, 0.007254, 0.005604, 0.003955, 0.002311, 0.000673,
    -0.000954, -0.002570, -0.004171, -0.005755, -0.007318, -0.008859, -0.010376, -0.011865, -0.013325, -0.014754,
    -0.016149, -0.017509, -0.018833, -0.020120, -0.021367, -0.022574, -0.023741, -0.024866, -0.025950, -0.026993,
    -0.027995, -0.028956, -0.029877, -0.030759, -0.031602, -0.032410, -0.033182, -0.033921, -0.034629, -0.035307,
    -0.035959, -0.036586, -0.037191, -0.037777, -0.038346, -0.038902, -0.039446, -0.039982, -0.040512, -0.041040,
    -0.041568, -0.042099, -0.042635, -0.043178, -0.043732, -0.044298, -0.044878, -0.045474, -0.046087, -0.046720,
    -0.047372, -0.048045, -0.048740, -0.049456, -0.050195, -0.050955, -0.051737, -0.052539, -0.053361, -0.054201,
    -0.055058, -0.055930, -0.056815, -0.057711, -0.058615, -0.059525, -0.060438, -0.061351, -0.062260, -0.063164,
    -0.064057, -0.064937, -0.065801, -0.066645, -0.067465, -0.068258, -0.069020, -0.069749, -0.070442, -0.071093,
    -0.071702, -0.072266, -0.072780, -0.073244, -0.073655, -0.074011, -0.074310, -0.074552, -0.074733, -0.074855,
)


def _resample_burp_table(table, n_steps):
    """Resample a fixed reference SIGNED b1(tau) table (e.g.
    _EBURP2_TABLE_B1/_REBURP_TABLE_B1, uniformly spaced over [0,1)) onto
    n_steps points via linear interpolation of the SIGNED profile -- NOT of
    separately-folded amp/phase, which would create spurious intermediate
    phase values (e.g. ~90 deg) right at a 0/180 sign-flip boundary.
    Returns (amp, phase_deg) in the same convention as the rest of
    generate_shape()."""
    table = np.asarray(table, dtype=float)
    tau_table = np.linspace(0.0, 1.0, len(table), endpoint=False)
    tau_out = np.linspace(0.0, 1.0, n_steps, endpoint=False)
    signed = np.interp(tau_out, tau_table, table)
    return np.abs(signed), np.where(signed < 0, 180.0, 0.0)


def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical machinery to LWG_Selective-Echo_H.py /
    LWG_Slice-Selective-PulseAcquire_H.py; duplicated here so this file
    stays fully self-contained.

    shape_name: 'GAUSSIAN', 'SINC', 'EBURP2'/'REBURP' (fixed Bruker-sourced
    tables, no coefficients needed), or 'EBURP1'/'BURP' (needs
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

    elif shape_name == 'EBURP2':
        amp, phase = _resample_burp_table(_EBURP2_TABLE_B1, n_steps)

    elif shape_name == 'REBURP':
        amp, phase = _resample_burp_table(_REBURP_TABLE_B1, n_steps)

    elif shape_name in ('EBURP1', 'BURP'):
        if not burp_coeffs_A:
            raise ValueError(
                "shape '{0}' requested but no Fourier A-coefficients were "
                "supplied. EBURP1 is defined in Geen, H. &amp; "
                "Freeman, R., J. Magn. Reson. 93, 93-141 (1991) as a "
                "truncated Fourier series -- enter coefficients as "
                "ExBurpCoeffsA/B or RefocusBurpCoeffsA/B, or use 'GAUSSIAN' "
                "/ 'SINC' / 'EBURP2' / 'REBURP' instead, which need no "
                "coefficients."
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
                          "SINC, EBURP2, REBURP, EBURP1/BURP (needs "
                          "BurpCoeffsA/B).".format(shape_name))

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
    lobes (EBURP1: ~4.9x too high -- see LWG_Selective-PulseAcquire_H.py's
    changelog for the full validation story, confirmed by Bloch-simulating
    the resulting 'default' RFAsh0 in leonmr/bloch_sim.py).

    Used once, at import time below, to turn a calibrated hard-90 into
    approximate starting RFAsh0/RFAsh1 defaults. VERIFY any shape/power
    combination with leonmr/bloch_sim.py (or a real nutation curve) before
    trusting it quantitatively -- this is still only an approximation.
    """
    peak = np.max(np.abs(amp_profile))
    if peak <= 0:
        raise ValueError("shape_integration_factor: shape has zero peak amplitude")
    signed_profile = np.asarray(amp_profile) * np.cos(np.deg2rad(phase_profile_deg))
    return float(np.mean(signed_profile) / peak)


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


# Per-axis max gradient strength (G/cm) for probes that can be fitted to
# this magnet. There is no way to detect which probe is mounted from
# software, so this is a MANUAL selector (Probe Parameter below) -- the
# operator must set it to match what is actually on the magnet. KEEP THIS
# TABLE IN SYNC with leonmr/xpulse_imaging.py's own GRADIENT_CALIBRATION
# dict (duplicated here rather than imported, since pulse programmes can't
# import external Python modules).
MAXGRAD_TABLE = {
    "HFX": {"x": 11.879, "y": 11.978, "z": 57.915},   # G/cm, measured/averaged 06/08/2026
    "LOWGAMMA": {"x": None, "y": None, "z": None},     # NOT YET CALIBRATED
}

def report_probe_gradient(P, comms, axis=None):
    """Log the max gradient strength (G/cm) for the currently-selected
    Probe and slice-select gradient axis (GradAxis). Since Probe is a
    Parameter, this also gets auto-recorded in the resulting JCAMP file's
    SpinFlow parameter block. Warns (does not fail) if the fitted probe's
    selected axis has no calibration yet."""
    probe_key = str(P.Probe).upper().replace('/', '').replace('-', '').replace(' ', '').replace('_', '')
    axis_key = str(axis if axis is not None else getattr(P, "Axis", "z")).strip().lower()
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


# ---- Approximate default RFAsh0/RFAsh1, calculated ONCE (at import time,
# not per run) from a calibrated hard 90 pulse on the HP port -- replaces
# the PowerCalcMethod machinery this file used to have (manual/shape/db
# selector, RefAmplitude_HP/HardPulseWidth/ExRotation/RefocusRotation/
# PowerAdjust_dB/Excitation_dB/Refocus_dB/DbConvention Parameters),
# simplified out 07/08/26. ALWAYS verify with a real nutation curve before
# trusting this quantitatively; it doesn't know your actual coil/sample.
#
#   RFAsh_x (LP-relative) ~= HardAmplitude * HardWidth * TargetRotation
#                            / (ShapeWidth * 90 * IntegrationFactor) / LPMaxFraction
#
_HARD_P90_WIDTH_US = 9.58          # this repo's usual calibrated hard-90 width (P90/P1Hard)
_HARD_P90_AMPLITUDE = 0.4          # ...at this HP-relative amplitude (RFA0)
_LP_MAX_FRACTION = 0.10            # measured LP-vs-HP max-power ratio (LPMaxFraction default, below)
_DEFAULT_EX_SHAPE_NAME = 'EBURP2'      # exact Bruker-sourced table -- see _EBURP2_TABLE_B1 above
_DEFAULT_REFOCUS_SHAPE_NAME = 'REBURP' # exact Bruker-sourced table -- see _REBURP_TABLE_B1 above
_DEFAULT_EX_SHAPE_WIDTH_US = 5000.0    # matches P90sh's own default, below
# P180sh's own default is DERIVED from P90sh's via the R_REBURP_180/
# R_EBURP2_90 ratio (see 'MATCHING THE EXCITATION AND REFOCUSING SLICES' in
# the file header), so the two slices are consistent even before
# AutoMatchSlice recomputes it at run time.
_DEFAULT_REFOCUS_SHAPE_WIDTH_US = round(_DEFAULT_EX_SHAPE_WIDTH_US * (R_REBURP_180 / R_EBURP2_90), 1)
# Only used as ExBurpCoeffsA/B/RefocusBurpCoeffsA/B's own Parameter
# defaults (i.e. if you switch ExShape/RefocusShape to 'EBURP1'/'BURP') --
# NOT used for EBURP2/REBURP, which are fixed tables, not coefficient
# reconstructions.
_DEFAULT_EXBURP_A = "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00"
_DEFAULT_EXBURP_B = "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01"
_DEFAULT_REFOCUSBURP_A = "0.49,-1.02,1.11,-1.57,0.83,-0.42,0.26,-0.16,0.10,-0.07,0.04,-0.03,0.01,-0.02,0.00,-0.01"
_DEFAULT_REFOCUSBURP_B = ""

# Precomputed 07/08/2026 from the formula above (generate_shape('EBURP2',
# 5000, ...)/('REBURP', 5870.4, ...) -> shape_integration_factor(...) =
# 0.0610296/~0.0680 -- REBURP's width here is the AUTO-MATCHED
# _DEFAULT_REFOCUS_SHAPE_WIDTH_US above, not a fixed 10000us, so RFAsh1
# differs from the plain Selective-Echo files' 0.096). Deliberately NOT
# evaluated at import time: SpinFlow executes this entire file just to
# LOAD the sequence into the parameter panel (before run() is ever
# called), so any on-the-fly shape synthesis at module scope runs on
# EVERY load, not just every scan -- a needless dependency for values
# that are meant to be swept by hand anyway. generate_shape()/
# shape_integration_factor() remain fully available and are used normally
# inside run(); only these one-time calibration defaults are now
# literals.
_DEFAULT_RFASH0 = 0.1256
_DEFAULT_RFASH1 = 0.1636


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
    # Defaults = exact Bruker-sourced tables: EXCITATION = E-BURP-2,
    # REFOCUSING = RE-BURP (see _EBURP2_TABLE_B1/_REBURP_TABLE_B1 above).
    ExcitationShape = Parameter("ExShape", _DEFAULT_EX_SHAPE_NAME, ParameterTypes.String, "Excitation (90) Shape [EBURP2(default, exact Bruker table)/EBURP1/GAUSSIAN/SINC/BURP] -- see ExBurpCoeffsA/B")
    RefocusShape = Parameter("RefocusShape", _DEFAULT_REFOCUS_SHAPE_NAME, ParameterTypes.String, "Refocusing (180) Shape [REBURP(default, exact Bruker table)/EBURP1/GAUSSIAN/SINC/BURP] -- see RefocusBurpCoeffsA/B")
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", _DEFAULT_EXBURP_A, ParameterTypes.String, "E-BURP-1 cosine coeffs A0..A8 (Geen &amp; Freeman'91 Tbl.2) -- for ExShape=EBURP1/BURP only")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", _DEFAULT_EXBURP_B, ParameterTypes.String, "E-BURP-1 sine coeffs B0..B8 (Geen &amp; Freeman'91 Tbl.2) -- for ExShape=EBURP1/BURP only")
    RefocusBurpCoeffsA = Parameter("RefocusBurpCoeffsA", _DEFAULT_REFOCUSBURP_A, ParameterTypes.String, "RE-BURP cosine coeffs A0..A15 (Geen &amp; Freeman'91 Tbl.8) -- for RefocusShape=EBURP1/BURP only")
    RefocusBurpCoeffsB = Parameter("RefocusBurpCoeffsB", _DEFAULT_REFOCUSBURP_B, ParameterTypes.String, "RE-BURP sine coeffs -- none published, leave empty")

    P90sh = Parameter("P90sh", _DEFAULT_EX_SHAPE_WIDTH_US, ParameterTypes.Double, "Shaped 90&#176; (excitation) Width [&#956;s] -- also sets shape res. + excitation slice width")
    # P180sh's displayed default already matches P90sh's default via the
    # R_REBURP_180/R_EBURP2_90 ratio (5.814/4.952), so the two slices are
    # consistent even before AutoMatchSlice recomputes it at run time.
    P180sh = Parameter("P180sh", _DEFAULT_REFOCUS_SHAPE_WIDTH_US, ParameterTypes.Double, "Shaped 180&#176; (refocusing) Width [&#956;s] -- IGNORED/recomputed from P90sh if AutoMatchSlice=1")
    AutoMatchSlice = Parameter("AutoMatchSlice", 1, ParameterTypes.Int32, "*** Auto-derive P180sh from P90sh to match slices [1=On(default,RECOMMENDED), 0=manual P180sh] ***")

    TXAmplitude90 = Parameter("RFAsh0", _DEFAULT_RFASH0, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0&#8230;1.0 of LP max] -- approx. from hard-90 calib., CALIBRATE via nutation before trusting", RFA, min=0.0, max=1.0)
    TXAmplitude180 = Parameter("RFAsh1", _DEFAULT_RFASH1, ParameterTypes.Double,
                               "Shaped 180&#176; TX Power [0&#8230;1.0 of LP max] -- approx. from hard-90 calib., CALIBRATE via nutation before trusting", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", _LP_MAX_FRACTION, ParameterTypes.Double, "LP Port Max as Fraction of HP [0&#8230;1.0] -- your measured value")

    # Slice-select gradient -- SAME gradient (G1, same axis) is used for
    # BOTH pulses, by design (see 'MATCHING THE EXCITATION AND REFOCUSING
    # SLICES' above) -- there is deliberately no separate G-for-180 knob.
    G1 = Parameter("G1", 1.0, ParameterTypes.Double, "Slice-Select Gradient Strength [-1.0&#8230;1.0] -- used for BOTH pulses")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Slice-Select Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z", ParameterTypes.String, "Gradient Axes")
    # Probe fitted to the magnet -- no automatic detection is possible, so
    # this must be set MANUALLY to match what's actually mounted. Used only
    # for logging/downstream Hz-to-mm conversion (see report_probe_gradient()
    # above and leonmr/xpulse_imaging.py) -- has no effect on this pp's own
    # timing/hardware calls. ProbeList pairs with Probe (same convention as
    # GradAxisList/GradAxis) to give SpinFlow the actual dropdown choices.
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet [HFX(default,calibrated)/Low Gamma(not yet calibrated)]")
    ProbeList = Parameter("ProbeList", "HFX,Low Gamma", ParameterTypes.String, "Probe Options")

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
    RValueRatio = R_REBURP_180 / R_EBURP2_90
    if int(P.AutoMatchSlice) != 0:
        P180sh_eff = P.P90sh * RValueRatio
        comms.log("AutoMatchSlice=1: P180sh overridden to {0:.2f} us (from "
                  "P90sh={1}us x R_REBURP_180/R_EBURP2_90={2:.4f}) so the "
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

    # ---- Excited slice width: CALCULATED AND RECORDED every run, for BOTH
    # pulses ------------------------------------------------------------
    ExBW_Hz, ExWidth_cm, GradTotal = compute_excited_width(
        P, comms, P.P90sh, R_EBURP2_90, "EXCITATION (90)")
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

    report_probe_gradient(P, comms)

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

            # ---- Excitation: E-BURP-2 + slice-select gradient + rephase --
            # (sequential composition -- see GradEventTime design note above)
            Delay(P.PreGrad)
            slew_rate_fn(1.0, abs(P.G1/P.RampTime))
            gradient_fn(P.RampTime, P.G1)
            shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, P.TXAmplitude90, P.TXEnableTime, Frequency, PulseOffsetHz, ExBurpA, ExBurpB)
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
            shaped_pulse(P180sh_eff, ph["PH2"], P.RefocusShape, P.TXAmplitude180, P.TXEnableTime, Frequency, PulseOffsetHz, RefocusBurpA, RefocusBurpB)
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
# 2. Claude - 07/08/26 - Simplified the power calculation: removed the
#    three-way PowerCalcMethod calculator and its RefAmplitude_HP/
#    HardPulseWidth/ExcitationRotation/RefocusRotation/PowerAdjust_dB/
#    Excitation_dB/Refocus_dB/DbConvention Parameters entirely. RFAsh0/
#    RFAsh1 now each have a single module-level computed default
#    (_DEFAULT_RFASH0/_DEFAULT_RFASH1), approximated from the calibrated
#    hard-90 and each shape's own (coherent) integration factor -- still
#    only a starting point, calibrate via nutation before trusting. Switched
#    both shapes from Fourier-coefficient reconstructions (E-BURP-1/RE-BURP)
#    to the exact 1000-point E-BURP-2/RE-BURP amplitude/phase tables
#    transcribed from real Bruker TopSpin shape-file exports
#    (_EBURP2_TABLE_B1/_REBURP_TABLE_B1, resampled via
#    _resample_burp_table()); EBURP1/BURP coefficient reconstruction kept as
#    a fallback option for the excitation shape. Fixed
#    shape_integration_factor() to use the COHERENT (signed) mean instead of
#    mean(|amp|) -- see LWG_Selective-PulseAcquire_H.py's changelog for the
#    full validation story. Replaced R_EBURP1_90/R_REBURP_180 (from-scratch
#    Bloch-simulated bandwidth-time products for the old EBURP1/RE-BURP
#    reconstructions, 4.9/5.9) with R_EBURP2_90=4.952/R_REBURP_180=5.814,
#    taken directly from each shape file's own ##$SHAPE_BWFAC metadata for
#    the exact tables now driving hardware -- the vendor's own numbers, not
#    a re-derivation. P180sh's own default is now DERIVED at import time
#    from P90sh's default via this new R-value ratio (previously hand-
#    computed and hardcoded as 6020.0), so the slice-matching logic in
#    run() stays self-consistent with the Parameter panel's displayed
#    defaults automatically.
# 3. Claude - 10/08/26 - FIXED time_calculation(): the reported 'Sequence
#    Time' was a bare P90sh+P180sh+Tau*2+points*DW guess that omitted
#    RecycleDelay (RD, by far the dominant term) AND the entire gradient/
#    rephase/crusher event (found while auditing every pp file's
#    time_calculation() after you reported the reported time never
#    reflecting real parameters). Rewrote as an exact term-for-term sum of
#    every Delay()/duration argument in the (purely sequential) per-scan
#    block, then algebraically simplified (via symbolic expansion) --
#    RefocusPulseWidth, ExGradHoldRephase (and hence RephaseFraction), and
#    CrusherBlockTime ALL cancel out of the total completely (crusher
#    timing is CrusherOn-independent: each crusher's own elapsed time is
#    exactly compensated by being subtracted out of the following TAU,
#    whether or not it actually runs), leaving RD + 2*Tau + P90sh/2 +
#    PreGrad + RampTime + Dead1 + 2.5*TXEnableTime + acquisition + a small
#    fixed remainder. Verified via the mock harness (4354296.0us at
#    default Parameters, matching the symbolic derivation exactly).
#
# 4. Claude - 18/08/26 - Added FULL Probe support at your request:
#    this file had GradAxis/GradAxisList (a working axis selector) but was
#    missing the Probe Parameter entirely -- the only gradient-using
#    imaging file in this repo where that was true (found while auditing
#    "everything that uses gradients for diffusion or imaging" against
#    the Probe/ProbeList pairing added elsewhere). Added Probe/ProbeList
#    (dropdown pairing, "HFX,Low Gamma", same convention as GradAxis/
#    GradAxisList), MAXGRAD_TABLE, and report_probe_gradient() (reads the
#    axis from P.Axis, matching the imaging family's own pattern -- e.g.
#    LWG_CPMG-Image-Echo_H.py), called once in run(). Also added Probe
#    and Dead1 to sequence_basic(). No change to timing/hardware calls --
#    verified via the mock harness (Probe/gradient calibration now logs
#    correctly).
#
# -----------------------------------------------------------------------------
