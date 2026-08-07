#-------------------------------------------------------------------------------
# Name:        LWG_Slice-Selective-PulseAcquire_X.py
# Purpose:     TRUE spatial SLICE-selective pulse-acquire (single shaped
#              E-BURP-2 excitation, gradient ON during the pulse), on the
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
# Revised:     07/08/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.1
#
# X-CHANNEL CALIBRATION -- IMPORTANT: this is a mechanical H->X port (all
# Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1). SF
# now defaults to 15.01 MHz as a placeholder -- P90sh/RFAsh0/GradCal_HzPerCm/
# R_EBURP2_90 below are still the 1H-tuned values/assumptions from the _H
# file and are almost certainly WRONG for whatever nucleus you actually put
# on the X channel (in particular, GradCal_HzPerCm and any gyromagnetic-
# ratio-dependent assumption must be re-measured/re-derived for your target
# nucleus -- R_EBURP2_90 itself is nucleus-independent, being a property of
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
#    envelope is: rephase area = 0.5 x (excitation-gradient area). E-BURP-2
#    is knowingly ASYMMETRIC in time (same family/rationale as E-BURP-1 --
#    see the design notes in LWG_1D-Image-Echo-Selective_H.py changelog item
#    4d), so 0.5 is a physically-motivated STARTING point, not an exact
#    value for this specific shape. RephaseFraction (Parameter, default 0.5)
#    is exposed so you can tune it empirically: on a real sample/phantom,
#    scan RephaseFraction and maximise the acquired FID's initial
#    amplitude -- that's your true optimum.
#  - EXCITED SLICE WIDTH: R_EBURP2_90 (the excitation bandwidth-time
#    product) is taken DIRECTLY from the real Bruker TopSpin EBurp2 shape
#    file's own ##$SHAPE_BWFAC= 4.952000E00 (see _EBURP2_TABLE_B1's own
#    comment for the source and cross-validation) -- no re-derivation
#    needed, since this is the vendor's own number for the exact shape now
#    driving hardware. (An earlier version of this file, when the default
#    shape was a hand-reconstructed E-BURP-1, Bloch-simulated its own R
#    value -- got R=4.96, reassuringly close to this file's current R=4.952
#    for the related E-BURP-2 shape, but that was for a different pulse and
#    has been superseded.) This is duration-independent (R = bandwidth x
#    duration is a fixed property of a given pulse SHAPE -- the Bloch
#    equation is invariant under jointly rescaling time and all rates, so R
#    does not need to be re-simulated for every P90sh you might use).
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

# ---- BURP bandwidth-time product -----------------------------------------
# R = (excitation bandwidth [Hz]) x (pulse duration [s]), FWHM definition.
# R is duration-independent for a fixed pulse SHAPE (rescaling time and all
# rates together leaves the Bloch equation, and therefore R, unchanged) --
# it does not need to be re-derived for whatever P90sh you actually use.
#
# Taken DIRECTLY from the same Bruker TopSpin EBurp2 shape file
# _EBURP2_TABLE_B1 above was transcribed from (##$SHAPE_BWFAC=
# 4.952000E00), rather than re-derived by Bloch-simulating an EBURP1
# Fourier-coefficient reconstruction (the previous approach, which got
# R=4.96 -- reassuringly close, but this is now the vendor's own number
# for the exact shape actually driving hardware, not a re-derivation of a
# different shape). See file header design notes for the original
# Bloch-simulation method, still valid as an independent cross-check if
# you want to re-verify (leonmr/bloch_sim.py can measure the FWHM of the
# 'EBURP2' excitation profile directly).
R_EBURP2_90 = 4.952

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
                "E-BURP-2 excitation with the slice-select gradient ON "
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

# E-BURP-2 (90 deg excitation), as a 1000-point SIGNED b1(tau) table
# (peak-normalised to 1.0, tau uniformly spaced over [0,1)) -- exported
# DIRECTLY from a real Bruker TopSpin shape file (JCAMP-DX 'Shape Data',
# ##$SHAPE_PARAMETERS= Type: EBurp2, ##$SHAPE_TOTROT= 90), supplied by
# L. Gordon 07/08/2026. Used instead of a hand-entered Fourier-coefficient
# reconstruction (as EBURP1/BURP still are, via ExBurpCoeffsA/B) because
# this is numerically EXACT to Bruker's own reference shape, with no risk
# of a transcribed-coefficient error driving real hardware. amp/phase are
# recovered from this table's sign (amp=|value|, phase=0 or 180) by
# _resample_burp_table() below -- see shape_integration_factor()'s
# docstring for why the SIGNED (not absolute-value) profile matters.
#
# CROSS-VALIDATED 07/08/2026: recomputing shape_integration_factor() on
# this exact table reproduces the shape file's own quoted
# ##$SHAPE_INTEGFAC= 6.102960E-02 to 6 significant figures -- independent
# confirmation both that this table was transcribed correctly AND that
# shape_integration_factor()'s coherent/signed-mean definition is the same
# one Bruker itself uses (see LWG_Selective-PulseAcquire_H.py's changelog
# for the full validation story).
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


def _resample_burp_table(table, n_steps):
    """Resample a fixed reference SIGNED b1(tau) table (uniformly spaced
    over [0,1)) onto n_steps points via linear interpolation of the SIGNED
    profile -- NOT of separately-folded amp/phase, which would create
    spurious intermediate phase values (e.g. ~90 deg) right at a 0/180
    sign-flip boundary. Returns (amp, phase_deg) in the same convention as
    the rest of generate_shape()."""
    table = np.asarray(table, dtype=float)
    tau_table = np.linspace(0.0, 1.0, len(table), endpoint=False)
    tau_out = np.linspace(0.0, 1.0, n_steps, endpoint=False)
    signed = np.interp(tau_out, tau_table, table)
    return np.abs(signed), np.where(signed < 0, 180.0, 0.0)


def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical machinery to LWG_Selective-Echo_H.py /
    LWG_1D-Image-Echo-Selective_H.py (see those files for the full
    docstring); duplicated here so this file stays fully self-contained.

    shape_name: 'GAUSSIAN', 'SINC', 'EBURP2' (fixed Bruker-sourced table,
    no coefficients needed), or 'EBURP1'/'BURP' (needs burp_coeffs_A)."""
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

    elif shape_name in ('EBURP1', 'BURP'):
        if not burp_coeffs_A:
            raise ValueError(
                "shape '{0}' requested but no Fourier A-coefficients were "
                "supplied (BurpCoeffsA is empty). EBURP1 is defined in "
                "Geen, H. &amp; Freeman, R., J. Magn. Reson. 93, 93-141 "
                "(1991) as a truncated Fourier series -- enter "
                "coefficients as ExBurpCoeffsA/B, or use 'GAUSSIAN' / "
                "'SINC' / 'EBURP2' instead, which need no coefficients."
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
                          "SINC, EBURP2, EBURP1/BURP (needs BurpCoeffsA/B)."
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

    # pulse_offset is in Hz (see PulseOffset Parameter); base_frequency
    # is in MHz -- *1.0e-6 converts Hz->MHz before adding. Omitting this
    # conversion sends the synth a frequency off by ~1e6x whenever
    # PulseOffset != 0, which fails hardware init (a real bug fixed
    # 06/08/2026 -- it was invisible with PulseOffset=0, the default).
    Channel2SetFrequency(1, base_frequency + pulse_offset*1.0e-6)
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

    Used once, at import time below, to turn a calibrated hard-90 into an
    approximate starting RFAsh0 default. VERIFY any shape/power
    combination with leonmr/bloch_sim.py (or a real nutation curve) before
    trusting it quantitatively -- this is still only an approximation.
    """
    peak = np.max(np.abs(amp_profile))
    if peak <= 0:
        raise ValueError("shape_integration_factor: shape has zero peak amplitude")
    signed_profile = np.asarray(amp_profile) * np.cos(np.deg2rad(phase_profile_deg))
    return float(np.mean(signed_profile) / peak)


# ---- Approximate default RFAsh0 from a calibrated hard 90 -----------------
# RFAsh0 (LP-relative) ~= HardAmplitude * HardWidth
#                          / (ShapeWidth * IntegrationFactor) / LPMaxFraction
# (target rotation is 90deg for both the hard reference and this excitation
# shape, so the two 90s in the full Bruker-style formula -- see
# LWG_Selective-Echo_H.py's changelog for the original derivation -- cancel
# out.) This is only a STARTING POINT -- calibrate via nutation on a real
# sample before trusting it quantitatively.
_HARD_P90_WIDTH_US = 9.58          # this repo's usual calibrated hard-90 width (P90/P1Hard)
_HARD_P90_AMPLITUDE = 0.4          # ...at this HP-relative amplitude (RFA0)
_LP_MAX_FRACTION = 0.10            # measured LP-vs-HP max-power ratio (LPMaxFraction default, below)
_DEFAULT_EX_SHAPE_NAME = 'EBURP2'  # exact Bruker-sourced table -- see _EBURP2_TABLE_B1 above
_DEFAULT_EX_SHAPE_WIDTH_US = 5000.0  # matches P90sh's own default, below
# Only used as ExBurpCoeffsA/B's own Parameter defaults (i.e. if you switch
# ExShape to 'EBURP1'/'BURP') -- NOT used for EBURP2, which is a fixed
# table, not a coefficient reconstruction.
_DEFAULT_EXBURP_A = "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00"
_DEFAULT_EXBURP_B = "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01"

_default_ex_amp, _default_ex_phase = generate_shape(
    _DEFAULT_EX_SHAPE_NAME, int(_DEFAULT_EX_SHAPE_WIDTH_US),
    parse_burp_coeffs(_DEFAULT_EXBURP_A), parse_burp_coeffs(_DEFAULT_EXBURP_B))
_DEFAULT_RFASH0 = round(
    _HARD_P90_AMPLITUDE * _HARD_P90_WIDTH_US
    / (_DEFAULT_EX_SHAPE_WIDTH_US * shape_integration_factor(_default_ex_amp, _default_ex_phase))
    / _LP_MAX_FRACTION, 4)


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
    R-value bandwidth-time product (see R_EBURP2_90 above) and the
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
    # E-BURP-2, exact Bruker-sourced table (see _EBURP2_TABLE_B1 above).
    ExcitationShape = Parameter("ExShape", _DEFAULT_EX_SHAPE_NAME, ParameterTypes.String, "Excitation (90) Shape [EBURP2(default, exact Bruker table)/EBURP1/GAUSSIAN/SINC/BURP] -- see ExBurpCoeffsA/B")
    ExBurpCoeffsA = Parameter("ExBurpCoeffsA", _DEFAULT_EXBURP_A, ParameterTypes.String, "E-BURP-1 Fourier cosine coeffs A0..A8 (Geen &amp; Freeman 1991, Table 2, nmax=8) -- for ExShape=EBURP1/BURP only")
    ExBurpCoeffsB = Parameter("ExBurpCoeffsB", _DEFAULT_EXBURP_B, ParameterTypes.String, "E-BURP-1 Fourier sine coeffs B0(unused)..B8 (Geen &amp; Freeman 1991, Table 2, nmax=8) -- for ExShape=EBURP1/BURP only")
    P90sh = Parameter("P90sh", _DEFAULT_EX_SHAPE_WIDTH_US, ParameterTypes.Double, "Shaped 90&#176; Pulse Width [&#956;s] -- also sets shape resolution (1 point/&#956;s) AND (with G1/GradCal_HzPerCm) the excited slice width -- see ExcitedWidth logged every run")

    TXAmplitude90 = Parameter("RFAsh0", _DEFAULT_RFASH0, ParameterTypes.Double,
                              "Shaped 90&#176; TX Power [0.0&#8230;1.0 of LP max] -- approx. from hard-90 calib., CALIBRATE via nutation before trusting", RFA, min=0.0, max=1.0)
    LPMaxFraction = Parameter("LPMaxFraction", _LP_MAX_FRACTION, ParameterTypes.Double, "LP Port Max Power as Fraction of HP Port [0.0&#8230;1.0] -- your measured value")

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

    # ---- Excited slice width: CALCULATED AND RECORDED every run ----------
    ExcitedBW_Hz, ExcitedWidth_cm, GradTotal = compute_excited_width(
        P, comms, P.P90sh, R_EBURP2_90, "excitation")

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
            shaped_pulse(P.P90sh, ph["PH1"], P.ExcitationShape, P.TXAmplitude90, P.TXEnableTime, Frequency, PulseOffsetHz, ExBurpA, ExBurpB)
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
# 2. Claude - 07/08/26 - Mechanical port of LWG_Slice-Selective-PulseAcquire_H.py's
#    v1.1 changes (see that file's own changelog entry 2 for the full
#    rationale): removed the three-way PowerCalcMethod calculator and its
#    RefAmplitude_HP/HardPulseWidth/ExcitationRotation/PowerAdjust_dB/
#    Excitation_dB/DbConvention Parameters; RFAsh0 now has a single
#    module-level computed default (_DEFAULT_RFASH0) approximated from the
#    calibrated hard-90 and the shape's own (coherent) integration factor --
#    still only a starting point, calibrate via nutation before trusting.
#    Switched the excitation shape from a Fourier-coefficient E-BURP-1
#    reconstruction to the exact 1000-point E-BURP-2 amplitude/phase table
#    transcribed from a real Bruker TopSpin shape-file export
#    (_EBURP2_TABLE_B1, resampled via _resample_burp_table()); EBURP1/BURP
#    coefficient reconstruction kept as a fallback option. Fixed
#    shape_integration_factor() to use the COHERENT (signed) mean instead of
#    mean(|amp|). Replaced R_EBURP1_90 (a from-scratch Bloch-simulated
#    bandwidth-time product for the old EBURP1 reconstruction) with
#    R_EBURP2_90=4.952, taken directly from the Bruker shape file's own
#    ##$SHAPE_BWFAC metadata for the exact EBURP2 table now driving hardware.
#    As with the rest of this X-channel port, all 1H-tuned numeric defaults
#    (P90sh/RFAsh0/GradCal_HzPerCm) still need re-calibrating for whatever
#    nucleus you actually put on the X channel -- see X-CHANNEL CALIBRATION
#    note at the top of this file.
#
# -----------------------------------------------------------------------------
