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
# Version:     1.3
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
#  - SHAPE: default is now E-BURP-2 (90 deg excitation), synthesised from a
#    1000-point SIGNED b1(tau) table exported directly from a real Bruker
#    TopSpin shape file (see _EBURP2_TABLE_B1's own comment for the source
#    metadata and cross-validation against Bruker's own quoted
#    SHAPE_INTEGFAC) -- NOT a hand-entered Fourier-coefficient
#    reconstruction, avoiding any risk of a transcription error driving
#    real hardware. EBURP1/BURP (Fourier-coefficient-based, via
#    ExBurpCoeffsA/B) remain available as ExShape options if needed.
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

# E-BURP-2 (90 deg excitation), as a 1000-point SIGNED b1(tau) table
# (peak-normalised to 1.0, tau uniformly spaced over [0,1)) -- exported
# DIRECTLY from a real Bruker TopSpin shape file (JCAMP-DX 'Shape Data',
# ##$SHAPE_PARAMETERS= Type: EBurp2, ##$SHAPE_TOTROT= 90,
# ##$SHAPE_TYPE= Excitation), supplied by L. Gordon 07/08/2026. Used
# instead of a hand-entered Fourier-coefficient reconstruction (as EBURP1/
# BURP below still are) because this is numerically EXACT to Bruker's own
# reference shape, with no risk of a transcribed-coefficient error driving
# real hardware. amp/phase are recovered from this table's sign (amp=
# |value|, phase=0 or 180) by _resample_burp_table() below -- see
# shape_integration_factor()'s docstring for why the SIGNED (not
# absolute-value) profile matters.
#
# CROSS-VALIDATED 07/08/2026: recomputing shape_integration_factor() on
# this exact table reproduces the shape file's own quoted
# ##$SHAPE_INTEGFAC= 6.102960E-02 to 6 significant figures (0.0610296) --
# independent confirmation both that this table was transcribed correctly
# AND that shape_integration_factor()'s coherent/signed-mean definition
# (changelog item 3) is the same one Bruker itself uses, not merely a
# plausible guess.
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
    """Resample a fixed reference SIGNED b1(tau) table (e.g.
    _EBURP2_TABLE_B1, uniformly spaced over [0,1)) onto n_steps points via
    linear interpolation of the SIGNED profile -- NOT of separately-folded
    amp/phase, which would create spurious intermediate phase values (e.g.
    ~90 deg) right at a 0/180 sign-flip boundary. Returns (amp, phase_deg)
    in the same convention as the rest of generate_shape()."""
    table = np.asarray(table, dtype=float)
    tau_table = np.linspace(0.0, 1.0, len(table), endpoint=False)
    tau_out = np.linspace(0.0, 1.0, n_steps, endpoint=False)
    signed = np.interp(tau_out, tau_table, table)
    return np.abs(signed), np.where(signed < 0, 180.0, 0.0)


def generate_shape(shape_name, n_steps, burp_coeffs_A=None, burp_coeffs_B=None):
    """Synthesise a (amplitude[0..1], phase[deg]) pair of length n_steps ON
    THE FLY -- identical to LWG_Selective-Echo_H.py's version; duplicated
    here so this file stays fully self-contained.

    shape_name: 'GAUSSIAN', 'SINC', 'EBURP2' (fixed Bruker-sourced table,
    no coefficients needed), or 'EBURP1'/'BURP' (needs burp_coeffs_A --
    see ExBurpCoeffsA/B Parameters)."""
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
                "supplied. EBURP1 is defined in Geen, H. &amp; "
                "Freeman, R., J. Magn. Reson. 93, 93-141 (1991) as a "
                "truncated Fourier series -- enter coefficients as "
                "ExBurpCoeffsA/B, or use 'GAUSSIAN' / 'SINC' / 'EBURP2' "
                "instead, which need no coefficients."
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
_DEFAULT_SHAPE_NAME = 'EBURP2'     # exact Bruker-sourced table -- see _EBURP2_TABLE_B1 above
_DEFAULT_SHAPE_WIDTH_US = 5000.0   # matches P90sh's own default, below
# Only used as ExBurpCoeffsA/B's own Parameter defaults (i.e. if you
# switch ExShape to 'EBURP1'/'BURP') -- NOT used for EBURP2, which is a
# fixed table, not a coefficient reconstruction.
_DEFAULT_EXBURP_A = "0.23,0.89,-1.02,-0.25,0.14,0.03,0.04,-0.03,0.00"
_DEFAULT_EXBURP_B = "0.00,-0.40,-1.42,0.74,0.06,0.03,-0.04,-0.02,0.01"

# Precomputed 07/08/2026 from the formula above (generate_shape('EBURP2',
# 5000, ...) -> shape_integration_factor(...) = 0.0610296, matching the
# EBurp2 shape file's own SHAPE_INTEGFAC to 6 s.f. -- see the changelog).
# Deliberately NOT evaluated at import time: SpinFlow executes this entire
# file just to LOAD the sequence into the parameter panel (before run() is
# ever called), so any on-the-fly shape synthesis at module scope runs on
# EVERY load, not just every scan -- a needless dependency for a value
# that's meant to be swept by hand anyway. generate_shape()/
# shape_integration_factor() remain fully available and are used normally
# inside run(); only this one-time calibration default is now a literal.
_DEFAULT_RFASH0 = 0.1256


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
# 4. Claude - 07/08/26 - SWITCHED DEFAULT SHAPE TO E-BURP-2: added
#    _EBURP2_TABLE_B1 (1000-point signed b1(tau) table) and
#    _resample_burp_table(), transcribed directly from a real Bruker
#    TopSpin EBurp2 shape-file export you supplied -- generate_shape()'s
#    'EBURP2' path now resamples this exact table (linear interpolation of
#    the SIGNED profile, avoiding spurious intermediate phase values at
#    sign-flip boundaries) instead of reconstructing from hand-entered
#    Fourier coefficients, removing any transcription-error risk for the
#    shape actually driving hardware. _DEFAULT_SHAPE_NAME switched from
#    'EBURP1' to 'EBURP2' -- _DEFAULT_RFASH0 recalculated automatically
#    (integration factor now 0.061030, vs EBURP1's 0.067516; both now
#    computed via the item-3 coherent/signed fix). CROSS-VALIDATION:
#    shape_integration_factor() on this table reproduces the Bruker shape
#    file's own quoted SHAPE_INTEGFAC (6.102960E-02) to 6 significant
#    figures -- confirms both the table transcription and the item-3 fix
#    against an independent, vendor-authoritative source. EBURP1/BURP
#    (coefficient-based) remain available as ExShape options.
#
# -----------------------------------------------------------------------------
