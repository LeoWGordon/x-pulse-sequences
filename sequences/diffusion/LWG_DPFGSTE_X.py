#-------------------------------------------------------------------------------
# Name:        LWG_DPFGSTE_X.py
# Purpose:     Double Pulsed-Field Gradient STimulated Echo (DPFGSTE) for
#              convection-compensated translational-diffusion measurement on
#              the {X} channel: two stimulated-echo diffusion
#              blocks in series (P1..P5, four G1 gradients, two DELTA/2
#              storage periods) followed by an LED eddy-current delay
#              (P6/P7) before acquisition. X-Pulse Broadband Benchtop NMR
#              Spectrometer. NO IMAGING / NO SLICE SELECTION.
#
# Reference:   Jerschow &amp; Mueller, J. Magn. Reson. 1997, 125, 372-375
#              (suppression of convection artifacts, double STE)
#              Stejskal &amp; Tanner, J. Chem. Phys. 1965, 42, 288 (b-value)
#              Price &amp; Kuchel, J. Magn. Reson. 1991, 94, 133 (trapezoid b)
#
# Authors:     Claude, for L. Gordon (DTU) -- ported from Oxford Instruments'
#              DPFGSTE_X.py (RJB), which remains the timing
#              reference for this file.
#
# Created:     15/09/2026
# Copyright:   (c) Oxford Instruments Magnetic Resonance, 2013-
# Version:     1.0
#
# status: draft
# validated_by: -
# validated_date: -
#
# Design notes:
#  - STRUCTURALLY A VERBATIM PORT of the vendor's DPFGSTE_X.py.
#    WHY YOU WOULD USE THIS OVER LWG_PFGSTE_X.py: convection. In a sample
#    with a vertical temperature gradient (anything heated or cooled, and
#    low-viscosity liquids even at ambient), convective flow adds a
#    coherent velocity phase that a single stimulated echo reads as extra
#    "diffusion" -- systematically inflating D, and often showing up as a
#    non-exponential or oscillating decay. Two STE blocks back to back, with
#    the second block's effective wavevector opposite to the first, cancel
#    that velocity phase to first order while the diffusion attenuation from
#    both blocks adds. The cost is four 90s more, three storage periods of
#    T1 decay instead of one, and a further ~2x signal loss.
#  - The LED (longitudinal eddy-current delay) tail, P6 - g8 homospoil - d2
#    - P7, parks magnetization along z while the last gradient's eddy
#    currents decay, so acquisition starts on a settled field. q is already
#    back to 0 by then, so it adds no diffusion weighting -- see
#    bvalue_segments().
#  - DELTA IS THE TOTAL diffusion time, split as DELTA/2 per storage period
#    (the vendor's own D74 semantics: d1 = DELTA/2 - fixed events, applied
#    twice). This differs from LWG_PFGSTE_X.py, where the single storage
#    period gets all of DELTA -- so the same DELTA in both files does NOT
#    give the same b. Compare the logged b-values, not the DELTA settings.
#  - GRADIENT POLARITY: all four diffusion gradients (g1, g3, g4, g7) use
#    the SAME physical sign and amplitude (G1); their EFFECTIVE signs
#    alternate +,-,-,+ via the storage pulses. Do NOT flip any of them by
#    hand. See bvalue_segments() for the full wavevector walk-through.
#  - X-CHANNEL PORT -- IMPORTANT, READ BEFORE THE FIRST RUN. The pulse
#    train is a mechanical H->X port of LWG_DPFGSTE_H.py (Channel1/Transmit1/
#    Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1), but the Parameter
#    block is NOT: it follows the vendor's own _X file, which differs from
#    its _H sibling by more than a channel swap.
#      * The observe-channel SHORT CODES are the vendor's _X spellings --
#        SFX, O1X, TxPPMX, P90X, RFAX -- not SF/O1/TxPPM/P90/RFA0. This
#        deliberately does NOT follow LWG_PGSTE_X.py/LWG_1D-Image-Echo_X.py,
#        which normalise them to the _H spellings: SpinFlow decides which
#        channel and nucleus a sequence is valid for from its parameter set
#        (see the "Sequence invalid for selection: ... Observed: 1H" errors
#        in docs/error-logs/SpinFlow_error_logs), and the vendor spelling is
#        the one known to load on the X channel. If you would rather have
#        the house spelling, change it deliberately and re-test loading --
#        do not assume it is cosmetic.
#      * NBlock (default 4) is present, because the vendor's _X CallBack1D
#        -- kept verbatim here -- genuinely accumulates across scans and
#        honours it, unlike the _H one.
#      * Filter, NS and P90X keep the vendor's _X defaults (which differ
#        from the _H file's), not the _H ones.
#    CALIBRATION: SF defaults to the vendor's 15.01 MHz and Nucleus to
#    13C, which are consistent with each other on this magnet (10.7084
#    MHz/T x ~1.40 T). P90X, RFAX and RA are still vendor defaults and are
#    almost certainly WRONG for whatever you actually put on the X channel.
#    Recalibrate SFX, P90X, RFAX, RA -- and set Nucleus, which b depends on
#    as gamma^2 -- before running this on a sample.
#  - HOUSE-STYLE CHANGES vs. the vendor file (timing is UNCHANGED; every
#    Delay()/duration expression is algebraically identical to the vendor's,
#    including its small -1/-2/-3/-5 us trim constants -- the point of this
#    port was to keep the vendor's proven timing and change only naming,
#    reporting and safety):
#      * PHASE-CYCLE PARAMETERS RENAMED, and this is a REAL BUG FIX, not
#        cosmetics: the vendor declares them as P1Phase = Parameter("PH1",
#        ...), i.e. Python attribute name != SpinFlow short code. run() then
#        builds its own PhaseListContainer per pulse, which happens to work
#        -- but the moment anything routes through PhasesManager (the house
#        convention used by every other sequence here) that mismatch raises
#        "KeyError: 'PH1'", which is exactly what this instrument did on
#        14/08/2026. Attribute name == short code == Phases.Incd() key for
#        every phase Parameter now (PH1/PH2/.../PHRX), and phases come from
#        a single PhasesManager(P). See the repo's pulse-programme-parameters
#        skill, rules 7 and 8.
#      * D74 (DiffDelta) -> DELTA and D71 (GradientOnTime) -> delta, so the
#        two Stejskal-Tanner quantities are named after the symbols they
#        actually are, matching LWG_PGSTE_H.py. sequence_basic() updated to
#        match. NOTE: renaming a short code orphans any existing .par file
#        for the vendor sequence -- SpinFlow will fall back to the defaults
#        below the first time you load this one.
#      * Gx_max/Gy_max/Gz_max REMOVED, replaced by the house Probe selector
#        + MAXGRAD_TABLE + FPX/FPY/FPZ. This is the "record the probe"
#        change you asked for, and it has a REAL NUMERICAL CONSEQUENCE:
#        the vendor defaults (Gx_max = Gy_max = 0.01 T/m, Gz_max = 0.05 T/m)
#        are ~12x SMALLER than this repo's own measured calibration
#        (HFX: x = 11.879, y = 11.978, z = 57.915 G/cm = 0.1188, 0.1198,
#        0.5792 T/m). Whichever is right, they cannot both be -- and since
#        b scales as G^2, the two differ by ~140x in b. The measured
#        MAXGRAD_TABLE values are used here because they are this repo's
#        own, are already what leonmr/xpulse_imaging.py uses for Hz->mm
#        conversion, and were measured on this magnet; the vendor numbers
#        look like untouched placeholder defaults. VALIDATE THIS against a
#        known-D standard (e.g. H2O/D2O at a known temperature) before
#        trusting any absolute b-value or D from this sequence.
#      * safe_delay() around every COMPUTED delay, so a DELTA (or LED delay)
#        too short for the fixed events inside it fails with a message
#        naming the delay, instead of a negative Delay() argument.
#      * apply_gradient() replaces the vendor's four-way copy-pasted
#        if axis == 1/2/3/else block at every gradient (7 copies in
#        DPFGSTE). It emits the identical instruction sequence --
#        SlewRate(1us), ramp up, plateau, ramp down, explicit off(1us),
#        total 2 + 2*RampTime + plateau -- and the GradAxis='none' branch is
#        a pure Delay() of that same total, so switching axis can never
#        change the timing. GradientMatrix is now diag(FPX,FPY,FPZ) (house
#        convention) instead of the vendor's per-axis 0/1 selector matrix.
#      * Probe/ProbeList, Nucleus/NucleusList, mains-lock trigger, and
#        duty-cycle guard rails added, per the rest of this family.
#
#  - B-VALUE REPORTING (the other thing you asked for). Three separate
#    routes, because the pulse programme CANNOT write a file where you
#    wanted it -- see the next note:
#      1. The b-value for the current G1 is computed at sequence start and
#         written to the SpinFlow log, along with the probe gradient
#         calibration and the gamma it used.
#      2. The whole b-value table for the GradList Parameter (your relative
#         gradient list) is logged too, one row per gradient -- so a
#         gradient-ramp experiment's b-values are all there in the log
#         before the first scan runs. Set ShowBTable = 0 to suppress.
#      3. Both go into the saved dataset's JCAMP metadata as jc_bvalue,
#         jc_bvalue_smm2, jc_gamma, jc_gradcal_Tpm, jc_gradlist and
#         jc_bvalue_list, alongside the vendor's existing jc_grad /
#         jc_bigDELTA / jc_smallDELTA keys. UNVERIFIED: whether the
#         vendor's JCamp1D_Diffusion.txt template silently ignores keys it
#         does not know. If a dataset fails to save or looks malformed,
#         delete the jc_data.extend(bvalue_jcamp) line in jcamp_meta() --
#         routes 1 and 2 are unaffected.
#    The b-value itself is NOT the textbook rectangular-pulse formula. It
#    is the general definition b = integral of q(t)^2 dt, with q(t) =
#    gamma * integral of g_eff dt built from the sequence's OWN timing
#    expressions -- so the gradient ramps enter as real trapezoids, the RF
#    pulse widths and settle delays are counted where they actually fall,
#    and the storage periods contribute frozen-q terms. See
#    bvalue_segments() and integrate_bvalue(). It was checked against both
#    standard closed forms and reproduces each to machine precision -- see
#    the changelog for the exact checks. WATCH THE delta CONVENTION when
#    you compare against a textbook: the quantity that belongs in those
#    formulae is the HALF-HEIGHT duration delta_eff = delta + RampTime (the
#    trapezoid's area is G*(delta + RampTime), one ramp time, not two), and
#    DELTA there means leading edge to leading edge -- start of one ramp-up
#    to start of the next -- which is NOT the same as this file's DELTA
#    Parameter. Using the plateau delta alone in a trapezoid formula
#    underestimates b by ~20% at the default RampTime=500 us, delta=4000 us.
#
#  - THE b-VALUE IS FACTORISED, SO IT SURVIVES A RECALIBRATION. b depends
#    on the gradient calibration, and that calibration is exactly the thing
#    currently in doubt (see the MAXGRAD_TABLE note above). So rather than
#    only recording a b that is hostage to today's table, the sequence also
#    reports the parts that are NOT:
#
#        b = gamma^2 * G_abs^2 * S0,   G_abs = G_rel * FP * GradMax
#
#      * S0 (jc_b_shape_s3, in s^3) is the integral of the UNIT-AMPLITUDE
#        waveform -- pure timing. No gamma, no probe, no calibration, no FP.
#        Change the gradient table, change the nucleus, change nothing but
#        those, and S0 is still correct.
#      * b/G_abs^2 (jc_b_per_G2) folds in gamma only. Multiply by
#        (G_rel*FP*GradMax_T_per_m)^2 with whatever calibration you trust
#        at analysis time.
#      * b/G_rel^2 (jc_b_per_Grel2) folds in everything including today's
#        GradMax, and is the number that multiplies your relative-gradient
#        list directly: b_i = (b/G_rel^2) * G_rel_i^2.
#
#    All three are CONSTANT across a whole gradient ramp -- one number plus
#    the relative-gradient list rebuilds every b in the experiment, which is
#    also why this works whether SpinFlow gives each array point its own
#    dataset or lumps the ramp into one.
#
#    WATCH THE POWER: the invariant is b/G_rel SQUARED. b/G_rel (first
#    power) is not constant -- across the default 0.05..0.95 ramp it varies
#    by the full 20x range of the ramp itself, so it is not a usable
#    normalisation.
#
#    GOING THE OTHER WAY: because b/G_rel^2 is constant, a ramp on a sample
#    of known D gives ln(S/S0) = -(b/G_rel^2)*D*G_rel^2, i.e. a straight
#    line in G_rel^2 whose slope MEASURES b/G_rel^2 with no calibration
#    assumed anywhere. Divide by gamma^2*S0*FP^2 and take the square root
#    and you have GradMax itself. That is the cleanest way to settle the
#    ~12x vendor-vs-table disagreement, and
#    tools/spinflow_bvalue_diff.py --solve-gmax <slope> --known-d <D>
#    does the arithmetic.
#
#  - GradMax PARAMETER: set it nonzero to override MAXGRAD_TABLE for the
#    selected axis without editing this file. Being a Parameter, it is saved
#    in the .par/JCAMP, so the b-value stays reproducible from the data.
#
#  - WHY THERE IS NO .diff FILE WRITTEN BY THIS SEQUENCE. You asked for the
#    b-value table to land in the experiment's parameter folder. It cannot
#    be done from here: run() does NOT execute on the Windows PC. The
#    instrument's own log (docs/error-logs/SpinFlow_error_logs) shows
#    SpinFlow copying the pulse programme to /tmp/tmpXXXXXX.py and running
#    it under /home/fire/PythonMonitor/ on the spectrometer console, while
#    the parameter files live at C:\Users\Public\Documents\SpinFlow\
#    Parameters\<serial>\User\<nucleus>\<Sequence>.par on the host. A
#    file written by run() would land on the console's filesystem, not next
#    to your parameters. Instead: run tools/spinflow_bvalue_diff.py on the
#    SpinFlow PC. It reads the saved <Sequence>.par, recomputes the table
#    with the IDENTICAL maths, and writes <Sequence>.diff beside it. It
#    does NOT re-implement any of this -- it imports bvalue_segments(),
#    integrate_bvalue(), MAXGRAD_TABLE and NUCLEUS_GAMMA_MHZ_PER_T straight
#    out of THIS file under a stub firebird, so there is nothing to keep in
#    sync: edit the maths here and the tool follows.
#
#  - NUCLEUS / GAMMA. b depends on gamma^2, so the b-value needs to know
#    what you are observing. The Nucleus Parameter (with NucleusList as its
#    dropdown, per the house GradAxis/GradAxisList convention) picks gamma
#    from a table; GammaOverride (gamma/2pi in MHz/T, 0 = use the table)
#    covers anything not in it. As a cross-check, the implied field
#    B0 = SF / (gamma/2pi) is logged, and a value outside 1.0-2.5 T warns --
#    that catches a Nucleus left set to the wrong isotope, which would
#    otherwise silently scale every b-value.
#
#  - Mains-lock trigger defaults OFF (UseMainsLock=0), per house convention.
#
# Changes/Modifications: At end of file.
#-------------------------------------------------------------------------------

from firebird import *
from firebird.applications import *
import numpy as np
import time

global BLP
BLP = 1

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
#   Callback to handle one-dimensional sequences, _X and _XH (NBlock compatible)
    def __init__(self, Params, comms, ReceiverFilter, jcamp_meta, template):

        self.comms = comms                                      # comms to the clients
        self.P = Params
        self.NS = self.P.NumScans
        self.NBlock = self.P.NumScans                           # default NBLOCK to NS
        self.jcamp_meta = jcamp_meta
        self.template = template
        self.start_time = time.time()
        self.end_time = 0
        if getattr(Params, "NBlock", 'nofound') != "nofound":
            self.NBlock = getattr(Params, "NBlock")             # if defined, set to the value defined
            if self.NBlock > self.NS:
                self.NBlock = self.NS
                self.comms.log("Error: NBlock > NS. Defaulting to NBlock == NS.")
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

        self.dummys = range(self.DS)                            # dummy scan
        self.first_scans, self.scan_to_save = get_permutation(self.DS, self.NS, self.NBlock)
        self.final_scan = (self.DS + self.NS) - 1               # last scan in sequence

    def process_data(self, scan, data):
        scan_data, clipped = data
        scaled_data = np.int32(np.round(scan_data / self.rcv_filter.gain))

        if (clipped.value == 1):                                # clipped
            self.comms.log(">>>>> !! Data is Clipped!! <<<<<")
        elif (clipped.value == 100):
            self.comms.log(">>>>> !! IPC failed!! <<<<<")

        if scan in self.first_scans or scan in self.dummys:
            self.acc_data = np.array(scaled_data)               # this assume dummy scan doesn't need accumulation
        else:
            self.acc_data += scaled_data

        if scan in self.dummys:                                 # sort out dummy scan
            self.comms.send_data(self.acc_data, scan+1, self.times, is_clipped=False,
                is_last_update=False, metadata={"DummyScan":"true"})
        elif (scan not in self.scan_to_save):
            self.comms.send_data(self.acc_data, scan+1, self.times, is_clipped=False,
                                 metadata={"JCAMP":self.jcamp_meta, "jc_template":"{0}".format(self.template), "nosave":"true", "Preview":4})
            if scan==self.final_scan:
                self.comms.send_data(self.acc_data, scan+1, self.times, is_clipped=False, is_last_update=(scan==self.final_scan),
                                metadata={"JCAMP":self.jcamp_meta, "jc_template":"{0}".format(self.template)})
        else:                                                   # scans_per_block'th scan, this data will be saved accumulator will be reset
            self.comms.send_data(self.acc_data, scan+1, self.times, is_clipped=False, is_last_update=(scan==self.final_scan),
                            metadata={"JCAMP":self.jcamp_meta, "jc_template":"{0}".format(self.template), "Preview":4})
            self.acc_data = None
        if scan == (self.DS+self.NS-1):
            self.end_time = time.time()
            self.comms.log("seqTime LWG_DPFGSTE_X Estimated: {0}, Real: {1}".format(1e-6*time_calculation(Parameters), self.end_time-self.start_time))


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
    silent/confusing failure if a computed delay would be negative. House
    convention; see LWG_PGSTE_H.py."""
    if value < 0:
        msg = ("FATAL TIMING ERROR computing '{0}': delay would be {1:.2f} us "
               "(negative). Increase DELTA, or decrease P90/PreGrad/delta/"
               "GradSpoil/RampTime/GradSettle, then retry.").format(label, value)
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
                  "({1:.2%}). Consider a longer RD or shorter P90, or confirm "
                  "with Oxford Instruments that this is within the "
                  "transmitter's rated duty cycle before running unattended."
                  .format(rf_duty, P.MaxRFDuty))
    if grad_duty > P.MaxGradDuty:
        comms.log("WARNING: estimated gradient duty cycle {0:.2%} exceeds "
                  "MaxGradDuty ({1:.2%}). Consider a longer RD or shorter "
                  "delta/GradSpoil, or confirm with Oxford Instruments that "
                  "this is within the gradient amplifier's rated duty cycle "
                  "before running unattended.".format(grad_duty, P.MaxGradDuty))


# Per-axis max gradient strength (G/cm at |G1| = 1.0, FP scaler = 1.0) for
# probes that can be fitted to this magnet. There is no way to detect which
# probe is mounted from software, so this is a MANUAL selector (Probe
# Parameter below) -- the operator must set it to match what is actually on
# the magnet. KEEP THIS TABLE IN SYNC with leonmr/xpulse_imaging.py's own
# GRADIENT_CALIBRATION dict (duplicated rather than imported, since pulse
# programmes can't import external Python modules -- same house convention
# used for safe_delay()/mains_lock_trigger() across this family).
# tools/spinflow_bvalue_diff.py needs no such syncing: it imports this table
# from this file rather than keeping its own copy.
MAXGRAD_TABLE = {
    "HFX": {"x": 11.879, "y": 11.978, "z": 57.915},   # G/cm, measured/averaged 06/08/2026
    "LOWGAMMA": {"x": 18.678, "y": 16.038, "z": 58.112},  # G/cm, PROVISIONAL avg 18-19/08/2026 (z=PGSE only; x/y=PGSTE; more replicates incl. smaller-delta PGSTE-z pending)
}

# Probe entries in MAXGRAD_TABLE whose calibration is PROVISIONAL -- real
# measured numbers, but not yet settled. Keyed like MAXGRAD_TABLE; the value
# is the caveat report_probe_gradient() surfaces at run time, so it reaches
# the SpinFlow console instead of living only in a comment nobody running the
# instrument will ever read. Delete a probe's entry once its calibration is
# confirmed. KEEP IN SYNC with leonmr/xpulse_imaging.py if that grows an
# equivalent.
PROVISIONAL_CALIBRATION = {
    "LOWGAMMA": "PROVISIONAL avg 18-19/08/2026 -- z from PGSE only, x/y from PGSTE; more replicates (incl. smaller-delta PGSTE-z) pending",
}

GAUSS_PER_CM_TO_T_PER_M = 0.01   # 1 G/cm = 1e-4 T / 1e-2 m = 1e-2 T/m


# Gyromagnetic ratios as gamma/2pi in MHz/T (IUPAC values; sign dropped
# where negative, since b depends on gamma^2 only). Used ONLY for the
# b-value report -- has no effect on timing or hardware calls.
NUCLEUS_GAMMA_MHZ_PER_T = {
    "1H": 42.5775, "2H": 6.5359, "6LI": 6.2661, "7LI": 16.5465,
    "11B": 13.6630, "13C": 10.7084, "15N": 4.3173, "17O": 5.7743,
    "19F": 40.0776, "23NA": 11.2688, "27AL": 11.1031, "29SI": 8.4654,
    "31P": 17.2351, "35CL": 4.1765, "119SN": 15.9660, "133CS": 5.6234,
    "195PT": 9.2920, "207PB": 8.8818,
}


def normalise_key(value):
    """Normalise a free-text Parameter string (Probe / Nucleus) to a table
    key: strip spaces, hyphens, slashes and underscores, then uppercase."""
    return str(value).upper().replace('/', '').replace('-', '').replace(' ', '').replace('_', '')


def get_gradient_functions(axisstr):
    """Map a 'x'/'y'/'z' GradAxis string to the (SlewRate, Gradient) hardware
    call pair for that physical axis; 'none' returns (None, None), which
    apply_gradient() renders as a pure Delay() of the same duration (the
    vendor sequences' own GradAxis='none' behaviour, preserved verbatim).
    Matches the imaging family's own get_gradient_functions() (e.g.
    LWG_CPMG-Image-Echo_H.py) otherwise."""
    if axisstr == 'none':
        return (None, None)
    elif axisstr == 'x':
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


def apply_gradient(slew_fn, grad_fn, ramp, plateau, amp):
    """One trapezoidal gradient pulse, EXACTLY as the vendor PFGSTE/DPFGSTE
    sequences emit it: SlewRate(1us) - ramp up(ramp) - plateau(Delay) -
    ramp down(ramp) - explicit off(1us). Total elapsed = 2 + 2*ramp +
    plateau us, identical in the 'none' (no gradient hardware) branch, so
    switching GradAxis never changes sequence timing."""
    if slew_fn is None:                      # GradAxis == 'none'
        Delay(2.0 + 2*ramp + plateau)
        return
    slew_fn(1.0, abs(amp/ramp))
    grad_fn(ramp, amp)
    Delay(plateau)
    grad_fn(ramp, 0)
    grad_fn(1.0, 0)                          # turn gradient off explicitly for safety


def report_probe_gradient(P, comms, axis=None):
    """Log AND return the max gradient strength (G/cm) for the currently-
    selected Probe and diffusion-gradient axis (GradAxis), or None if
    unknown. Since Probe is a Parameter, it is auto-recorded in the
    resulting JCAMP file's SpinFlow parameter block -- letting
    leonmr/xpulse_imaging.py (and tools/spinflow_bvalue_diff.py) recover the
    gradient calibration later without the operator having to remember which
    probe was fitted for a given experiment. Warns (does not fail) if the
    fitted probe's selected axis has no calibration yet."""
    axis_key = str(axis if axis is not None else getattr(P, "Axis", "z")).strip().lower()
    if axis_key == 'none':
        comms.log("Probe: {0}. GradAxis='none' -- no diffusion gradient is "
                  "applied, b = 0.".format(P.Probe))
        return None
    if P.GradMaxOverride > 0.0:
        # Operator-supplied calibration wins over the built-in table, and is
        # itself a Parameter, so it lands in the saved .par/JCAMP alongside
        # everything else needed to reproduce the b-value.
        comms.log("Probe: {0}, but GradMax={1} G/cm OVERRIDES the built-in "
                  "table for the {2}-axis (= {3:.4f} T/m at |G1|=1.0, FP "
                  "scaler=1.0).".format(P.Probe, P.GradMaxOverride, axis_key,
                                        P.GradMaxOverride*GAUSS_PER_CM_TO_T_PER_M))
        return P.GradMaxOverride
    probe_key = normalise_key(P.Probe)
    if probe_key not in MAXGRAD_TABLE:
        comms.log("WARNING: unrecognised Probe='{0}'. Known probes: {1}. "
                  "b-value and Hz->mm conversion are undefined for this "
                  "probe.".format(P.Probe, sorted(MAXGRAD_TABLE)))
        return None
    maxgrad = MAXGRAD_TABLE[probe_key].get(axis_key)
    if maxgrad is None:
        comms.log("WARNING: Probe='{0}', axis='{1}' has NO gradient "
                  "calibration yet -- b-value is undefined for this data "
                  "until it is measured.".format(P.Probe, axis_key))
        return None
    comms.log("Probe: {0}. Max gradient strength ({1}-axis) = {2} G/cm "
              "= {3:.4f} T/m at |G1|=1.0, FP scaler=1.0."
              .format(P.Probe, axis_key, maxgrad, maxgrad*GAUSS_PER_CM_TO_T_PER_M))
    caveat = PROVISIONAL_CALIBRATION.get(probe_key)
    if caveat:
        comms.log("WARNING: Probe='{0}' gradient calibration is "
                  "PROVISIONAL: {1}. b-values and Hz-to-mm "
                  "conversion from this run inherit that "
                  "uncertainty.".format(P.Probe, caveat))
    return maxgrad


def axis_fp_scaler(P):
    """The GradientMatrix diagonal element actually applied to the selected
    gradient axis -- part of the real output amplitude, so it MUST enter the
    b-value (see design notes: a wrong FP scaler silently rescales b)."""
    return {"x": P.XGradNorm, "y": P.YGradNorm, "z": P.ZGradNorm}.get(
        str(P.Axis).strip().lower(), 0.0)


def gradient_amplitude_T_per_m(maxgrad_G_per_cm, fp_scaler, g_rel):
    """Absolute encoding-gradient amplitude [T/m] for a relative (DAC
    fraction) gradient value. UNITS: maxgrad is in G/cm (MAXGRAD_TABLE's own
    unit, matching leonmr/xpulse_imaging.py), g_rel and fp_scaler are
    dimensionless, and GAUSS_PER_CM_TO_T_PER_M does the single conversion --
    see the pulse-programme-parameters skill's rule 6 on checking units at
    every arithmetic combination."""
    return abs(g_rel) * fp_scaler * maxgrad_G_per_cm * GAUSS_PER_CM_TO_T_PER_M


def gamma_rad_per_s_per_T(P, comms):
    """Gyromagnetic ratio [rad s^-1 T^-1] for the Nucleus Parameter, or the
    GammaOverride Parameter (MHz/T) if it is nonzero. Also logs the magnetic
    field implied by (SF, gamma) as a cross-check -- if the Nucleus setting
    does not match what is really on the observe channel, the implied field
    comes out wrong by the ratio of the two gyromagnetic ratios, which is a
    much easier error to spot than a silently wrong b-value."""
    if P.GammaOverride > 0.0:
        gamma_over_2pi = P.GammaOverride
        source = "GammaOverride Parameter"
    else:
        key = normalise_key(P.Nucleus)
        if key not in NUCLEUS_GAMMA_MHZ_PER_T:
            comms.log("WARNING: unrecognised Nucleus='{0}'. Known: {1}. "
                      "b-value NOT calculated -- set Nucleus, or set "
                      "GammaOverride to gamma/2pi in MHz/T."
                      .format(P.Nucleus, sorted(NUCLEUS_GAMMA_MHZ_PER_T)))
            return None
        gamma_over_2pi = NUCLEUS_GAMMA_MHZ_PER_T[key]
        source = "Nucleus='{0}'".format(P.Nucleus)
    implied_B0 = P.FrequencyBaseX / gamma_over_2pi
    comms.log("Gyromagnetic ratio from {0}: gamma/2pi = {1} MHz/T "
              "-> SF={2} MHz implies B0 = {3:.4f} T ({4:.2f} MHz for 1H)."
              .format(source, gamma_over_2pi, P.FrequencyBaseX, implied_B0,
                      implied_B0*NUCLEUS_GAMMA_MHZ_PER_T["1H"]))
    if implied_B0 < 1.0 or implied_B0 > 2.5:
        comms.log("WARNING: implied B0 = {0:.4f} T is outside the plausible "
                  "range for this magnet (~1.4 T for a 60 MHz X-Pulse, ~2.1 T "
                  "for an 90 MHz HyB / X-Pulse90). Nucleus/SF are probably "
                  "inconsistent -- the b-value below is very likely wrong."
                  .format(implied_B0))
    return 2.0e6 * np.pi * gamma_over_2pi


# 3-point Gauss-Legendre nodes/weights on [-1, 1]. Over one segment of the
# waveform below the dephasing wavevector q(t) is a QUADRATIC polynomial in
# t, so q(t)^2 is quartic (degree 4) and 3-point Gauss-Legendre integrates
# it EXACTLY (exact to degree 5) -- no sampling/discretisation error, and no
# large arrays to build on the console.
_GL3_NODES = (-0.7745966692414834, 0.0, 0.7745966692414834)
_GL3_WEIGHTS = (0.5555555555555556, 0.8888888888888888, 0.5555555555555556)


def integrate_bvalue(segments, gamma):
    """Exact b-value of an EFFECTIVE gradient waveform.

    segments : list of (duration_us, g_start_T_per_m, g_end_T_per_m),
               piecewise-LINEAR in g, in time order, covering the whole
               per-scan timeline from excitation to acquisition. "Effective"
               means the coherence-pathway sign has already been folded into
               g (see the bvalue_segments() docstring in each sequence).
    gamma    : gyromagnetic ratio [rad s^-1 T^-1].

    Returns b = integral of q(t)^2 dt in s/m^2, where q(t) = gamma *
    integral of g_eff dt is the dephasing wavevector [rad/m].

    This is the general Stejskal-Tanner b-value definition, NOT the
    rectangular-pulse closed form -- it therefore includes the gradient
    RAMPS exactly (a trapezoid, not a rectangle), the RF pulse widths, and
    the frozen-q storage periods, with no approximation beyond "RF pulses
    and coherence transfers are instantaneous at their midpoints". Reduces
    ANALYTICALLY to b = (gamma*G*delta)^2 * (DELTA - delta/3) in the
    zero-ramp limit, with DELTA measured leading-edge to leading-edge --
    verified numerically against that closed form (see changelog)."""
    q = 0.0                                   # rad/m, running wavevector
    b = 0.0                                   # s/m^2
    for dur_us, g0, g1 in segments:
        T = dur_us * 1.0e-6                   # us -> s
        if T <= 0.0:
            continue
        a0 = q                                # q(t) = a0 + a1*t + a2*t^2
        a1 = gamma * g0
        a2 = gamma * (g1 - g0) / (2.0 * T)
        acc = 0.0
        for node, weight in zip(_GL3_NODES, _GL3_WEIGHTS):
            t = 0.5 * T * (node + 1.0)
            qt = a0 + a1*t + a2*t*t
            acc += weight * qt * qt
        b += 0.5 * T * acc
        q = a0 + a1*T + a2*T*T
    return b


def bvalue_shape_factor(P):
    """The purely-geometric part of the b-value: S0 = integral of
    (integral g_hat dt)^2 dt over the UNIT-AMPLITUDE effective waveform,
    in s^3.

    b factorises exactly as

        b = gamma^2 * G_abs^2 * S0                                       (1)

    because every duration in bvalue_segments() is independent of the
    gradient amplitude, so scaling g scales q linearly and b quadratically.
    S0 therefore contains ONLY this sequence's delays -- no gyromagnetic
    ratio, no probe, no gradient calibration, no FP scaler. It is the number
    to keep if you expect any of those to be revised later: b can be rebuilt
    from it at any time, with whatever calibration you trust then.

    Verified numerically: b(G)/G^2 is constant to ~3e-16 across G_rel =
    0.05..1.0, and (1) reproduces the directly-integrated b to the same
    precision (see changelog).

    NOTE this is a SQUARE law. A "b per unit relative gradient" normalised
    as b/G_rel (first power) is NOT constant across a gradient ramp -- it
    varies by the full 20x range of the ramp itself. b/G_rel^2 is the
    invariant.
    """
    return integrate_bvalue(bvalue_segments(P, 1.0), 1.0)


def report_bvalue(P, comms, gamma):
    """Log the b-value for the CURRENT G1, plus the full b-value table for
    every relative gradient in the GradList Parameter, and return the JCAMP
    metadata strings that carry those numbers into the saved dataset.

    NOTE (see design notes): run() executes on the spectrometer console, NOT
    on the Windows PC that owns the SpinFlow Parameters folder -- a pulse
    programme therefore CANNOT write a .diff file into the experiment's own
    parameter folder. The numbers are instead (a) logged here and (b)
    written into the saved JCAMP metadata; run tools/spinflow_bvalue_diff.py
    on the SpinFlow PC to turn the saved .par into a .diff file sitting next
    to it."""
    jc = []
    if str(P.Axis).strip().lower() == 'none':
        report_probe_gradient(P, comms)
        # The shape factor is a property of the TIMING, so it is still
        # meaningful (and worth recording) even though no gradient fires.
        jc.append("jc_bvalue=0.0")
        jc.append("jc_bvalue_smm2=0.0")
        jc.append("jc_b_shape_s3={0}".format(bvalue_shape_factor(P)))
        jc.append("jc_b_per_Grel2=0.0")
        return jc
    if gamma is None:
        comms.log("b-value NOT calculated (no gyromagnetic ratio).")
        return jc
    maxgrad = report_probe_gradient(P, comms)
    if maxgrad is None:
        comms.log("b-value NOT calculated (no gradient calibration for this "
                  "Probe/GradAxis).")
        return jc
    fp = axis_fp_scaler(P)
    G = gradient_amplitude_T_per_m(maxgrad, fp, P.G1)
    b = integrate_bvalue(bvalue_segments(P, G), gamma)
    comms.log("b-value at G1={0} (FP scaler {1} -> {2:.5f} T/m on the {3}-axis): "
              "b = {4:.6g} s/m^2 = {5:.6g} s/mm^2"
              .format(P.G1, fp, G, P.Axis, b, b*1.0e-6))
    jc.append("jc_bvalue={0}".format(b))
    jc.append("jc_bvalue_smm2={0}".format(b*1.0e-6))
    jc.append("jc_gamma={0}".format(gamma))
    jc.append("jc_gradcal_Tpm={0}".format(
        gradient_amplitude_T_per_m(maxgrad, fp, 1.0)))

    # ---- calibration-independent factorisation --------------------------
    # b = gamma^2 * G_abs^2 * S0, with G_abs = G_rel * FP * GradMax. These
    # three coefficients are CONSTANT across a whole gradient ramp, so one
    # of them plus your relative-gradient list rebuilds every b in the
    # experiment -- and the first two survive a change of gradient
    # calibration entirely. See bvalue_shape_factor().
    S0 = bvalue_shape_factor(P)
    b_per_G2 = gamma*gamma*S0
    b_per_Grel2 = b_per_G2 * (fp*maxgrad*GAUSS_PER_CM_TO_T_PER_M)**2
    comms.log("b-value factorisation (b = gamma^2 * G_abs^2 * S0, all three "
              "constant over a gradient ramp):")
    comms.log("   S0 (timing only, no gamma/no calibration) = {0:.9g} s^3"
              .format(S0))
    comms.log("   b / G_abs^2  = {0:.9g} s/m^2 per (T/m)^2   [needs only gamma]"
              .format(b_per_G2))
    comms.log("   b / G_rel^2  = {0:.9g} s/m^2               [uses GradMax={1} "
              "G/cm, FP={2}]".format(b_per_Grel2, maxgrad, fp))
    comms.log("   -> b = (b/G_rel^2) * G_rel^2, or (b/G_abs^2) * "
              "(G_rel * FP * GradMax_T_per_m)^2. NOTE the SQUARE: b/G_rel "
              "(first power) is not constant across a ramp.")
    jc.append("jc_b_shape_s3={0}".format(S0))
    jc.append("jc_b_per_G2={0}".format(b_per_G2))
    jc.append("jc_b_per_Grel2={0}".format(b_per_Grel2))
    jc.append("jc_gradmax_Gpcm={0}".format(maxgrad))
    jc.append("jc_fp_scaler={0}".format(fp))

    if int(P.ShowBTable) == 0:
        return jc
    try:
        rel_list = [float(x) for x in str(P.GradList).split(",") if x.strip() != ""]
    except ValueError:
        comms.log("WARNING: GradList='{0}' is not a comma-separated list of "
                  "numbers -- b-value table skipped.".format(P.GradList))
        return jc
    if not rel_list:
        return jc
    comms.log("---- b-value table (GradList, Probe='{0}', axis='{1}', "
              "Nucleus='{2}') ----".format(P.Probe, P.Axis, P.Nucleus))
    comms.log("  #    G_rel      G [T/m]      b [s/m^2]    b [s/mm^2]")
    b_list = []
    for i, rel in enumerate(rel_list):
        Gi = gradient_amplitude_T_per_m(maxgrad, fp, rel)
        bi = integrate_bvalue(bvalue_segments(P, Gi), gamma)
        b_list.append(bi)
        comms.log("  {0:<3d} {1:>8.4f}  {2:>10.5f}  {3:>13.6g} {4:>13.6g}"
                  .format(i+1, rel, Gi, bi, bi*1.0e-6))
    comms.log("---- end b-value table ----")
    jc.append("jc_gradlist={0}".format(",".join("{0}".format(x) for x in rel_list)))
    jc.append("jc_bvalue_list={0}".format(",".join("{0:.6g}".format(x) for x in b_list)))
    return jc

def mixing_delay(P):
    """The vendor DPFGSTE's 'd1' -- the remainder of each HALF of the DELTA
    diffusion period (there are two storage periods, P2->P3 and P4->P5),
    left over after the fixed events inside it. Kept ALGEBRAICALLY IDENTICAL
    to the vendor sequence's own expression; only the Parameter names
    changed (DiffDelta -> DELTA, GradientOnTime -> delta)."""
    return ((P.DELTA/2.0) - 2*(P.TXEnableTime + P.P90X)
            - P.delta - P.GradSpoil
            - 2*(P.PreGrad + 2*P.RampTime + P.GradSettle))


def eddy_delay(P):
    """The vendor DPFGSTE's 'd2' -- the remainder of the longitudinal
    eddy-current (LED) delay between P6 (third storage) and P7 (read)."""
    return (P.EddyCurrentDelay
            - (P.PreGrad + 2*P.RampTime + P.GradSpoil + P.GradSettle))


def bvalue_segments(P, G):
    """The EFFECTIVE gradient waveform of one scan, as a list of
    (duration_us, g_start, g_end) with g in T/m -- fed to integrate_bvalue().

    "Effective" = the coherence-pathway sign is already folded in, so that
    q(t) = gamma * integral(g_eff) dt is the real dephasing wavevector. The
    double stimulated echo winds q through FOUR encoding gradients, all of
    the same PHYSICAL sign (G1), whose EFFECTIVE signs alternate because
    each storage/restore 90 pair keeps only the cosine-modulated component
    and so discards the sign of the phase accumulated before it:

        g1 (+) : q  0  -> +k     encode, block 1
        [storage P2->P3, q frozen at +k for DELTA/2 -- this is block 1's
         diffusion time; the g2 homospoil acts on longitudinal
         magnetization and does not touch the stored helix]
        g3 (-) : q +k  ->  0     decode, block 1
        g4 (-) : q  0  -> -k     encode, block 2 (contiguous with g3, so
                                 same effective sign -- nothing between
                                 them flips it)
        [storage P4->P5, q frozen at -k for DELTA/2 -- block 2's diffusion
         time; g6 homospoil, again longitudinal]
        g7 (+) : q -k  ->  0     decode, block 2
        [storage P6->P7 = the LED eddy-current delay, with q ALREADY at 0
         -- which is exactly why it costs no extra diffusion weighting]

    k = gamma*G*(delta + RampTime): the trapezoid's area is G*(delta +
    RampTime), ONE ramp time, not two (equal-area rectangle equivalence,
    matching the vendor's own jc_smallDELTA = RampTime + GradientOnTime).

    Because the two blocks wind q to +k and -k respectively, a velocity
    (flow/convection) term accumulated in block 1 is cancelled in block 2
    while the diffusion attenuation adds -- the whole point of the double
    STE (J. Magn. Reson. 1997, 125, 372-375).

    Every duration below is the same expression run() actually emits, so the
    b-value tracks the real sequence timing rather than an idealised model.
    """
    R   = P.RampTime
    d   = P.delta
    PG  = P.PreGrad
    GS  = P.GradSettle
    GSp = P.GradSpoil
    TX  = P.TXEnableTime
    P90 = P.P90X
    d1  = mixing_delay(P)
    # SetBasePhase(3) + BlankingOn(1) + TXEnable + P90 + BlankingOff(1)
    pulse = 3.0 + 1.0 + TX + P90 + 1.0

    def trapezoid(amp):
        return [(1.0, 0.0, 0.0),        # SlewRate instruction
                (R,   0.0, amp),        # ramp up
                (d,   amp, amp),        # plateau (Stejskal-Tanner delta)
                (R,   amp, 0.0),        # ramp down
                (1.0, 0.0, 0.0)]        # explicit gradient-off instruction

    segments = []
    segments.append((PG - 2.0, 0.0, 0.0))
    segments.extend(trapezoid(G))                 # g1, encode block 1
    segments.append((GS - 5.0, 0.0, 0.0))
    segments.append((pulse, 0.0, 0.0))            # P2, store 1
    segments.append((PG - 2.0, 0.0, 0.0))
    segments.append((2.0 + 2*R + GSp, 0.0, 0.0))  # g2 homospoil (longitudinal)
    segments.append((GS - 5.0, 0.0, 0.0))
    segments.append((d1, 0.0, 0.0))               # DELTA/2 remainder
    segments.append((pulse, 0.0, 0.0))            # P3, restore 1
    segments.append((PG - 2.0, 0.0, 0.0))
    segments.extend(trapezoid(-G))                # g3, decode block 1
    segments.append((GS - 1.0, 0.0, 0.0))
    segments.append((PG - 1.0, 0.0, 0.0))
    segments.extend(trapezoid(-G))                # g4, encode block 2
    segments.append((GS - 5.0, 0.0, 0.0))
    segments.append((pulse, 0.0, 0.0))            # P4, store 2
    segments.append((PG - 2.0, 0.0, 0.0))
    segments.append((2.0 + 2*R + GSp, 0.0, 0.0))  # g6 homospoil (longitudinal)
    segments.append((GS - 5.0, 0.0, 0.0))
    segments.append((d1, 0.0, 0.0))               # DELTA/2 remainder
    segments.append((pulse, 0.0, 0.0))            # P5, restore 2
    segments.append((PG - 1.0, 0.0, 0.0))
    segments.extend(trapezoid(G))                 # g7, decode block 2
    # q is back to 0 from here on (LED storage + read pulse + acquisition),
    # so nothing after this point adds to b.
    return segments


def time_calculation(P):
    FilterFile = ChooseFilter(P.Filter)
    ReceiverFilter = Filter(FilterFile)

    # Term-by-term sum of every Delay()/duration argument in run()'s
    # (purely sequential) per-scan body, algebraically simplified: seven
    # 90s, seven gradient blocks (four encode/decode at G1 plus three
    # homospoils), seven PreGrad waits, seven GradSettle waits, 2*d1 and
    # d2 -- after substituting d1 and d2 the GradSpoil terms cancel
    # exactly and what survives is the expression below. NOTE: the vendor
    # DPFGSTE's own formula has 230 here; the correct constant is 231 (the
    # vendor's is 1 us short per scan -- see changelog).
    t_scanTime = ( 231 + P.RecycleDelay + 3*(P.TXEnableTime + P.P90X)
                  + 2*(P.PreGrad + 2*P.RampTime + P.delta + P.GradSettle)
                  + P.DELTA + P.EddyCurrentDelay + P.Dead1
                  + ReceiverFilter.dead_time
                  + P.ReceiverPoints*ReceiverFilter.dwell )

    t_acqTime = 10555 + t_scanTime * (P.NumScans + P.DS)

    return t_acqTime

def sequence_description():

    seq_desc = ("Double Pulsed-Field Gradient STimulated Echo (DPFGSTE) on "
                "the {X} channel -- two stimulated-echo diffusion "
                "blocks in series, so flow/convection phase cancels while "
                "diffusion attenuation adds, with a longitudinal "
                "eddy-current (LED) delay before the read pulse. Vendor "
                "DPFGSTE_X.py timing, in house style: probe-aware "
                "gradient calibration, Stejskal-Tanner delta/DELTA naming, "
                "and a b-value reported to the log and the saved JCAMP.")

    return seq_desc

def sequence_basic():

    basic = "NS,RD,NP,Filter,DELTA,delta,G1,GradAxis"

    return basic

@ParameterBlock
class Parameters:

    # Sequence name and basic parameter list for SpinFlow
    Sequence = Parameter("Sequence", "LWG_DPFGSTE_X", ParameterTypes.String, "Sequence Name")

    # General acquisition
    FrequencyBaseX = Parameter("SFX", 15.01, ParameterTypes.Double, "X Base Freq [MHz]")
    FrequencyOffsetX = Parameter("O1X", 0.0, ParameterTypes.Double, "X Freq Offset [Hz]")
    TxPPMX = Parameter("TxPPMX", 0.0, ParameterTypes.Double, "X TX Freq Offset [ppm]")

    ReceiverPoints = Parameter("NP", 32768, ParameterTypes.Int32, "Acquisition Points")
    ReceiverAttenuation = Parameter("RA", 30, ParameterTypes.Int32, "RX Attenuation [0&#8230;77dB]")
    Filter = Parameter("Filter", "10000", ParameterTypes.String, "f2 Spectral Window [Hz]")

    # Scans -- NumScans should be a MULTIPLE OF 8 to complete the phase
    # cycle below (PH4/PHRX are 8-step; PH1 is 4-step and PH3 2-step,
    # each repeating within one 8-step cycle).
    # DS keeps the vendor's default of 4, which is NOT a whole cycle: the
    # dummy scans are discarded, so this only means the steady state is
    # approached from a fixed partial cycle rather than a whole one. Raise
    # DS to 8 if you would rather it were a whole cycle.
    NumScans = Parameter("NS", 16, ParameterTypes.Int32, "Scans")
    NBlock = Parameter("NBlock", 4, ParameterTypes.Int32, "Blocks (NS/NBlock=integer)")
    DS = Parameter("DS", 4, ParameterTypes.Int32, "Dummy Scans")

    # Hardware and other standard delays
    Dead1 = Parameter("Dead1", 100.0, ParameterTypes.Double, "X Probe Ringdown Time [&#956;s]")
    TXEnableTime = Parameter("TXEnable", 20.0, ParameterTypes.Double, "TX Enable Time [&#956;s]")
    RecycleDelay = Parameter("RD", 2000000, ParameterTypes.Int32, "Relaxation Delay [s]", RD, min=100000, max=1800000000)

    # Hard pulses -- seven identical 90s: P1 excitation, P2/P4/P6 storage,
    # P3/P5 restore, P7 read. No 180 anywhere in this sequence.
    P90X = Parameter("P90X", 10.0, ParameterTypes.Double, "X 90&#176; Pulse Width [&#956;s]")
    TXAmplitudeX = Parameter("RFAX", 0.4, ParameterTypes.Double, "X TX Power [0.0&#8230;1.0]", RFA, min=0.0, max=1.0)

    # Gradients -- the SAME G1 amplitude and physical sign drives all FOUR
    # diffusion-encoding pulses (g1/g3/g4/g7; their EFFECTIVE signs alternate
    # via the storage pulses -- see bvalue_segments()). G2/G3/G4 are the
    # three homospoils, one per storage period. GradAxis picks which coil is
    # driven; FPX/FPY/FPZ are the only per-axis scalers, applied via
    # GradientMatrix() in run() -- and they scale the real output amplitude,
    # so they enter the b-value too.
    G1 = Parameter("G1", 1.00, ParameterTypes.Double, "Diffusion Gradient Strength [-1.0&#8230;1.0]")
    G2 = Parameter("G2", -0.13, ParameterTypes.Double, "Homospoil Strength, 1st storage [-1.0&#8230;1.0]")
    G3 = Parameter("G3", -0.17, ParameterTypes.Double, "Homospoil Strength, 2nd storage [-1.0&#8230;1.0]")
    G4 = Parameter("G4", -0.15, ParameterTypes.Double, "Homospoil Strength, LED storage [-1.0&#8230;1.0]")
    delta = Parameter("delta", 4000.0, ParameterTypes.Double, "Gradient Duration [&#956;s] -- Stejskal-Tanner &#948;, plateau only")
    RampTime = Parameter("D70", 500.0, ParameterTypes.Double, "Gradient Ramp Time [&#956;s]")
    GradSpoil = Parameter("D72", 1000.0, ParameterTypes.Double, "Gradient Homospoil Duration [&#956;s]")
    GradSettle = Parameter("D73", 500.0, ParameterTypes.Double, "Gradient Settling Duration [&#956;s]")
    PreGrad = Parameter("D75", 500.0, ParameterTypes.Double, "Pre-Gradient Time (pulse to gradient start) [&#956;s]")
    Axis = Parameter("GradAxis", "z", ParameterTypes.String, "Diffusion Gradient Axis")
    AxisList = Parameter("GradAxisList", "x,y,z,none", ParameterTypes.String, "Gradient Axes")
    XGradNorm = Parameter("FPX", 1.0, ParameterTypes.Double, "X Grad Scaler [0.0&#8230;1.0]")
    YGradNorm = Parameter("FPY", 1.0, ParameterTypes.Double, "Y Grad Scaler [0.0&#8230;1.0]")
    ZGradNorm = Parameter("FPZ", 1.0, ParameterTypes.Double, "Z Grad Scaler [0.0&#8230;1.0]")

    # Probe fitted to the magnet -- no automatic detection is possible, so
    # this must be set MANUALLY to match what is actually mounted. Recorded
    # in the saved .par/JCAMP, which is what lets the b-value be recomputed
    # downstream (tools/spinflow_bvalue_diff.py, leonmr/xpulse_imaging.py).
    Probe = Parameter("Probe", "HFX", ParameterTypes.String, "Probe fitted to magnet [HFX(default,calibrated)/Low Gamma(provisional)]")
    ProbeList = Parameter("ProbeList", "HFX,Low Gamma", ParameterTypes.String, "Probe Options")
    # Supply a fresh gradient calibration per run without editing this file.
    # Nonzero wins over MAXGRAD_TABLE for the selected GradAxis. Because it
    # is a Parameter, whatever you set is recorded in the .par/JCAMP, so the
    # b-value stays reproducible from the saved data alone.
    GradMaxOverride = Parameter("GradMax", 0.0, ParameterTypes.Double, "Max grad [G/cm] at |G1|=1.0, 0 = use Probe table")

    # Diffusion time -- TOTAL across BOTH storage periods (each gets
    # DELTA/2), matching the vendor's own D74 semantics. ARRAY G1 (not
    # DELTA) for a standard diffusion ramp.
    DELTA = Parameter("DELTA", 20000.0, ParameterTypes.Double,
                      "Diffusion Time &#916; [&#956;s], TOTAL -- split as &#916;/2 per storage period")
    EddyCurrentDelay = Parameter("D80", 5000.0, ParameterTypes.Double, "Eddy Current (LED) Delay T&#8337; [&#956;s]")

    # b-value reporting -- log/JCAMP only, no effect on timing or hardware.
    Nucleus = Parameter("Nucleus", "13C", ParameterTypes.String, "Observed nucleus -- for b-value only, no timing effect")
    NucleusList = Parameter("NucleusList", "1H,2H,6Li,7Li,11B,13C,15N,17O,19F,23Na,27Al,29Si,31P,35Cl,119Sn,133Cs,195Pt,207Pb", ParameterTypes.String, "Nucleus Options")
    GammaOverride = Parameter("GammaOverride", 0.0, ParameterTypes.Double, "&#947;/2&#960; [MHz/T], 0 = use Nucleus table")
    GradList = Parameter("GradList", "0.05,0.11,0.17,0.23,0.29,0.35,0.41,0.47,0.53,0.59,0.65,0.71,0.77,0.83,0.89,0.95", ParameterTypes.String, "Relative gradient list for the b-value table [comma-sep]")
    ShowBTable = Parameter("ShowBTable", 1, ParameterTypes.Int32, "Log full b-value table for GradList [0=Off,1=On(default)]")

    # Duty-cycle guard rails
    MaxRFDuty = Parameter("MaxRFDuty", 0.05, ParameterTypes.Double, "Max RF Duty Cycle Warning Threshold [0.0&#8230;1.0]")
    MaxGradDuty = Parameter("MaxGradDuty", 0.10, ParameterTypes.Double, "Max Gradient Duty Cycle Warning Threshold [0.0&#8230;1.0]")

    # Mains-lock trigger -- OFF by default (see mains_lock_trigger()).
    UseMainsLock = Parameter("UseMainsLock", 0, ParameterTypes.Int32, "Emit Mains-Lock Trigger Before Sequence [0=Off(default),1=On]")
    MainsLockChannel = Parameter("MainsLockChannel", 2, ParameterTypes.Int32, "Mains-Lock Trigger Channel [1-3, unconfirmed for X-Pulse]")

    # Phases -- vendor DPFGSTE_X.py's own cycle, unchanged in VALUE.
    # Attribute name == SpinFlow short code == PhasesManager dict key for all
    # eight (the vendor files' P1Phase = Parameter("PH1", ...) form is the
    # exact mismatch that caused a real KeyError on this instrument -- see
    # the pulse-programme-parameters skill, rule 8).
    PH1 = Parameter("PH1", "0,90,180,270", ParameterTypes.String, "P1 (excitation) RF-Pulse Phase")
    PH2 = Parameter("PH2", "0", ParameterTypes.String, "P2 (1st storage) RF-Pulse Phase")
    PH3 = Parameter("PH3", "180,270", ParameterTypes.String, "P3 (1st restore) RF-Pulse Phase")
    PH4 = Parameter("PH4", "180,180,180,180,0,0,0,0", ParameterTypes.String, "P4 (2nd storage) RF-Pulse Phase")
    PH5 = Parameter("PH5", "0", ParameterTypes.String, "P5 (2nd restore) RF-Pulse Phase")
    PH6 = Parameter("PH6", "0", ParameterTypes.String, "P6 (LED storage) RF-Pulse Phase")
    PH7 = Parameter("PH7", "0", ParameterTypes.String, "P7 (read) RF-Pulse Phase")
    PHRX = Parameter("PHRX", "0,0,180,180,180,180,0,0", ParameterTypes.String, "Acquisition Phase")

def run(comms):

    P = Parameters

    FPX = P.XGradNorm
    FPY = P.YGradNorm
    FPZ = P.ZGradNorm

    Matrix = np.array([
            [FPX,  0,  0],
            [0  ,FPY,  0],
            [0  ,  0,FPZ]
            ])

    Frequency = P.FrequencyBaseX + (P.FrequencyOffsetX*1.0e-6) + (P.FrequencyBaseX*P.TxPPMX*1.0e-6)
    Channel2SetFrequency(10.0, Frequency)
    Channel2RestartSynth(10.0)

    ReceiverFilter = Filter(P.Filter)
    ReceiverTime = ReceiverFilter.dwell*(P.ReceiverPoints+0.99)

    Phases = PhasesManager(P)
    Phases.Reset()

    points = P.ReceiverPoints
    times = np.arange(0, points*ReceiverFilter.dwell, ReceiverFilter.dwell) / 1.0e6

    GradientSlewRateFn, GradientFn = get_gradient_functions(P.Axis)

    # ---- Probe / gyromagnetic ratio / b-value report ----------------------
    # report_bvalue() logs the probe's gradient calibration, this scan's
    # b-value, and the whole b-value table for GradList, and hands back the
    # JCAMP strings that carry those numbers into the saved dataset.
    gamma = gamma_rad_per_s_per_T(P, comms)
    bvalue_jcamp = report_bvalue(P, comms, gamma)

    # Absolute gradient amplitude [T/m] for the vendor's jc_grad key --
    # looked up silently here, since report_bvalue() has already logged it.
    if P.GradMaxOverride > 0.0 and str(P.Axis).strip().lower() != 'none':
        maxgrad = P.GradMaxOverride
    else:
        maxgrad = MAXGRAD_TABLE.get(normalise_key(P.Probe), {}).get(str(P.Axis).strip().lower())
    grad = 0.0 if maxgrad is None else gradient_amplitude_T_per_m(maxgrad, axis_fp_scaler(P), P.G1)

    def jcamp_meta():
        global BLP
        if BLP == 1:
            me_mod = 4
            jc_blp = -2*int(round((((ReceiverFilter.dead_time/2.0)+P.Dead1)/ReceiverFilter.dwell),0))
        else:
            me_mod = 0
            jc_blp = 0
        jc_data = ["jc_sf={0}".format(Frequency),
                "jc_sf_txppm={0}".format(Frequency+(P.FrequencyBaseX*P.TxPPMX*1.0e-6)),
                "ME_mod={0}".format(me_mod),
                "jc_blp={0}".format(jc_blp),
                "jc_grad={0}".format(grad),
                "jc_axis={0}".format(P.Axis),
                "jc_bigDELTA={0}".format(P.DELTA / 1000),
                "jc_smallDELTA={0}".format((P.delta + P.RampTime) / 1000),
                "jc_ns=1"
                    ]
        jc_data.extend(bvalue_jcamp)
        return jc_data

    comms.log("Sequence Time: {0}".format(time_calculation(Parameters)))

    template = "JCamp1D_Diffusion"      # JCAMP-DX Template = "JCamp1D_Diffusion.txt"
    seqAcqu = CallBack1D(Parameters, comms, ReceiverFilter, jcamp_meta(), template)

    def RecvCallback(acqData, scan):
        seqAcqu.process_data(scan, acqData)

    d1 = mixing_delay(P)
    d2 = eddy_delay(P)

    with sequential:

        # Mains-lock trigger, if enabled -- must come before the first pulse
        # event (manual 3.7.19). OFF by default; see mains_lock_trigger().
        mains_lock_trigger(P, comms)

        Transmit2SelectPort(1.0, 1)
        Transmit2LPEnable(1.0, 0)
        Delay(10000)    # This is required for changing the relay state from tune mode

        Transmit2SetScale(5.0, P.TXAmplitudeX)

        Receiver2Preamp(128, P.ReceiverAttenuation)
        Receiver2Filter(200, ReceiverFilter)

        # Gradient Setup
        GradientMatrix(200, Matrix.T)

        Phases.Reset()

    # ---- Duty-cycle warning ----------------------------------------------
    estimate_duty_cycles(P, 7*P.P90X,
                         4*(2*P.RampTime + P.delta) + 3*(2*P.RampTime + P.GradSpoil), comms)

    for seqScans in range(P.NumScans+P.DS):
        with sequential_main(seqScans, P.NumScans+P.DS):

            if seqScans == P.DS:
                Phases.Reset()

            ph = Phases.Incd()

            Channel2SetFrequency(10.0, Frequency)
            Channel2SetBasePhase(10.0, 0)
            Channel2RestartSynth(10.0)

            Receiver2FilterFlush(200, ReceiverFilter)

            # RD
            Delay(P.RecycleDelay-5.0e5-4.0)

            # P1 -- excitation 90
            Channel2SetBasePhase(3.0, ph["PH1"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            Delay(P.PreGrad-2.0)

            # g1 -- diffusion ENCODE gradient, block 1 [G1, duration delta]
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.delta, P.G1)

            Delay(P.GradSettle-5.0)

            # P2 -- 1st storage 90
            Channel2SetBasePhase(3.0, ph["PH2"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            Delay(P.PreGrad-2.0)

            # g2 -- homospoil during the 1st storage period [G2]
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.GradSpoil, P.G2)

            Delay(P.GradSettle-5.0)

            # d1 -- 1st DELTA/2 remainder
            safe_delay(d1, "1st DELTA/2 remainder (d1)", comms)

            # P3 -- 1st restore 90
            Channel2SetBasePhase(3.0, ph["PH3"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            Delay(P.PreGrad-2.0)

            # g3 -- diffusion DECODE gradient, block 1 [G1]
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.delta, P.G1)

            Delay(P.GradSettle-1.0)

            Delay(P.PreGrad-1.0)

            # g4 -- diffusion ENCODE gradient, block 2 [G1]
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.delta, P.G1)

            Delay(P.GradSettle-5.0)

            # P4 -- 2nd storage 90
            Channel2SetBasePhase(3.0, ph["PH4"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            Delay(P.PreGrad-2.0)

            # g6 -- homospoil during the 2nd storage period [G3]
            # (vendor's own gradient numbering; there is no g5)
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.GradSpoil, P.G3)

            Delay(P.GradSettle-5.0)

            # d1 -- 2nd DELTA/2 remainder
            safe_delay(d1, "2nd DELTA/2 remainder (d1)", comms)

            # P5 -- 2nd restore 90
            Channel2SetBasePhase(3.0, ph["PH5"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            Delay(P.PreGrad-1.0)

            # g7 -- diffusion DECODE gradient, block 2 [G1]
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.delta, P.G1)

            Delay(P.GradSettle-5.0)

            # P6 -- LED storage 90 (parks magnetization along z while the
            # gradient eddy currents decay -- q is already 0 here, so this
            # period costs no extra diffusion weighting)
            Channel2SetBasePhase(3.0, ph["PH6"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            Delay(P.PreGrad-2.0)

            # g8 -- homospoil during the LED delay [G4]
            apply_gradient(GradientSlewRateFn, GradientFn, P.RampTime, P.GradSpoil, P.G4)

            Delay(P.GradSettle-5.0)

            # d2 -- LED (eddy-current) delay remainder
            safe_delay(d2, "LED eddy-current delay remainder (d2)", comms)

            # P7 -- read 90
            Channel2SetBasePhase(3.0, ph["PH7"])
            Transmit2BlankingOn(1.0)
            Delay(P.TXEnableTime)
            Transmit2(P.P90X)
            Transmit2BlankingOff(1.0)

            # ACQU
            Channel2SetBasePhase(1.0, 0)
            Receiver2Phase(1.0, ph["PHRX"])
            Delay(P.Dead1-3.0)
            Delay(ReceiverFilter.dead_time)
            Receiver2(P.ReceiverPoints*ReceiverFilter.dwell, P.ReceiverPoints)

        Delay(5.0e5)
        WriteToHardware(seqScans, P.NumScans+P.DS)
        TX1.setup_receive(points, 32, scan=seqScans, circular=True)
        start(seqScans, P.NumScans+P.DS)
        TX1.wait_for_data(seqScans, P.NumScans+P.DS, RecvCallback)
        comms.log("Loop = %s exec time = %s compile %s" % (seqScans, get_single_scan_execution_time(), Get_Compilation_Time()))

# -----------------------------------------------------------------------------
# Changes/Modifications (Initials - Date - Description):
#
# 1. Claude - 15/09/26 - Initial version: house-style port of the vendor
#    DPFGSTE_X.py you added to sequences/vendor/, at your request ("take
#    those vendor sequences as a starting point and update them with the
#    house style"). EVERY Delay()/duration expression is algebraically
#    identical to the vendor's -- re-derived term by term and re-checked
#    against the vendor file line by line, since the whole reason for
#    starting from the vendor sequence is that its timing works on this
#    instrument. What changed: phase-cycle Parameters renamed so attribute
#    name == SpinFlow short code == PhasesManager key (PH1/.../PHRX) and
#    routed through a single PhasesManager(P) -- this is the "KeyError:
#    'PH1'" class of bug the vendor form causes the moment it meets the
#    house phase convention, see the pulse-programme-parameters skill rules
#    7-8; D74/D71 renamed to DELTA/delta; Gx_max/Gy_max/Gz_max replaced by
#    Probe + MAXGRAD_TABLE + FPX/FPY/FPZ; safe_delay() on every computed
#    delay; apply_gradient() replacing the copy-pasted per-axis gradient
#    blocks (identical instruction stream, including the GradAxis='none'
#    pure-Delay branch); mains-lock trigger and duty-cycle guard rails
#    added.
# 2. Claude - 15/09/26 - b-VALUE REPORTING added (the second half of your
#    request). integrate_bvalue() computes b = integral q(t)^2 dt over the
#    effective gradient waveform that bvalue_segments() builds from this
#    file's OWN timing expressions, using 3-point Gauss-Legendre per
#    segment -- which is EXACT here, not approximate, because q(t) is
#    quadratic within a segment so q^2 is quartic (Gauss-Legendre with 3
#    nodes is exact to degree 5). Verification actually run, not asserted:
#      (a) rectangular pair vs b = (gamma*G*delta)^2*(DELTA - delta/3),
#          2000 randomised (G, delta, DELTA) cases -> max relative error
#          5.6e-16;
#      (b) trapezoidal pair vs the Price-Kuchel closed form
#          b = gamma^2*G^2*[de^2*(DELTA - de/3) + ramp^3/30 - de*ramp^2/6]
#          with de = delta + ramp (HALF-HEIGHT duration -- this convention
#          is the whole trick; using the bare plateau instead is off by
#          ~20% at RampTime=500 us) and DELTA leading-edge to leading-edge,
#          3000 randomised (G, delta, ramp, DELTA) cases -> max relative
#          error 1.2e-15;
#      (c) this file's OWN bvalue_segments() output vs that same closed
#          form, over RampTime 50/200/500 us x delta 1/4 ms x DELTA 20/60
#          ms -> relative error <= 4.2e-16 everywhere (check (c) was run against LWG_PFGSTE_H's single-block segments; this file's four-gradient waveform has no simple closed form, so it rests on checks (a),(b),(d) plus the wavevector walk-through in bvalue_segments());
#      (d) exact integrator vs a brute-force 0.05 us-step numerical
#          integration of the identical waveform -> agreement to 3.4e-11,
#          i.e. limited by the brute-force method, not by this one.
#    Reported three ways: the current G1's b to the log, the whole GradList
#    table to the log, and both into the saved JCAMP metadata.
# 3. Claude - 15/09/26 - NOTE ON WHERE THE .diff FILE GOES. You asked for
#    it in the experiment's parameter folder. run() executes on the
#    spectrometer console (/home/fire/PythonMonitor, sequence copied to
#    /tmp/tmpXXXXXX.py -- visible in docs/error-logs/SpinFlow_error_logs),
#    not on the Windows PC that holds C:\Users\Public\Documents\SpinFlow\
#    Parameters\<serial>\User\<nucleus>\, so a file written from here
#    would land on the wrong machine. tools/spinflow_bvalue_diff.py does it
#    from the SpinFlow PC instead, reading the saved .par and writing
#    <Sequence>.diff beside it with the same maths.
# 4. Claude - 15/09/26 - GRADIENT CALIBRATION DISCREPANCY FLAGGED, not
#    silently resolved: the vendor's Gx_max/Gy_max/Gz_max defaults (0.01,
#    0.01, 0.05 T/m) are ~12x smaller than this repo's own measured
#    MAXGRAD_TABLE (HFX x/y/z = 11.879/11.978/57.915 G/cm = 0.1188/0.1198/
#    0.5792 T/m), i.e. ~140x in b. MAXGRAD_TABLE is used here (it is this
#    repo's own measurement, on this magnet, and is already what
#    leonmr/xpulse_imaging.py uses), but absolute b-values from this file
#    are UNVALIDATED until checked against a known-D standard.
# 5. Claude - 15/09/26 - b FACTORISED into calibration-independent parts,
#    at your request ("is there another way to get b ... then I can use the
#    relative gradients file and the known max G per probe"). Two answers,
#    one correction:
#      (a) CORRECTION: b/G_rel is NOT the right normalisation -- b goes as
#          G^2, so b/G_rel still varies across a ramp by the ramp's own
#          dynamic range (measured: 4.50e8 -> 9.00e9 over G_rel 0.05..1.0,
#          a factor of 20). b/G_rel^2 is exactly constant (identical to 10
#          significant figures across the same range).
#      (b) The sequence now logs and records three constants --
#          S0 = jc_b_shape_s3 (timing only), jc_b_per_G2 = gamma^2*S0
#          (calibration-free), and jc_b_per_Grel2 (today's calibration
#          folded in) -- so b can be rebuilt from the relative-gradient
#          list alone, or rebuilt under a different GradMax later.
#          Verified: b(G)/G^2 constant to 3.5e-16 over G_rel 0.05..1.0;
#          (b/G_rel^2)*G_rel^2 reproduces every row of the logged table;
#          and a GradMax override changes b/G_rel^2 by exactly the
#          calibration ratio squared while leaving S0 and b/G_abs^2
#          bit-identical.
#      (c) It is NOT possible to write this back into a real SpinFlow
#          Parameter from run(). The whole comms API is log/logs/send_data/
#          lock_status/error -- there is no parameter write-back, and no
#          precedent for one anywhere in this repo. The JCAMP metadata is
#          the only channel out, and it is per-dataset, so it is already
#          "per slice". Every INPUT needed to recompute b is a Parameter
#          and so is in the .par already, which is what lets
#          tools/spinflow_bvalue_diff.py rebuild the table off a saved .par.
#    Also added: a GradMax Parameter (0 = use MAXGRAD_TABLE) to supply a
#    fresh calibration per run, honoured by both the sequence and the tool
#    (--gmax); and tools/spinflow_bvalue_diff.py --solve-gmax <slope>
#    --known-d <D>, which inverts a measured decay on a known-D sample to
#    the gradient calibration it implies (round-trip tested in both
#    directions against the vendor and MAXGRAD_TABLE numbers).
# 6. Claude - 15/09/26 - Mechanical H->X port of LWG_DPFGSTE_H.py
#    (Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/
#    Receiver2/TX1, SF default -> 15.01 MHz placeholder, Nucleus
#    default -> 13C placeholder, NucleusList widened to the full
#    gamma table). Keeps the vendor _X CallBack1D verbatim (it
#    accumulates across scans and honours NBlock, which the _H one
#    does not) rather than reusing the _H callback the way
#    LWG_PGSTE_X.py does -- the vendor pair genuinely differs here
#    and the vendor behaviour is what you reported works.
# 7. Claude - 15/09/26 - LOW GAMMA CALIBRATION SURFACED TO THE OPERATOR.
#    MAXGRAD_TABLE's LOWGAMMA row has held real numbers since 18-19/08/2026,
#    but the Probe Parameter's description still read "not yet calibrated"
#    and nothing warned at run time -- so selecting Low Gamma silently used a
#    provisional z (PGSE-only) and x/y (PGSTE) calibration for quantitative
#    work while SpinFlow's panel said it did not exist, and the caveat lived
#    only in a source comment nobody at the instrument would read. Added
#    PROVISIONAL_CALIBRATION (keyed like MAXGRAD_TABLE, value = the caveat)
#    and a report_probe_gradient() WARNING that prints it to the SpinFlow
#    console when such a probe is selected; Probe's description now reads
#    "Low Gamma(provisional)". No calibration NUMBERS changed, no timing
#    changed, and HFX is unaffected (verified: the warning fires for Low
#    Gamma and not for HFX, across all 22 gradient sequences under the mock
#    harness).
#
# -----------------------------------------------------------------------------
