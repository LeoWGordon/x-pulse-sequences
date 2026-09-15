# Validation note: LWG_DPFGSTE_X

- **File:** `sequences/diffusion/LWG_DPFGSTE_X.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample)
- **Ported from:** `sequences/vendor/DPFGSTE_X.py`

## What was checked
- Syntax (`ast.parse`) and CRLF line endings confirmed on all four files.
- Full mock-hardware-harness run (`tools/pp_visualizer.py`) through `run()`
  at default Parameters, and with `GradAxis` set to each of `x`, `y`, `z`
  and `none` -- all four axis settings produce an **identical** total scan
  time (2080754.0 us at the harness's placeholder Filter), confirming
  `apply_gradient()`'s `none` branch is timing-neutral.
- `time_calculation()` verified against the harness trace rather than
  asserted: with `NumScans=1, DS=0` the traced one-scan total minus
  `t_acqTime - 10555` comes out at **exactly 10555.0 us** (the `with
  sequential:` prologue), i.e. the per-scan formula is exact to the
  microsecond. NOTE: the vendor `DPFGSTE_X.py`'s own formula uses 230 here; the correct constant is 231, confirmed by the prologue landing on exactly 10555 us with 231 and 10556 us with 230. The vendor formula is 1 us short per scan.
- Phase-cycle Parameters go through `PhasesManager(P)` with attribute name
  == SpinFlow short code == `Phases.Incd()` key, so the `KeyError: 'PH1'`
  that the vendor's `P1Phase = Parameter("PH1", ...)` form produces on this
  instrument (see `docs/error-logs/SpinFlow_error_logs`, 14/08/2026) cannot
  recur. The hardened mock `PhasesManager` exercises this.
- A deliberately too-short `DELTA` (100 us) raises the expected
  `safe_delay()` FATAL TIMING ERROR (`ValueError`) naming the delay,
  instead of emitting a negative `Delay()`.
- Bad inputs degrade to a warning, never a crash: unrecognised `Probe`,
  unrecognised `Nucleus`, and a non-numeric `GradList` each log a specific
  warning and skip only the b-value, leaving the sequence runnable.
  `ShowBTable=0` suppresses the table.
- `GammaOverride` correctly overrides the `Nucleus` table lookup, and the
  implied-B0 cross-check (`SF / (gamma/2pi)`) fires its warning when the
  two are inconsistent.

## b-value verification

`integrate_bvalue()` computes `b = integral q(t)^2 dt` over the effective
gradient waveform that `bvalue_segments()` builds from the sequence's own
timing expressions, using 3-point Gauss-Legendre per segment (exact here:
`q(t)` is quadratic within a segment, so `q^2` is quartic). Checks actually
run, not asserted:

| Check | Result |
|---|---|
| Rectangular pair vs `b = (gamma*G*delta)^2 (DELTA - delta/3)`, 2000 random cases | max rel. error **5.6e-16** |
| Trapezoidal pair vs Price-Kuchel `b = gamma^2 G^2 [de^2(DELTA - de/3) + ramp^3/30 - de*ramp^2/6]`, `de = delta + ramp`, 3000 random cases | max rel. error **1.2e-15** |
| Exact integrator vs brute-force 0.05 us-step numerical integration of the identical waveform | agreement to **3.4e-11** (limited by the brute-force method) |
| This file's four-gradient waveform | no simple closed form exists; rests on the three checks above plus the wavevector walk-through documented in `bvalue_segments()` |

| `tools/spinflow_bvalue_diff.py` output vs the sequence's own logged table, non-default parameter set | agreement to **1.3e-6**, limited only by the 6-significant-figure log formatting |

### Factorisation (calibration-independent b)

`b = gamma^2 * G_abs^2 * S0` with `G_abs = G_rel * FP * GradMax`. The sequence
logs and records all three constants: `S0` (`jc_b_shape_s3`, timing only),
`b/G_abs^2` (`jc_b_per_G2`, needs only gamma), and `b/G_rel^2`
(`jc_b_per_Grel2`, today's calibration folded in). Checks run:

| Check | Result |
|---|---|
| `b(G)/G^2` constant over `G_rel` = 0.05..1.0 | constant to **3.5e-16** |
| `b = (gamma^2 S0) * G_abs^2` vs directly integrated b | **3.5e-16** |
| `(b/G_rel^2) * G_rel_i^2` rebuilds every row of the logged table | **<= 3.9e-6** (the 6-s.f. log formatting floor) |
| `GradMax` override (57.915 -> 5.0 G/cm) | `b/G_rel^2` scales by exactly `(5/57.915)^2 = 0.0075`; `S0` and `b/G_abs^2` bit-identical |
| `--solve-gmax` round trip, both directions | recovers 57.9150 and 4.9999 G/cm from the slopes those calibrations would produce |

`b/G_rel` (**first** power) is *not* a usable normalisation -- it varies by the
ramp's own 20x dynamic range (4.50e8 -> 9.00e9 over `G_rel` 0.05..1.0). The
invariant is `b/G_rel` **squared**.

**Watch the `delta` convention** when comparing against a textbook: the
quantity that belongs in those closed forms is the half-height duration
`delta_eff = delta + RampTime` (the trapezoid's area is `G*(delta +
RampTime)` -- one ramp time, not two), and their `DELTA` means leading edge
to leading edge, which is *not* this file's `DELTA` Parameter. Using the
bare plateau instead understates b by ~20% at the default
`RampTime=500 us, delta=4000 us`.

## Sequence-specific notes

- Timing is a verbatim port of vendor `DPFGSTE_X.py`, including its small
  -1/-2/-3/-5 us trim constants.
- The four diffusion gradients (`g1`, `g3`, `g4`, `g7`) all share one
  physical sign and amplitude (`G1`); their *effective* signs alternate
  `+,-,-,+` through the storage pulses, which is what cancels flow phase
  while diffusion attenuation adds. `bvalue_segments()` documents the full
  wavevector walk. Do not hand-flip any of them.
- The LED tail (`P6` - `g8` homospoil - `d2` - `P7`) runs with `q` already
  back to 0, so it costs no extra diffusion weighting -- worth confirming
  experimentally by varying `EddyCurrentDelay` and checking b-independence
  of the decay.
- `DELTA` is the TOTAL diffusion time, split `DELTA/2` per storage period
  (the vendor's own `D74` semantics), unlike `LWG_PFGSTE_X.py`.

## Results / evidence
- Not yet run on real hardware. Everything above is from the mock
  harness and from offline numerical checks.

## Known limitations or open questions
- **The `_X` Parameter block follows the vendor's `_X` file, not the `_H`
  one.** The observe-channel short codes are `SFX`/`O1X`/`TxPPMX`/`P90X`/
  `RFAX` (not `SF`/`O1`/`TxPPM`/`P90`/`RFA0`), `NBlock` is present because
  the vendor `_X` `CallBack1D` genuinely honours it, and `Filter`/`NS`/`P90X`
  keep the vendor `_X` defaults. This deliberately diverges from
  `LWG_PGSTE_X.py` and `LWG_1D-Image-Echo_X.py`, which normalise those codes
  to the `_H` spellings. **Open question:** whether SpinFlow actually cares.
  It does decide channel/nucleus validity from the parameter set (see the
  `Sequence invalid for selection: ... Observed: 1H` errors in
  `docs/error-logs/SpinFlow_error_logs`), so the vendor spelling was taken as
  the safer default -- but if the existing house `_X` sequences load and run
  fine on the X channel with `SF`, this can be normalised.
- **X-channel RF parameters are not calibrated.** `SFX` keeps the vendor
  default (15.01 MHz) and `Nucleus` defaults to `13C`, which are consistent
  with each other on this magnet (10.7084 MHz/T x ~1.40 T). `P90X`, `RFAX`
  and `RA` are still vendor defaults and are almost certainly wrong for
  whatever is actually on the X channel. Recalibrate all of them -- and set
  `Nucleus`, which b depends on as `gamma^2` -- before running on a sample.
- **Not yet run on hardware.** Status is `draft`.
- **Absolute b-values are unvalidated.** The vendor sequences' own
  `Gx_max`/`Gy_max`/`Gz_max` defaults (0.01, 0.01, 0.05 T/m) are ~12x
  smaller than this repo's measured `MAXGRAD_TABLE` (HFX x/y/z =
  11.879/11.978/57.915 G/cm = 0.1188/0.1198/0.5792 T/m). Since `b` scales
  as `G^2`, the two disagree by ~140x in b. `MAXGRAD_TABLE` is used here,
  because it is this repo's own measurement on this magnet and is already
  what `leonmr/xpulse_imaging.py` uses for Hz->mm conversion -- but this
  needs settling against a known-D standard (e.g. H2O/D2O at a known
  temperature) before any absolute b or D from this sequence is trusted.
  Relative b-values within one gradient ramp are unaffected by this, and
  because the b-value is factorised (above), absolute b can be rebuilt from
  already-acquired data under a corrected calibration without re-running
  anything -- `jc_b_shape_s3` and `jc_b_per_G2` do not depend on `GradMax`.
  `tools/spinflow_bvalue_diff.py --solve-gmax <slope> --known-d <D>` inverts a
  ramp on a known-D sample to the calibration it implies, which is the
  intended way to close this out.
- **`jc_bvalue` and friends in the JCAMP metadata are unverified.** Whether
  the vendor's `JCamp1D_Diffusion.txt` template silently ignores keys it
  does not know has not been confirmed. If a dataset fails to save or looks
  malformed, delete the `jc_data.extend(bvalue_jcamp)` line in
  `jcamp_meta()`; the logged table and the `.diff` file are unaffected.
- **b cannot be written back into a SpinFlow Parameter.** The whole `comms`
  API is `log`/`logs`/`send_data`/`lock_status`/`error` -- there is no
  parameter write-back, and no precedent for one in the repo. b reaches the
  saved data through the JCAMP metadata (which is per-dataset, so already
  per-slice) and through the log. Every *input* needed to recompute b is a
  Parameter and so is in the `.par`, which is what lets the tool rebuild the
  table offline.
- **`GradList` is operator-entered.** A pulse programme cannot see the
  array SpinFlow is sweeping, so `GradList` must be set to match the
  gradient array by hand for the logged table (and the `.diff` file) to
  correspond to what is actually being run. The per-experiment `G1` and its
  b-value are always correct regardless.
- **Renaming `D74`/`D71` to `DELTA`/`delta` orphans existing `.par` files**
  for the vendor sequence -- SpinFlow will fall back to the defaults in
  this file the first time it is loaded. Check the parameter panel before
  the first real run.
- Duty-cycle thresholds (`MaxRFDuty`, `MaxGradDuty`) remain conservative
  placeholders, not vendor-confirmed amplifier ratings.
- `MainsLockChannel` remains unconfirmed for the X-Pulse; the trigger is
  off by default.

- Needs on the instrument: confirm the echo appears where expected and
  phases cleanly; confirm a low-`G1` point matches an unweighted FID in
  integral; confirm `ln(S/S0)` vs `-b` is linear on a known-D standard and
  that the fitted D matches the literature value -- that single measurement
  settles the gradient-calibration question above.
