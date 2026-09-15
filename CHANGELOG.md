# Changelog

All notable changes to this repo are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- `sequences/diffusion/LWG_PFGSTE_H.py` / `_X.py` and
  `sequences/diffusion/LWG_DPFGSTE_H.py` / `_X.py` (draft): house-style ports
  of the vendor `PFGSTE`/`DPFGSTE` sequences Leo added to `sequences/vendor/`
  (which they found work better on this instrument than `LWG_PGSTE_*.py`).
  **Timing is a verbatim port** -- every `Delay()`/duration expression is
  algebraically identical to the vendor's, trim constants included -- with
  the house layer on top: phase-cycle Parameters renamed so attribute name ==
  SpinFlow short code == `PhasesManager` key (the vendor's
  `P1Phase = Parameter("PH1", ...)` form is what produced the real
  `KeyError: 'PH1'` on this instrument on 14/08/2026), `D74`/`D71` renamed to
  `DELTA`/`delta`, `Gx_max`/`Gy_max`/`Gz_max` replaced by the `Probe`
  selector + `MAXGRAD_TABLE` + `FPX`/`FPY`/`FPZ`, `safe_delay()` on every
  computed delay, `apply_gradient()` collapsing the vendor's copy-pasted
  per-axis gradient blocks (7 copies in DPFGSTE) into one timing-identical
  helper, plus the mains-lock trigger and duty-cycle guard rails. `DPFGSTE`'s
  `time_calculation()` constant corrected from the vendor's 230 to 231 (the
  vendor formula is 1 us short per scan; confirmed against the harness).
- **b-value reporting** in all four: `integrate_bvalue()` computes
  `b = integral q(t)^2 dt` over the effective gradient waveform that
  `bvalue_segments()` builds from each sequence's *own* timing expressions --
  trapezoidal ramps, RF pulse widths and frozen-q storage periods included,
  no idealised-rectangle assumption -- via 3-point Gauss-Legendre, which is
  exact here rather than approximate. Verified to machine precision against
  both the rectangular Stejskal-Tanner form and the Price-Kuchel trapezoid
  form (see the validation notes for the numbers). Reported three ways: the
  current `G1`'s b to the SpinFlow log, the whole `GradList` table to the
  log, and both into the saved JCAMP metadata. New `Nucleus`/`NucleusList`/
  `GammaOverride` Parameters supply gamma, with an implied-B0 cross-check
  that warns when `Nucleus` and `SF` disagree.
- **b factorised into calibration-independent parts**, so an absolute
  b-value is not hostage to the gradient table: `b = gamma^2 * G_abs^2 * S0`
  with `G_abs = G_rel * FP * GradMax`. All four sequences log and record `S0`
  (`jc_b_shape_s3`, pure timing -- no gamma, no probe, no calibration),
  `b/G_abs^2` (`jc_b_per_G2`, gamma only) and `b/G_rel^2`
  (`jc_b_per_Grel2`, today's calibration folded in). All three are constant
  across a gradient ramp, so one of them plus the relative-gradient list
  rebuilds every b in the experiment -- and the first two survive a
  recalibration, letting already-acquired data be re-costed without
  re-running it. Note the power: `b/G_rel` (first power) is *not* constant,
  varying by the ramp's own 20x range; the invariant is `b/G_rel` squared.
  New `GradMax` Parameter (0 = use `MAXGRAD_TABLE`) supplies a calibration
  per run without editing the file, and is recorded in the `.par`.
- `tools/spinflow_bvalue_diff.py --solve-gmax <slope> --known-d <D>`: inverts
  a measured diffusion decay on a known-D sample to the gradient calibration
  it implies (`slope` of `ln(S/S0)` vs `G_rel^2`, so no calibration is assumed
  anywhere). Round-trip tested in both directions against the vendor and
  `MAXGRAD_TABLE` numbers. This is the intended way to settle the open
  question below. `--gmax` re-costs an existing experiment under a different
  calibration.
- `tools/spinflow_bvalue_diff.py`: writes the b-value table as a `.diff`
  file into a SpinFlow experiment's own parameter folder. This has to be a
  host-side tool rather than something the sequence does itself -- `run()`
  executes on the spectrometer console (`/home/fire/PythonMonitor`, sequence
  copied to `/tmp/tmpXXXXXX.py`, visible in
  `docs/error-logs/SpinFlow_error_logs`), not on the Windows PC that owns
  `C:\Users\Public\Documents\SpinFlow\Parameters\...`. It imports
  `bvalue_segments()`/`integrate_bvalue()`/`MAXGRAD_TABLE`/the gamma table
  straight out of the pulse programme under a stub `firebird`, so there is
  no duplicated maths to keep in sync, and overlays a saved `.par` on the
  sequence's own defaults (both XML and `key=value` `.par` shapes are
  auto-detected).
- `docs/validation-notes/LWG_PFGSTE_H.md`, `_X.md`, `LWG_DPFGSTE_H.md`,
  `_X.md`.
- Initial repo scaffold: README, LICENSE, validation status legend, folder
  structure for diffusion / relaxation / imaging / calibration sequences.
- `sequences/imaging/LWG_CPMG-Image-Echo_H.py` (draft): multi-echo (CPMG)
  1D image-echo sequence, for separating T2 from diffusion contributions
  near a liquid-liquid interface. v2.0: applies NECH refocusing pulses per
  scan (timing/phase-cycle copied verbatim from the vendor's own
  `CPMG_H.py`) and images only the last one (diffusion-compensated
  gradient geometry copied from `LWG_1D-Image-Echo-GradAfterRefocus_H.py`)
  -- run as a quasi-2D sweep of NECH to build a T2(x) curve with short,
  fixed echo spacing. Supersedes a v1.0 draft that instead tried to
  acquire every echo from a single excitation (no precedent for that
  anywhere in this repo; not carried forward). See the file's own header
  and `docs/validation-notes/LWG_CPMG-Image-Echo_H.md` for what's still
  unverified -- not yet run on hardware.
- `sequences/diffusion/LWG_PGSE-Image_H.py` / `_X.py` (draft): a "PFG spin
  echo with an imaging gradient during acquisition" -- merges
  `LWG_PGSE-WET_H.py`'s plain two-pulse Stejskal-Tanner PGSE skeleton with
  `LWG_1D-Image-UTE_H.py`'s "Hold" gradient-during-acquisition pattern,
  giving a diffusion-weighted 1D profile where the b-value (DiffGrad/
  GradientOnTime/TAU) and the spatial resolution/FOV (ReadGrad/NP/Filter)
  are set by fully independent gradient channels/amplitudes -- unlike
  `diffprof.lwg`/`LWG_diffprof_2_H.py`, where one gradient amplitude
  serves both roles. Cross-checked against
  `LWG_1D-Image-Echo-GradAfterRefocus_H.py` and its Bruker translation
  `gradafterrefocus.lwg`, at Leo's request, for house conventions and
  timing methodology; deliberately keeps a plain hard 180 (not the
  composite 90/270 used by the imaging family) and a one-sided-k-space
  UTE-style read gradient (not GradAfterRefocus's bipolar pair), per the
  design notes in the file itself. A matching Bruker pulse programme
  (`pgse_img.lwg`) was written alongside it, at
  `/opt/topspin4.5.0/exp/stan/nmr/lists/pp/user/pgse_img.lwg` (outside
  this repo -- lives with the instrument's other user pulse programmes,
  e.g. `diffprof.lwg`), reusing `diffprof.lwg`'s dephase/rephase gradient
  syntax and `gradafterrefocus.lwg`'s `d24`/l-register-logging
  conventions. `time_calculation()` verified via the mock hardware
  harness against an independent term-by-term re-implementation (exact
  match) and against a deliberately negative-delay case (correctly
  raises). See `docs/validation-notes/LWG_PGSE-Image_H.md` and `_X.md`
  for what's still unverified -- not yet run on hardware.
- `.gitattributes` marking `sequences/**/*.py` as `-text`, so git never
  converts line endings in either direction: the CRLF bytes committed are the
  bytes checked out on every platform (including a Windows clone, where Git
  for Windows defaults to `core.autocrlf=true`) and the bytes GitHub serves
  for a raw download, which is often how a sequence reaches the instrument
  PC. Deliberately not `text eol=crlf`, which would store LF and break raw
  downloads; deliberately no `* text=auto` catch-all. Note this guards
  against *git* converting endings -- it cannot stop an editor rewriting a
  file on save, which is what caused the flip above.

### Changed
- **Low Gamma gradient calibration is now surfaced to the operator.**
  `MAXGRAD_TABLE`'s `LOWGAMMA` row has held real measured numbers since
  18-19/08/2026 (x/y/z = 18.678/16.038/58.112 G/cm), but the `Probe`
  Parameter's description -- the only thing visible in SpinFlow's panel --
  still read `Low Gamma(not yet calibrated)`, and nothing warned at run time.
  Selecting Low Gamma therefore used a provisional z (PGSE-only) and x/y
  (PGSTE) calibration for quantitative b-values and Hz-to-mm conversion while
  the panel said it did not exist, with the caveat living only in a source
  comment. Across all 22 gradient sequences: added `PROVISIONAL_CALIBRATION`
  (keyed like `MAXGRAD_TABLE`, value = the caveat, so the *code* can act on
  it rather than a comment), a `report_probe_gradient()` WARNING that prints
  it to the SpinFlow console, and the description now reads
  `Low Gamma(provisional)`. No calibration numbers changed and no timing
  changed; verified under the mock harness that the warning fires for Low
  Gamma and not for HFX in all 22, both probes, 44 runs clean.
- `sequences/**/*.py` restored to CRLF. Sixteen files (4 diffusion, 12
  imaging) had been saved as LF at some point, which the X-Pulse compiler
  rejects with no useful error. Restoring them collapsed a ~14,800-line
  working-tree diff to the 16 lines that were a real change (the LOWGAMMA
  numbers above). All 80 sequence files are now pure CRLF.

### Open questions
- **Gradient calibration disagreement, unresolved.** The vendor PFGSTE/
  DPFGSTE `Gx_max`/`Gy_max`/`Gz_max` defaults (0.01, 0.01, 0.05 T/m) are
  ~12x smaller than this repo's measured `MAXGRAD_TABLE` (HFX x/y/z =
  11.879/11.978/57.915 G/cm = 0.1188/0.1198/0.5792 T/m) -- ~140x in `b`,
  which scales as `G^2`. The measured table is used, and the discrepancy is
  flagged in each file's design notes rather than silently resolved. A single
  measurement on a known-D standard settles it, and
  `tools/spinflow_bvalue_diff.py --solve-gmax` does the arithmetic. Because b
  is factorised (above), data acquired before this is resolved can be
  re-costed afterwards without re-running it.
- `LWG_PGSTE_H.py`'s design-note approximation for the effective gradient
  duration, `delta = P.delta + 2*RampTime`, is off by one ramp time: a
  symmetric trapezoid's area is `G*(delta + RampTime)`, so `delta_eff =
  delta + RampTime`. (The vendor sequences' own `jc_smallDELTA` already uses
  `RampTime + GradientOnTime`, i.e. the correct one.) Comment-only in that
  file -- no code uses the wrong value -- but worth correcting there.
