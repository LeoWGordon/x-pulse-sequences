# Validation note: LWG_Hahn-Echo_H

- **File:** `sequences/relaxation/LWG_Hahn-Echo_H.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample)

## What was checked
- Written by Claude (AI-assisted), not yet run on hardware.
- Basic Hahn (spin) echo: hard-90 - TAU - hard-180 - TAU - echo, single
  acquisition, no gradients/imaging. Built as the direct comparison
  counterpart to `LWG_CPMG_H.py` (see that file's own note) -- arraying
  TAU directly, vs. CPMG's fixed-TAU/arrayed-NECH approach.
- `pulse()`/`_first_gap()`/the post-refocus-wait formula and the 8-step
  Meiboom-Gill phase cycle (PH1/PH2/PHRX) are copied verbatim from the
  vendor's own `CPMG_H.py` (via `sequences/imaging/LWG_CPMG-Image-Echo_H.py`'s
  already-documented provenance chain), not re-derived.
- `time_calculation()` derived exactly term-by-term by hand (every
  `Delay()`/duration-argument instruction in `run()`'s per-scan body
  summed and algebraically simplified), then independently cross-checked
  two ways:
  1. Numerically against a from-scratch Python re-implementation of the
     same exact term-by-term sum (not the simplified closed form) across
     2000 random parameter combinations (N/A here, N=1 only) -- max
     floating-point discrepancy ~2e-9us.
  2. Against `LWG_CPMG_H.py`'s own `time_calculation()` at `NECH=1` --
     this sequence IS that file's NECH=1 case (identical instruction
     sequence), and the two formulas were confirmed to give the EXACT
     same result (32085564.64us at matching default Tau=2000us),
     both symbolically (by hand) and by actually running both files'
     `time_calculation()` through the mock harness.
- `safe_delay()`'s negative-delay guard rail confirmed to actually fire
  (raises `ValueError` with a clear message) when TAU is set too short
  for `_first_gap()` to stay positive (tested with TAU=1.0us).
- Python syntax checked (module loads/runs cleanly via `runpy`), no raw
  un-escaped `&` in Parameter strings, file saved with CRLF line endings.
- Ran cleanly through the mock hardware harness at default Parameters.

## Results / evidence
- None yet -- draft, unrun on real hardware.

## Known limitations or open questions
1. Never run on the instrument -- timing formulas are algebraically and
   numerically self-consistent (see above) but have not been checked
   against a real echo/oscilloscope trace.
2. As with any Hahn-echo-with-TAU-arrayed measurement, results at long
   TAU will include some diffusion attenuation on top of true T2 decay
   (see the file's own "HAHN ECHO vs. CPMG" design note) -- this is
   expected/by design, not a bug, but should be kept in mind when
   interpreting a fitted T2 from this sequence alone.
3. Default `TAU=2000us`/`RD=2000000us` (2s) are starting points, not
   calibrated for any particular sample -- adjust RD to >=5x the sample's
   expected T1, and choose a TAU array spanning several T2 half-lives.
