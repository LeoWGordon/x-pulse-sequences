# Validation note: LWG_CPMG-Image-Echo_H

- **File:** `sequences/imaging/LWG_CPMG-Image-Echo_H.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample)

## What was checked
- Written by Claude (AI-assisted), not yet run on hardware.
- v1.0 (superseded) called `Receiver1()` `NumEchoes` times per scan to
  acquire every echo from one excitation -- no precedent for this anywhere
  in the repo.
- v2.0 (current) rebuilt around the vendor's own `CPMG_H.py` (supplied by
  L. Gordon, 2026-08-11), which images only the LAST of `NECH` refocusing
  pulses per scan and gets the T2 decay curve from a quasi-2D sweep of
  `NECH` across separate scans. The `NECH-1` blind refocusing pulses'
  timing formulas and 8-step phase cycle are copied verbatim from
  `CPMG_H.py`; the final (imaged) echo's gradient geometry and RF timing
  are copied verbatim from `LWG_1D-Image-Echo-GradAfterRefocus_H.py`.
  There is now only one `parallel` (gradient vs RF) block per scan, for
  the final echo only -- no multi-Receiver1-per-scan assumption, no
  hand-derived inter-echo branch-padding.
- Timing formulas checked algebraically by hand and against a standalone
  Python numeric sanity check at default Parameter values, for
  `NECH = 1, 2, 4, 8, 16`: all `safe_delay()` arguments positive, echo TE
  increases monotonically and by approximately `2*TAU` per additional
  echo once `NECH` is large enough that the loop dominates -- consistent
  with the intended CPMG timing.
- Python syntax checked (`ast.parse`), no raw un-escaped `&` in Parameter
  strings, file saved with CRLF line endings.

## Results / evidence
- None yet -- draft, unrun.

## Known limitations or open questions
1. The transition point from the last blind pulse into the final,
   gradient-carrying echo is new -- not present in `CPMG_H.py` (which has
   no gradient at all) or, as a multi-pulse train, in
   `LWG_1D-Image-Echo-GradAfterRefocus_H.py` (which only ever has one
   refocusing pulse). Check gradient/RF synchronisation on the final echo
   specifically, e.g. by comparing the profile shape/position at
   `NECH=1` against a same-Parameters run of
   `LWG_1D-Image-Echo-GradAfterRefocus_H.py`, where the two sequences'
   final-echo timing should be nearly identical.
2. The final echo's post-refocus RF wait deliberately does NOT reuse
   `CPMG_H.py`'s own pre-acquisition formula (which uses
   `ReceiverFilter.group_delay`, tuned for a plain spectroscopy FID, not
   for centring a gradient echo) -- it reuses
   `LWG_1D-Image-Echo-GradAfterRefocus_H.py`'s formula instead. This
   reasoning is written out in the file header ("WHY THE FINAL ECHO
   DOESN'T REUSE CPMG_H.py's OWN PRE-ACQUISITION FORMULA") but has not
   been checked against real hardware timing.
3. Duty-cycle estimate: RF duty scales with `NECH` (`P90` once +
   `NECH x P180`); gradient duty does not (only the final echo drives the
   gradient coil). Check `estimate_duty_cycles()`'s logged values against
   actual instrument specs before running long trains (large `NECH`)
   unattended.
4. Intended use: run as a quasi-2D sweep of `NECH` (e.g. 1, 2, 4, 8, 16,
   ...) at a small, fixed `TAU`, to get a T2(x) decay curve with short,
   fixed echo spacing (`2*TAU`) that's less diffusion-attenuated than the
   existing `diffprof.lwg`-style single-echo `TAU` array (where TE is
   increased by increasing `TAU` directly). Comparing the two T2(x) maps
   (this sequence vs. the diffusion-weighted array) is itself part of
   what still needs validating once real data exists from both.
