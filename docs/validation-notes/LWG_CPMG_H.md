# Validation note: LWG_CPMG_H

- **File:** `sequences/relaxation/LWG_CPMG_H.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample)

## What was checked
- Written by Claude (AI-assisted), not yet run on hardware.
- Basic Carr-Purcell-Meiboom-Gill (CPMG) multi-echo train: hard-90, then
  NECH hard-180 refocusing pulses at fixed short spacing 2*TAU, single
  acquisition after the LAST refocusing pulse only, no gradients/imaging.
  Built as the direct comparison counterpart to `LWG_Hahn-Echo_H.py`
  (arraying NECH at fixed short TAU, vs. Hahn echo's arrayed-TAU
  approach) at the user's request.
- Structurally identical to `sequences/imaging/LWG_CPMG-Image-Echo_H.py`'s
  own NECH-1 "blind" refocusing-pulse loop, but with that file's
  imaging/gradient machinery removed from the final echo too -- this file
  just acquires a plain FID after the final pulse, exactly like the
  vendor's own (non-imaging) `CPMG_H.py` (Created 10/09/2013, Author AS,
  supplied by L. Gordon 11/08/2026). `pulse()`/`_first_gap()`/
  `_subsequent_gap()`/the post-refocus-wait formula/the 8-step
  Meiboom-Gill phase cycle (PH1/PH2/PHRX) are all copied verbatim from
  that already-documented provenance chain, not re-derived.
- `time_calculation()` derived exactly term-by-term by hand for general
  NECH, then cross-checked three ways:
  1. Numerically against a from-scratch Python re-implementation of the
     exact (unsimplified) term-by-term sum across 2000 random trials
     spanning NECH in {1,2,3,4,8,16,32} and randomised Tau/TXEnableTime/
     P90/P180/dead_time/group_delay/Dead1/ReceiverTime -- max
     floating-point discrepancy ~1.9e-9us, i.e. exact agreement.
  2. At NECH=1, confirmed to reduce to EXACTLY `LWG_Hahn-Echo_H.py`'s own
     independently-derived `time_calculation()` formula (both by hand and
     by running both files' actual `time_calculation()` through the mock
     harness with matching Tau=2000us: both return 32085564.64us).
  3. Ran the actual file (not just the formula) through the mock harness
     at NECH=1,2,4,8,16,32 -- all completed without error, scan times
     scale monotonically and sensibly with NECH.
- `safe_delay()`'s negative-delay guard rail confirmed to fire correctly
  for an impossible TAU (shared code path with LWG_Hahn-Echo_H.py, tested
  there).
- Python syntax checked (module loads/runs cleanly via `runpy`), no raw
  un-escaped `&` in Parameter strings, file saved with CRLF line endings.

## Results / evidence
- None yet -- draft, unrun on real hardware.

## Known limitations or open questions
1. Never run on the instrument -- timing formulas are algebraically and
   numerically self-consistent (see above, including exact agreement with
   both an independent re-derivation and the sibling Hahn-echo file) but
   have not been checked against a real multi-echo trace.
2. RF duty cycle scales with NECH (P90 once + NECH x P180 per scan) --
   `estimate_rf_duty_cycle()` warns if `MaxRFDuty` is exceeded, but that
   threshold is a conservative, editable PLACEHOLDER, not a vendor-
   confirmed transmitter rating (see the pulse-programme-parameters
   skill). Check with Oxford Instruments before running a large-NECH
   train unattended.
3. Default `TAU=500us` (short, fixed -- intentional) / `RD=2000000us` (2s)
   / `NECH=8` are starting points, not calibrated for any particular
   sample. TAU should be set as short as the hardware/probe ring-down
   allows; array NECH (e.g. 1,2,4,8,16,32,...) to build the decay curve.
