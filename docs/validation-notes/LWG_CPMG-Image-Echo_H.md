# Validation note: LWG_CPMG-Image-Echo_H

- **File:** `sequences/imaging/LWG_CPMG-Image-Echo_H.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample)

## What was checked
- Written by Claude (AI-assisted), not yet run on hardware.
- Timing formulas (inter-echo branch-sync padding, per-echo TE) checked
  algebraically by hand and against a standalone Python numeric sanity check
  (default Parameter values -> all `safe_delay()` arguments positive, echo
  times monotonically increasing and physically plausible) -- NOT checked
  against real hardware timing or a reference measurement.
- Python syntax checked (`ast.parse`), no raw un-escaped `&` in Parameter
  strings, file saved with CRLF line endings.

## Results / evidence
- None yet -- draft, unrun.

## Known limitations or open questions
1. **Biggest open question:** this sequence calls `Receiver1()` `NumEchoes`
   times within a single scan and correspondingly
   `TX0.setup_receive(NumEchoes*ReceiverPoints, ...)`. No other sequence in
   this repo does this (checked via grep across `sequences/**/*.py` before
   writing this file) -- there is no local precedent or firebird framework
   source available to confirm this actually delivers one
   `NumEchoes*ReceiverPoints`-length array per scan (echo-major) rather than
   failing to compile, only capturing one echo, or something else.
   `CallBack1D.process_data()` defensively logs an error if the returned
   array's LENGTH doesn't match expectations, but cannot detect echoes
   silently overwriting each other in a correctly-sized buffer. **First
   thing to check on hardware:** run with `NumEchoes=2`, `NS=1`, `DS=0`,
   plot both halves of the returned array, and confirm they look like two
   distinct, correctly-timed spin echoes before trusting anything from this
   sequence.
2. Inter-echo timing (the padding that keeps the gradient and RF branches
   synchronised echo-to-echo) was derived algebraically from
   `LWG_1D-Image-Echo_H.py`/`LWG_1D-Image-Echo-GradAfterRefocus_H.py`'s own
   validated per-echo formulas, not copied from a working multi-echo
   reference (none exists in this repo). Check the actual profile shape at
   each echo (should stay a clean, positive, correctly-positioned spin
   echo, not drifting or truncated) across the full echo train, not just
   the first one or two.
3. Duty-cycle estimate for RF and gradient scales with `NumEchoes` -- for
   large `NumEchoes` this could be a real constraint; check
   `estimate_duty_cycles()`'s logged values against actual instrument specs
   (obtained from Oxford Instruments) before running long trains
   unattended.
4. Motivation/intended use: to separate T2 from diffusion contamination in
   the octanol-water interfacial T2(x) measurements (see project chat,
   2026-08-11) -- the per-echo gradient geometry is deliberately
   diffusion-compensated (both dephase and readout lobes on the post-
   refocus side, back-to-back), so this should give a "purer" T2 decay than
   the diffprof.lwg-style pseudo-2D TAU array. That comparison (this
   sequence's T2(x) vs. the diffusion-weighted array's apparent T2(x)) is
   itself part of what still needs validating once real data exists.
