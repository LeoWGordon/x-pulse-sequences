# Validation note: LWG_Hahn-Echo_X

- **File:** `sequences/relaxation/LWG_Hahn-Echo_X.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample, target nucleus)

## What was checked
- Written by Claude (AI-assisted), not yet run on hardware.
- Mechanical X-channel port of `LWG_Hahn-Echo_H.py` (all
  Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1,
  per the vendor's own H/X sequence-pair convention) -- see that file's
  own validation note for the full timing-formula derivation/cross-checks
  (identical algebra, only the hardware channel differs).
- Confirmed `time_calculation()` gives the identical numeric result to
  `LWG_Hahn-Echo_H.py` at matching default Parameters (32085564.64us at
  Tau=2000us) -- expected, since only the channel-selector hardware calls
  differ, not the timing arithmetic.
- Confirmed `time_calculation()` matches `LWG_CPMG_X.py`'s own
  `time_calculation()` at NECH=1 (32085564.64us, matching Tau) -- same
  cross-check as the H-channel pair.
- Ran cleanly through the mock hardware harness at default Parameters
  (which includes Channel2/Transmit2/Receiver2/TX1 stubs).
- Python syntax checked, no raw un-escaped `&` in Parameter strings, file
  saved with CRLF line endings.

## Results / evidence
- None yet -- draft, unrun on real hardware.

## Known limitations or open questions
1. Never run on the instrument.
2. **X-CHANNEL CALIBRATION**: `SF` defaults to a placeholder 15.01 MHz;
   `P90`/`P180`/`RFA0` (TXAmplitude)/`RA` (ReceiverAttenuation) are still
   the 1H-calibrated values ported from `LWG_Hahn-Echo_H.py` and are
   almost certainly WRONG for whatever nucleus is actually on the X
   channel. Recalibrate SF/P90/P180/RFA0/RA for the target nucleus (e.g.
   via a simple pulse-acquire nutation curve) before running this on a
   sample -- see the X-CHANNEL CALIBRATION note at the top of the file.
3. Same diffusion-attenuation caveat as `LWG_Hahn-Echo_H.py` applies
   (arraying TAU directly couples some diffusion attenuation into the
   apparent T2 at long TAU) -- see that file's note and
   `LWG_CPMG_X.py` for the comparison sequence.
