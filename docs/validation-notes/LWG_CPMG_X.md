# Validation note: LWG_CPMG_X

- **File:** `sequences/relaxation/LWG_CPMG_X.py`
- **Status:** draft
- **Validated by:**
- **Date:**
- **Hardware:** (instrument, probe, sample, target nucleus)

## What was checked
- Written by Claude (AI-assisted), not yet run on hardware.
- Mechanical X-channel port of `LWG_CPMG_H.py` (all
  Channel1/Transmit1/Receiver1/TX0 -> Channel2/Transmit2/Receiver2/TX1,
  per the vendor's own H/X sequence-pair convention) -- see that file's
  own validation note for the full timing-formula derivation/cross-checks
  (identical algebra, only the hardware channel differs).
- Ran the actual file through the mock hardware harness (which includes
  Channel2/Transmit2/Receiver2/TX1 stubs) at default Parameters (NECH=8)
  -- completed without error, scan time (32148556.64us) matches
  `LWG_CPMG_H.py`'s own result at the same default Parameters exactly.
- Confirmed `time_calculation()` at NECH=1 matches `LWG_Hahn-Echo_X.py`'s
  own `time_calculation()` exactly (32085564.64us at matching Tau=2000us)
  -- same cross-check as the H-channel pair.
- Python syntax checked, no raw un-escaped `&` in Parameter strings, file
  saved with CRLF line endings.

## Results / evidence
- None yet -- draft, unrun on real hardware.

## Known limitations or open questions
1. Never run on the instrument.
2. **X-CHANNEL CALIBRATION**: `SF` defaults to a placeholder 15.01 MHz;
   `P90`/`P180`/`RFA0` (TXAmplitude)/`RA` (ReceiverAttenuation) are still
   the 1H-calibrated values ported from `LWG_CPMG_H.py` and are almost
   certainly WRONG for whatever nucleus is actually on the X channel.
   Recalibrate SF/P90/P180/RFA0/RA for the target nucleus before running
   this on a sample -- see the X-CHANNEL CALIBRATION note at the top of
   the file.
3. Same RF duty-cycle caveat as `LWG_CPMG_H.py` applies (duty scales with
   NECH; `MaxRFDuty` threshold is a conservative placeholder, not a
   vendor-confirmed rating) -- check before running a large-NECH train
   unattended.
