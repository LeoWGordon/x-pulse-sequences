# Selective-Pulse Power Calibration

For `LWG_1D-Image-Echo-Selective_H/X.py` and `LWG_Selective-Echo_H/X.py`.

## What was wrong, and what changed

The v1.0–v3.0 selective-echo sequences ran the shaped pulses on the **high-power (HP)** transmit port at `RFAsh0`/`RFAsh1` = 0.10, an uncalibrated number that was never validated against a real flip angle. At v4.0 the sequences switch to the **low-power (LP)** port instead, matching the vendor's own long-shaped-pulse sequence (`WET-FID_H.py`). LP tops out at roughly 10% of HP's full scale (your measured figure, now the `LPMaxFraction` parameter), so `RFAsh0`/`RFAsh1` now mean something different from before — they're a fraction of LP's own, already-reduced, maximum. A badly wrong flip angle from the old HP/uncalibrated combination is the most likely explanation for the selective image "burning a hole" rather than forming a clean echo: instead of a coherent 90–180 spin echo, an over- or under-rotated pulse pair mostly saturates or dephases signal at the target frequency, which looks like a notch rather than a refocused peak.

`LWG_Selective-Echo_H.py` / `_X.py` are new, gradient-free pulse programmes built specifically to re-calibrate this in isolation, before combining the shaped pulses with imaging gradients again.

## Three ways to find the right power

### 1. Nutation calibration (the ground truth — always do this at the end)

This is the standard, reliable way to calibrate any shaped pulse on real hardware, and it doesn't depend on anything from the Geen & Freeman paper beyond the shapes themselves.

1. Acquire an ordinary hard-pulse spectrum of your sample first, so you know where your target resonance sits relative to `SF`+`O1`.
2. In `LWG_Selective-Echo_H.py` (or `_X.py`), set `PulseOffset` to that resonance's offset in Hz.
3. Temporarily reduce the refocusing pulse's effect — set `RefShape` to `'GAUSSIAN'` with a very low `RFAsh1`, or just watch the FID right after the excitation pulse rather than the refocused echo — and sweep `RFAsh0` (or `P90sh`) upward from a low starting value. Signal amplitude should rise, peak (90°), fall through zero (180°), and reach a negative peak (270°) if you keep going. The zero-crossing at 180° pins down the calibration more precisely than the amplitude maximum alone, so it's worth sweeping past 90° the first time.
4. With the excitation pulse now fixed at its calibrated 90° condition, sweep `RFAsh1` (or `P180sh`) and watch the **refocused echo** amplitude (not the FID right after the pulse) — it should maximise at the true 180° condition.
5. Copy the calibrated `PulseOffset`, `P90sh`, `P180sh`, `RFAsh0`, `RFAsh1` across to `LWG_1D-Image-Echo-Selective_H.py` / `_X.py`'s identically-named parameters.

Repeat for both channels (`_H` and `_X`) and for each nucleus/resonance you plan to image — the calibration is specific to the shape, duration, port, and target frequency.

Whichever calculation route below you use to get a *starting point*, always finish with this nutation check — neither calculation knows your actual probe, coil, or sample.

### 2. Shape-integral calculation (recommended starting point — this is your Bruker macro, translated)

You supplied a working TopSpin pulse programme that calculates a selective pulse's power level from a calibrated hard pulse, the shape's target rotation, and the shape's "integration factor":

```
cnst11 = -10*log10(plw1)
cnst2  = cnst11 - 20*(log10((p1*totrot2)/(p2*90.0)) - log10(integfac2)) + cnst0
spw2   = pow(10,(-cnst2)/10)
```

That's a power-domain (Watts, 10·log10) calculation. X-Pulse's `RFAsh0`/`RFAsh1` are a linear **amplitude** scale (0–1, not Watts), so the translated version works directly in amplitude terms — no squaring/square-rooting needed:

```
scale_needed = RefAmplitude_HP × (P1Hard × TargetRotation) / (ShapeDuration × 90 × IntegFactor) × 10^(−PowerAdjust_dB/20)
```

- `P1Hard` / `RefAmplitude_HP` — your calibrated **hard** 90° pulse's width and relative amplitude on the HIGH-power port (Bruker's `p1`/`plw1`).
- `TargetRotation` — the rotation the shape is designed to produce: `ExRotation=90` for E-BURP-1 (excitation), `RefRotation=180` for RE-BURP (refocusing) — matches Bruker's `totrot2`.
- `IntegFactor` — the shape's average/peak amplitude ratio (Bruker's `integfac2`, what its shape tool calls "Integ. Factor"). **This is computed automatically** from the actual E-BURP-1/RE-BURP envelope the sequence synthesises — you don't need to look this up anywhere, unlike the Bruker version where it comes from a shape-tool lookup.
- `PowerAdjust_dB` — an optional manual fine-tune (dB, amplitude convention), equivalent to Bruker's `cnst0`. Leave at 0 for the first pass; only nudge it after comparing against a real nutation curve.

The result is referenced to the HIGH-power port (same port the hard pulse was calibrated on); the sequence automatically divides by `LPMaxFraction` to re-express it on the LOW-power port's own smaller scale before it goes into `RFAsh0`/`RFAsh1`.

To use it: set `PowerCalcMethod='shape'`, fill in `P1Hard`/`RefAmplitude_HP` from your normal hard-pulse calibration, leave `ExRotation`/`RefRotation` at their defaults (90/180) unless you're using different shapes, and run. The sequence logs the computed integration factors and resulting `RFAsh0`/`RFAsh1` (and warns if either exceeds 1.0). As a sanity check with typical values (5 ms E-BURP-1 excitation, 10 ms RE-BURP refocusing, `P1Hard`=9.58 µs hard pulse at `RefAmplitude_HP`=0.40, `LPMaxFraction`=10%): the calculation lands on `RFAsh0`≈0.023 and `RFAsh1`≈0.047 — much lower than the 0.30 placeholder default, which fits: a multi-millisecond shaped pulse needs far less peak power than a ~10 µs hard pulse to produce a comparable rotation.

**Then verify against §1's nutation curve.** This calculation gets you close on physical grounds, but only a real nutation curve on your actual sample/probe confirms it.

### 3. Legacy direct dB attenuation (only if you have a dB figure that ISN'T derived from a hard-pulse ratio)

Both `LWG_Selective-Echo_H.py` and `_X.py` include a small built-in calculator, `db_to_lp_relative_scale()`. Set `PowerCalcMethod='db'` (plus `RefAmplitude_HP`, `Excitation_dB`, `Refocus_dB`, `DbConvention`) to override `RFAsh0`/`RFAsh1` with a value computed directly from a dB attenuation figure, instead of the shape-integral route above. Use this only if you have a dB number that ISN'T already expressible as a hard-pulse ratio (§2 is the better route if it is, since it needs no external table at all).

**I do not have Geen & Freeman's specific dB/power table** — you gave me the Fourier-coefficient tables (E-BURP-1 Table 2, RE-BURP Table 8), not a power table, so I have nothing from the paper to plug in here. This calculator is general-purpose machinery for *your own* dB figure and *your own* calibrated reference point, not something pre-filled with literature values.

The math:

```
relative_scale = reference_relative_scale × 10^(−dB / divisor)
```

- `divisor = 20` for an **amplitude/voltage** dB convention (dB defined on B1, the usual convention in NMR pulse literature — flip angle is linear in B1, so 6 dB attenuation halves B1).
- `divisor = 10` for a **power** dB convention (dB defined on deposited power directly, as on many RF component datasheets).

Getting the convention wrong is a factor-of-2-in-dB error in the result (equivalently, a square-root error in relative amplitude), so check which one your source actually uses before trusting a number.

Because the shaped pulses run on the LP port while a dB figure from the literature is normally referenced to a HP hard-pulse calibration, the function does one more step:

```
lp_relative_scale = (reference_relative_scale × 10^(−dB / divisor)) / LPMaxFraction
```

To use it:

1. Calibrate an ordinary hard 90° pulse on the HP port as you normally would (this gives you `RefAmplitude_HP`, e.g. your usual `RFA0`).
2. Set `Excitation_dB` / `Refocus_dB` to the attenuation (relative to that hard-pulse calibration) your source specifies for the shaped pulse.
3. Set `DbConvention` to `'amplitude'` or `'power'` to match your source.
4. Set `PowerCalcMethod='db'` and run — the sequence logs the computed `RFAsh0`/`RFAsh1` values (and warns if the result exceeds 1.0, meaning LP can't actually reach that power).
5. **Verify the result against the nutation procedure above** before trusting it for real measurements — a calculated number is only as good as the dB figure and reference calibration you fed it, neither of which this file can independently check.

## Quick reference: parameters involved

| Parameter | Meaning |
|---|---|
| `PulseOffset` | Hz offset applied only during the shaped pulses (not to acquisition) — the dedicated place to target/calibrate a resonance |
| `P90sh` / `P180sh` | Shaped pulse durations (µs); also set the shape's time resolution |
| `RFAsh0` / `RFAsh1` | Relative TX power for excitation/refocusing, 0–1 of the **LP** port's own max |
| `LPMaxFraction` | Your measured LP-vs-HP power ratio, used by both calculation routes |
| `PowerCalcMethod` | `'manual'` (default): use `RFAsh0`/`RFAsh1` as entered. `'shape'` (recommended): calculate from your hard-pulse calibration + each shape's own integration factor (§2). `'db'`: legacy direct dB-attenuation route (§3) |
| `P1Hard`, `RefAmplitude_HP`, `ExRotation`, `RefRotation`, `PowerAdjust_dB` | Inputs to the shape-integral calculator, only used when `PowerCalcMethod='shape'` |
| `Excitation_dB`, `Refocus_dB`, `DbConvention` | Inputs to the legacy dB calculator, only used when `PowerCalcMethod='db'` |
