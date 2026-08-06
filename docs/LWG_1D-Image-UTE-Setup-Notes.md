# LWG_1D-Image-UTE_H/X — Setup Notes

Practical setup guide for `LWG_1D-Image-UTE_H.py` / `_X.py`. Read the file's own header first — this covers what to actually set, in what order, before your first real acquisition.

## What's different about this sequence

No refocusing pulse. The read gradient ramps up and settles **before** the excitation pulse, the pulse fires with the gradient already flat, and acquisition starts as soon as fixed hardware dead times allow. There's no TAU — TE is essentially the sum of `P90/2 + PostPulseOverhead + Dead1 + ReceiverFilter.dead_time`. It also acquires one-sided k-space (k=0 at the pulse, out to +kmax), so the reconstruction is a magnitude-mode projection, not a phased spectrum — fine for locating short-T2/T2* signal, not a substitute for a phased profile if you need one.

## 1. Hardware prerequisites

- Probe suited to your target nucleus and short-T2/T2* species, tuned/matched at `SF`+`O1`.
- Gradient coil connected and calibrated (see §3) on whichever axis you'll use.
- Know your transmitter's real dead-time/ring-down behaviour at short pulse widths — this sequence pushes right up against it.

## 2. TE floor: determining Dead1 and dead_time

TE here is not a knob you set directly — it falls out of four terms, logged every run as `MinimumTE`:

```
MinimumTE = P90/2 + PostPulseOverhead(2us, fixed) + Dead1 + ReceiverFilter.dead_time
```

- **`Dead1`** (probe ring-down) is a real hardware property, not a free parameter to shrink at will. Determine it empirically: run a simple hard-pulse FID with acquisition starting increasingly close to the pulse and watch for ringing/receiver overload artefacts near the start of the FID. The shortest `Dead1` that stays clean is your floor — don't set it lower than what you've verified on your own probe.
- **`ReceiverFilter.dead_time`** is set by which `Filter` (spectral width) you choose — wider filter = shorter dead_time = shorter TE, at the cost of more noise bandwidth. Check the logged `MinimumTE` for a couple of `Filter` settings and pick the narrowest one that still gets you the TE you need.
- **`P90`** — a smaller flip angle (shorter pulse and/or lower `RFA0`) both shortens TE slightly (TE includes `P90/2`) and permits faster repetition. Tune to your T1/SNR needs, not just to minimise TE — too small a flip angle costs you signal.

Practical order: pick `Filter` first (bandwidth/noise trade-off), then measure `Dead1` on your probe, then choose the smallest `P90` your SNR budget tolerates, then read off `MinimumTE` from the log and confirm it's short enough for your species' T2*.

## 3. Gradient calibration workflow

`FPX`/`FPY`/`FPZ` are the *only* place per-axis gradient calibration is applied (via `GradientMatrix()`), separate from `G1` (the logical -1..1 strength knob you vary per acquisition). Calibrate each axis you intend to use once, then leave `FPX/FPY/FPZ` fixed:

1. Run a phantom of known geometry with `G1` at a moderate, safe value and `GradAxis` set to the axis under test.
2. Compare the reconstructed profile's spatial extent to the phantom's known physical dimension.
3. Adjust `FPX`/`FPY`/`FPZ` (for that axis) until the reconstructed extent matches the true dimension.
4. Repeat per axis — the three scalers are independent.
5. Re-check periodically or after any hardware/coil change; gradient coil calibration can drift.

`D70` (RampTime), `D73` (GradSettle), and `D75` (PreGrad) don't affect TE (PreGrad happens well before the pulse, and the gradient is already flat by the time the pulse fires) — set them long enough for a clean, eddy-current-settled plateau, short only where you're pushing scan time down.

## 4. Phantom validation before real samples

1. Confirm gradient calibration (§3) first.
2. Run the UTE sequence on the same phantom used for calibration and check the profile width/position matches expectation on all axes you plan to use.
3. Compare against a spin-echo profile (`LWG_1D-Image-Echo_H/X.py`) of the same phantom if your phantom's T2 is long enough to give a clean echo — the two profiles should agree on geometry; differences beyond expected T2*/susceptibility effects suggest a timing or calibration issue worth chasing down before moving to short-T2 samples.
4. Only then move to your actual short-T2/T2* sample.

## 5. Duty-cycle safety recap

`estimate_duty_cycles()` logs estimated RF and gradient duty cycle every run and warns if either exceeds `MaxRFDuty` (default 5%) / `MaxGradDuty` (default 10%). Two things worth knowing:

- These thresholds are **conservative placeholders**, not a vendor-confirmed rating for the X-Pulse RF/gradient amplifiers — confirm real limits with Oxford Instruments before relying on them for unattended runs.
- The gradient here is held on for the **entire scan** (pulse + dead time + acquisition), not just a short lobe like in the spin-echo sequences — so gradient duty cycle is easy to hit with short `RD` or many `NP`. If you get a gradient-duty warning, lengthen `RD`, reduce `NP`, or reduce `G1` before ignoring it.

## 6. Mains lock

`UseMainsLock=0` by default, same as the other imaging sequences in this set — and *especially* worth leaving off here specifically, since a mains-lock trigger adds a variable pre-pulse delay (up to one AC half-cycle) directly onto this sequence's whole point, which is a short, tight, reproducible TE. Only enable it if you've confirmed your X-Pulse needs mains-synchronous triggering, and double check `MainsLockChannel` — the manual doesn't list a confirmed channel for X-Pulse.

## Quick checklist before your first real scan

1. Probe tuned/matched.
2. `Filter` chosen for your noise/TE trade-off.
3. `Dead1` measured on your probe, not guessed.
4. `P90`/`RFA0` set for adequate flip angle at that pulse width.
5. `MinimumTE` (logged) short enough for your species.
6. `FPX`/`FPY`/`FPZ` calibrated per axis against a known phantom.
7. Phantom profile validated, ideally cross-checked against a spin-echo profile.
8. Duty-cycle warnings clear (or explicitly accepted with vendor confirmation).
9. `UseMainsLock=0` unless you have a specific reason to change it.
