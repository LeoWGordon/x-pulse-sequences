---
name: pulse-programme-parameters
description: Conventions for writing/editing Parameter(...) blocks in this repo's X-Pulse pulse programmes (SpinFlow XML safety, panel-width-aware descriptions, refocusing naming). Load whenever adding, editing, or reviewing a Parameters block or Parameter() call in sequences/**/*.py.
---

# Pulse programme Parameter conventions

Rules learned the hard way (a real crash) from working on this repo's `sequences/**/*.py` files. Apply all of these whenever touching a `Parameter(...)` call or a file's `Parameters` block.

## 1. File encoding: CRLF, always

The X-Pulse pulse-sequence compiler only accepts Windows CRLF line endings. LF-saved files fail to compile with no useful error message (see `LWG_1D-Image-Echo_H.py`'s own changelog item 10g for the debugging story). After creating or editing any `sequences/**/*.py` file, convert/verify CRLF before considering the edit done:

```python
data = open(path, "rb").read()
data = data.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")
open(path, "wb").write(data)
```

Check with `file <path>` — it must report "with CRLF line terminators".

## 2. Escape XML special characters in Parameter strings — this WILL crash the sequence otherwise

SpinFlow transmits every `Parameter(...)`'s name/value/description as an XML packet. A raw, un-escaped `&` in any of those strings (e.g. `"Geen & Freeman 1991"`) produces a "badly formed XML packet" warning, SpinFlow silently **drops that Parameter** (it does not just mangle the text), and if `run()` relies on that Parameter's value (e.g. parsing burp coefficients), the sequence then crashes deep inside `run()` with a confusing, unrelated-looking traceback (`RuntimeError: Sequence init failed` etc.) — the actual cause is upstream and easy to miss.

- Never write a raw `&` in a `Parameter(...)` call's name, default value, or description. Use `&amp;`.
- The house style already uses HTML entities for other special characters — keep using them: `&#956;` (µ), `&#176;` (°), `&#8230;` (…).
- This applies to description strings specifically (what actually gets XML-transmitted); raw `&` in plain `#` comments or docstrings elsewhere in the file is harmless but should still be avoided/cleaned up for consistency if you're already touching that text.
- Before finishing an edit, grep the file for stray ampersands: a raw `&` not already followed by `#\d+;` or `amp;` is a bug.

## 3. Keep Parameter descriptions SHORT and front-loaded — SpinFlow's panel has very little room

SpinFlow's parameter panel does not have space to show a full paragraph per parameter. Long, discursive descriptions (the kind that read well in a docstring) get cut off or make the panel unusable. When writing a description:

- Lead with the essential, scannable info: units, valid range, and default/critical dependency (e.g. `"used if PowerCalcMethod='shape'"`) — in that order, front to back, since truncation cuts from the right.
- Cut it down to roughly one short clause. If you need to explain *why*, put that explanation in a `#` comment above the `Parameter(...)` line or in a docstring — not in the description string itself.
- Compare:
  - Too long: `"Reference HARD-pulse relative amplitude on the HIGH-power port for a KNOWN, calibrated flip angle (e.g. your normal RFA0) -- used if PowerCalcMethod='shape' or 'db'"`
  - Right-sized: `"Ref. hard-pulse rel. amplitude [0…1] on HP port, known flip angle -- for shape/db methods"`
- This is a real, applied rule now — see `sequences/selective/LWG_Selective-PulseAcquire_H.py` for a worked example of the target style across a whole Parameters block.

## 4. Never use bare "Ref" as shorthand for "refocusing"

This repo already has genuine *reference*-value parameters (`RefAmplitude_HP`, `P1Hard` — a calibration anchor point) sitting in the same Parameters block as refocusing-pulse parameters. If refocusing ones are also abbreviated "Ref" (`RefShape`, `RefRotation`, `RefBurpCoeffsA/B`), the two meanings collide and become genuinely hard to tell apart at a glance in SpinFlow's parameter list or in code.

- Always spell out **"Refocus"** for anything about the refocusing pulse: `RefocusShape`, `RefocusRotation`, `RefocusBurpCoeffsA`/`RefocusBurpCoeffsB`, `Refocus_dB`, `RefocusPulseWidth`, `RefocusAmpProfile`, `RefocusIntegFactor`, etc. — both in the Parameter's SpinFlow-facing key string (first arg to `Parameter(...)`) AND in the Python attribute/variable name, so the two never drift apart the way `RefocusShape = Parameter("RefShape", ...)` used to.
- Reserve bare "Ref"/"Reference" exclusively for genuine reference/calibration-anchor values (a known, already-calibrated pulse you're comparing against — like `RefAmplitude_HP`).
- If you're renaming an existing key (not just adding a new one), rename the Python attribute name AND every `P.<name>` usage together, and use word-boundary matching so you don't clobber unrelated identifiers that merely contain "Ref" as a substring (e.g. don't touch `RefAmplitude_HP` while renaming `RefShape`).

## 5. Naming symmetry

Keep a Parameter's Python attribute name and its SpinFlow-facing string key either identical or at least non-misleading relative to each other — don't let them drift apart the way `RefocusShape = Parameter("RefShape", ...)` did (clear in Python, ambiguous in SpinFlow). If you rename one, rename the other to match.

## 6. Check units at every arithmetic combination of two Parameters — a zero default can hide a real bug for months

Real bug found 06/08/2026, present in every `shaped_pulse()`-using file since the shaped-pulse machinery was first written: `Frequency` (the sequence's base frequency) is computed and always carried in **MHz** (`P.FrequencyBase + P.FrequencyOffset*1.0e-6 + ...` — note `FrequencyOffset`, a Hz-valued Parameter, is correctly converted with `*1.0e-6` before being added). But `shaped_pulse()`'s hardware call added `base_frequency + pulse_offset` directly — `pulse_offset` (from the `PulseOffset` Parameter) is in **Hz**, un-converted. A 500 Hz offset silently became `59.7 MHz + 500 = 559.7` "MHz", far outside the hardware's valid range, and the sequence failed to initialize with a completely generic, unrelated-looking error (`RuntimeError: Sequence init failed.`) deep inside the firmware/broker layer — nothing in the Python traceback pointed at a units mismatch.

This bug was invisible for as long as `PulseOffset` stayed at its default of `0.0` (0 Hz still equals 0 MHz, so the missing conversion never mattered) — it only surfaced the first time someone actually followed the documented calibration procedure and set a real, nonzero offset. **A default of 0 on either operand of an addition/subtraction is not evidence the units are consistent — it just means the bug hasn't been exercised yet.**

- Whenever two Parameters (or a Parameter and a hardware call's fixed argument) are combined arithmetically, explicitly check they're in the same unit before combining, and comment the conversion inline if one is needed (see the `pulse_offset*1.0e-6` fix for the exact style).
- Don't trust "it worked in testing" as proof a formula is unit-correct if every test happened to use a zero/default value for one of the operands. Test with a real nonzero value for every Parameter that's supposed to matter before considering a sequence validated.

## 7. Never name a non-cycled Parameter attribute ending in "Phase" — `PhasesManager` auto-discovers it and assumes it's a comma-separated list

Real crash found 14/08/2026 (full traceback supplied by the user from the instrument itself), in `LWG_PGSTE-WET_H.py`/`LWG_PGSE-WET_H.py`'s first version: a WET module's fixed pulse-phase Parameter was named `WetPhase = Parameter("WetPhase", 0.0, ParameterTypes.Double, ...)`. `PhasesManager(P)` — called as the very first line of `run()`, before any sequence-specific code — auto-discovers **every** Parameter attribute whose name ends in the literal substring `"Phase"` (matching the legitimate `P1Phase`/`P2Phase`/`P3Phase`/`RXPhase` convention used throughout this repo) and unconditionally tries to parse its value as a comma-separated integer list via `PhaseListContainer(value.split(","))`. It does this **regardless of the Parameter's declared type** — a `Double`/`Int32`-typed Parameter whose name merely happens to end in "Phase" gets swept in exactly the same as a real phase-cycle string, and crashes with `AttributeError: 'float' object has no attribute 'split'` (or `'int' object has no attribute 'split'`) deep inside `PhasesManager.__init__`, before `run()`'s own logic even starts.

- **Never end a Parameter's Python attribute name in "Phase"** unless it genuinely holds a comma-separated, `PhasesManager`-cycled integer-degree list (like `P1Phase`/`RXPhase`). For any other angle/phase-like value (e.g. a single fixed transmit phase for a non-cycled shaped pulse), pick a name that doesn't end in that substring — e.g. `WetAngle` (short code `WETANG`), not `WetPhase`.
- This is a **name-pattern match, not a type check** — `RephaseFraction` (short code `RephaseFrac`) is safe precisely because its attribute name ends in `"Fraction"`, not `"Phase"`, even though it contains "phase" as a substring elsewhere in the word. The ENDING is what matters, matching how `PhasesManager` appears to scan for its known suffix.
- The mock hardware harness's `PhasesManager` (see §9's reference skeleton) has been hardened to reproduce this exact check — it raises the same `AttributeError` locally for any non-string `"...Phase"`-suffixed attribute, so this class of bug is now caught before deployment rather than on the instrument. If you rebuild the mock harness from scratch, make sure to carry this check over (see `LWG_PGSTE-WET_H.py`'s or `LWG_PGSE-WET_H.py`'s changelog for the exact patched mock code).
