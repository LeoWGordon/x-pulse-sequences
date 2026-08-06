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
