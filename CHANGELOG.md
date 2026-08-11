# Changelog

All notable changes to this repo are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Initial repo scaffold: README, LICENSE, validation status legend, folder
  structure for diffusion / relaxation / imaging / calibration sequences.
- `sequences/imaging/LWG_CPMG-Image-Echo_H.py` (draft): multi-echo (CPMG)
  1D image-echo sequence, for separating T2 from diffusion contributions
  near a liquid-liquid interface. v2.0: applies NECH refocusing pulses per
  scan (timing/phase-cycle copied verbatim from the vendor's own
  `CPMG_H.py`) and images only the last one (diffusion-compensated
  gradient geometry copied from `LWG_1D-Image-Echo-GradAfterRefocus_H.py`)
  -- run as a quasi-2D sweep of NECH to build a T2(x) curve with short,
  fixed echo spacing. Supersedes a v1.0 draft that instead tried to
  acquire every echo from a single excitation (no precedent for that
  anywhere in this repo; not carried forward). See the file's own header
  and `docs/validation-notes/LWG_CPMG-Image-Echo_H.md` for what's still
  unverified -- not yet run on hardware.
