# Changelog

All notable changes to this repo are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Initial repo scaffold: README, LICENSE, validation status legend, folder
  structure for diffusion / relaxation / imaging / calibration sequences.
- `sequences/imaging/LWG_CPMG-Image-Echo_H.py` (draft): multi-echo (CPMG)
  1D image-echo sequence -- NumEchoes gradient-refocused spin echoes per
  scan from a single excitation, each echo's dephase/readout gradient
  lobes kept diffusion-compensated (back-to-back, post-refocus) to give a
  T2(x) decay minimally confounded by diffusion, for separating T2 from
  diffusion contributions near a liquid-liquid interface. See the file's
  own header and `docs/validation-notes/LWG_CPMG-Image-Echo_H.md` for the
  (currently unverified) assumptions this relies on -- not yet run on
  hardware.
