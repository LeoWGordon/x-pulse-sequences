# x-pulse-sequences

A living collection of pulse sequences for the Oxford Instruments X-Pulse Broadband
Benchtop NMR Spectrometer, written for its firebird (Python 2.7) pulse-sequence
framework.

## Purpose

This repo is a working library of X-Pulse pulse programmes — new sequences as
they're developed, existing ones as they're refined, and the validation record
for each. It's meant to grow over time rather than represent a single finished
release.

## Provenance and validation

Some sequences in this repo are initially drafted with AI assistance (Claude,
Anthropic) — for example, generating a first version from a description of the
desired pulse train, or adapting an existing sequence. **No AI-assisted sequence
is run on the instrument, or trusted for data collection, until it has been
reviewed and validated by a human.** This is stated here explicitly rather than
left implicit, so that provenance is never ambiguous to anyone reading or reusing
a sequence from this repo.

Each sequence carries a validation status (see below), tracked in its own header
comment and summarized in the table for its category folder. `docs/validation-notes/`
holds a longer note per sequence describing what was checked, on what hardware,
and by whom.

### Validation status legend

| Status | Meaning |
|---|---|
| 🟡 `draft` | Written or AI-assisted; not yet run on hardware. Timing/logic not independently checked. |
| 🔵 `bench-tested` | Has been run on the instrument and produces sensible output, but results/edge cases have not been fully checked against theory or a reference sequence. |
| 🟢 `validated` | Reviewed line-by-line and confirmed correct by a named person; results checked against expected behavior or a reference measurement. Safe to use for data collection. |

A sequence's status is recorded as a comment near the top of the `.py` file, e.g.:

```python
# status: validated
# validated_by: <name>
# validated_date: <date>
```

## Structure

```
sequences/
├── diffusion/       PFG diffusion sequences (stimulated echo, profiling, etc.)
├── relaxation/      T1, T2, and related relaxation measurements
├── imaging/         Gradient-echo / spin-echo imaging sequences
├── selective/       Frequency selective sequences
└── calibration/     Pulse width, power, and other calibration sequences
docs/
└── validation-notes/  One note per validated (or in-progress) sequence
```

## Contributing / updating status

When adding a new sequence:
1. Place it in the relevant `sequences/` subfolder.
2. Mark it `draft` in its header.
3. As it's tested and reviewed, update the status comment and add/update a note
   in `docs/validation-notes/`.

## License

See [LICENSE](LICENSE).
