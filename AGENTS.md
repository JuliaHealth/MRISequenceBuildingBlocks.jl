# Migrated sequence-builder guidance

## Purpose

This directory is a clean staging area for high-level MRI sequence construction
on top of KomaMRI. KomaMRI owns event and sequence mechanics; this code should
express recognizable sequence-design concepts such as Cartesian readouts,
balanced rewinders, spin-echo timing, preparation modules, and acquisition
ordering.

The directory is intentionally not a Julia package. It is expected to move into
a separate high-level package later, so keep it self-contained, free of
top-level execution, and easy to transplant.

## Current contents

| File | Responsibility |
| --- | --- |
| `gradient_design.jl` | Shared minimum-time trapezoid timing. |
| `cartesian.jl` | Reusable Cartesian line ordering and shot grouping. |
| `excitation.jl` | Centered and custom slice-selective sinc excitations. |
| `gre.jl` | Cartesian GRE and spoiled-GRE readout kernels. |
| `epi.jl` | Partial-Fourier, multishot Cartesian EPI readout kernel. |
| `bSSFP.jl` | Balanced Cartesian readout kernel and full linear bSSFP train. |
| `spinEcho.jl` | Refocusing block and readout-independent spin-echo assembly. |
| `grase.jl` | Minimum-time GRASE assembly from refocused EPI echo groups. |
| `tse.jl` | Linear or center-out turbo spin-echo acquisition assembly. |
| `spi.jl` | Three-dimensional SPI builder and Cartesian sampling orders. |
| `preprocess.jl` | Materialization of gradients that cross block boundaries. |
| `utils.jl` | Small sequence-composition utilities. |

Use this include order when loading the directory directly:

```julia
include("gradient_design.jl")
include("utils.jl")
include("cartesian.jl")
include("preprocess.jl")
include("excitation.jl")
include("gre.jl")
include("epi.jl")
include("bSSFP.jl")
include("spinEcho.jl")
include("grase.jl")
include("tse.jl")
include("spi.jl")
```

## Design philosophy

### Model sequence concepts directly

Prefer functions named after meaningful design operations. A helper such as
minimum-time lobe design or k-space-center timing is appropriate; a generic
calculation framework used by only one caller is not. Keep short, single-use
logic inside its caller unless extracting it makes the sequence easier to
understand.

Do not introduce plan objects or configuration types merely to replace
closures. The current callable kernels are concise and useful. A callable
kernel type may be worthwhile later if inspection, serialization, or dispatch
on the kernel itself becomes important.

### Compose constructors

Reusable readout designers should return a small kernel interface:

- `readout(args...)` constructs a concrete sequence for one encoding;
- `center_time(args...)` reports the k-space-center time relative to that
  readout's beginning;
- optional data such as `encodings` should be included only when a higher-level
  builder actually consumes it.

Higher-level builders should accept excitation, preparation, or readout
constructors where that keeps acquisition-specific decisions outside the
builder. Do not encode construction state in `seq.DEF`; definitions are written
to Pulseq files and should remain limited to necessary interoperability or
acquisition metadata.

### Use KomaMRI's public constructors

Prefer KomaMRI `make_` and `build_` functions whenever they preserve the needed
waveform and timing:

- use `make_` for an event that will be combined with other events in a block;
- use `build_` for a complete one-block sequence;
- use direct event constructors only for custom sampled RF/gradient waveforms,
  deliberate cross-block overlap, or behavior the public constructor cannot
  express.

Do not add compatibility wrappers around a KomaMRI constructor unless the
wrapper adds sequence-design meaning. Small floating-point differences, such
as one ULP, are acceptable.

### Use SI internally

Public design values and internal calculations use SI scalars:

- time: seconds;
- distance and FOV: metres;
- gradient amplitude: tesla per metre;
- gradient moment: tesla-seconds per metre;
- RF amplitude: tesla;
- frequency: hertz;
- phase and flip angle: radians.

Attach Unitful quantities at KomaMRI constructor boundaries. A function that
deliberately mirrors a Unitful KomaMRI API, such as
`build_centered_sinc_pulse`, may retain that external interface.

### Treat timing and moments as invariants

Sequence output must respect RF, gradient, ADC, and block rasters as well as
scanner amplitude and slew limits. In particular:

- full Cartesian GRE readouts must acquire `kx=0`;
- a centered readout must put `kx=0` at its temporal center;
- bSSFP must have zero net gradient moment over a TR and place k-space center
  halfway between adjacent RF centers when full readout is used;
- partial Fourier must retain the k-space-center sample;
- spin-echo timing is measured between RF centers and the actual readout-center
  sample;
- GRASE requires equal odd-length M0-refocused EPI groups, with every group
  center placed on its corresponding spin echo and `ky=0` acquired at the
  derived effective TE;
- TSE requires equal-duration M0-refocused Cartesian lines, uniform echo
  spacing, unique phase-encoding coverage, and `ky=0` acquired on the selected
  echo; every echo train is either triggered or padded to a requested TR;
- SPI phase encoding must be rewound before the next TR, with the requested
  spoiler moment added to that rewind;
- EPI gradients that intentionally cross block boundaries must be materialized
  before Pulseq export.

Do not hide raster corrections in unexplained tolerances. Name the physical or
format constraint that requires the correction.

### Keep labels explicit

Callers use centered Cartesian indices. Pulseq labels use zero-based array
indices. Current conventions are:

- `LIN`: phase-encoding line;
- `PAR`: partition;
- `SLC`: third SPI coordinate;
- `REV`: readout polarity;
- `SEG`: EPI segment or polarity state;
- `NAV`: navigator state.

Only add `AVG`, `SET`, `PHS`, or other labels when they describe an actual
acquisition dimension. Do not use `AVG` to distinguish ADC dwell samples within
one event.

## Editing workflow

1. Read `README.md` and `TODO.md` before changing public behavior.
2. When recovering functionality from the original project, inspect the live
   source working copy, including uncommitted and untracked files. Use it as
   design evidence, not as code to copy wholesale.
3. Make the smallest change that expresses the requested sequence concept.
4. Update the affected docstring and public-surface documentation with any API
   change.
5. Construct a small representative sequence and run `check_timing` and
   `check_hw_limits`.
6. Verify sequence-specific invariants numerically. When export behavior is in
   scope, write and read back a small Pulseq file as an additional check.

Keep the Julia REPL alive during iterative work. Use the project environment and
reload edited files with Revise or `include` as appropriate.

## Deliberate boundaries

- Do not add generated `.seq` files, plotting scripts, scanner-specific recipes,
  or top-level sequence generation to these library files.
- Gropt-based implementations are outside this migration.
- The original tests are not migration inputs; add validation in the eventual
  package around the clean public behavior.
- Phantom conversion and Koma simulation setup are adjacent workflows, not
  sequence-construction primitives.
- Avoid accumulating large diagnostic dictionaries in `seq.DEF`.
