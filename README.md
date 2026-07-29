# Migrated sequence builders

This directory is a staging area for high-level sequence construction code. It
is intentionally not a Julia package yet. Numeric design calculations use SI
units internally; Unitful quantities are added only when calling KomaMRI pulse
constructors.

## Suggested include order

```julia
include("gradient_design.jl")
include("utils.jl")
include("cartesian.jl")
include("preprocess.jl")
include("excitation.jl")
include("gre.jl")
include("PC.jl")
include("epi.jl")
include("bSSFP.jl")
include("spinEcho.jl")
include("se_epi.jl")
include("grase.jl")
include("tse.jl")
include("spi.jl")
```

`gre.jl`, `epi.jl`, `bSSFP.jl`, and `spi.jl` use `_lobe_timing` from
`gradient_design.jl`. `bSSFP.jl` also uses `cartesian_line_order`, the GRE
kernel, and `rf_excitation`. `PC.jl` uses the Cartesian ordering, excitation,
and spoiled-GRE builders.
The public EPI kernel uses `materialize_gradients` from `preprocess.jl`.
`se_epi.jl`, `grase.jl`, and `tse.jl` use the refocusing helpers from
`spinEcho.jl`.

## Public surfaces

| File | Main functions |
| --- | --- |
| `cartesian.jl` | `cartesian_line_order` |
| `excitation.jl` | `build_centered_sinc_pulse`, `slice_selective_sinc`, `rf_sinc` |
| `gre.jl` | `gre_readout_kernel`, `sgre_base` |
| `PC.jl` | `build_pc`, `build_pc_encoding_scans` |
| `epi.jl` | `epi_readout_kernel`, `build_epi_navigator` |
| `bSSFP.jl` | `bssfp_readout_kernel`, `build_cartesian_bssfp` |
| `spinEcho.jl` | `build_refocusing_block`, `build_spin_echo` |
| `se_epi.jl` | `build_se_epi` |
| `grase.jl` | `build_grase` |
| `tse.jl` | `build_tse` |
| `spi.jl` | `spi_linear_order`, `spi_center_out_order`, `build_spi` |
| `preprocess.jl` | `materialize_gradients` |
| `utils.jl` | `rf_excitation` |

GRE, EPI, and bSSFP use a small kernel contract: `readout(args...)` builds one
concrete readout and `center_time(args...)` reports the k-space-center time
relative to its beginning. `build_spin_echo` consumes this contract without
branching on the readout family.

`epi_readout_kernel` partitions phase-encoding lines into interleaved or
contiguous echo groups. It designs common worst-case prephaser, blip, and
rewinder timing across those groups, returns each group's centered line indices
as `lines`, and rewinds every group to zero gradient zeroth moment by default.
`echo_center_time` reports the temporal center of a group; for equal odd-length
groups this is the same central-line ADC time for every shot.
`build_epi_navigator` adds an independent excitation to alternating `ky=0`
readouts designed by the same EPI kernel. Its explicit x rewinder permits any
positive navigator-line count, and its required `post_delay` separates the
navigator from imaging. The final navigator line carries `AVG=1`; earlier lines
carry `AVG=0`. The resulting sequence can be prepended once to an EPI
acquisition or supplied once at sequence start to `build_grase` or `build_tse`.

`build_se_epi` assembles every group from an EPI kernel as a separate spin-echo
shot. It supports both single-shot and multishot kernels, derives one feasible
common TE when none is requested, and either triggers every shot or appends a
fixed post-shot delay.

`build_cartesian_bssfp` can split one two-dimensional image across triggered
heartbeats. Every heartbeat repeats the ADC-disabled flip-angle ramp before
acquiring at most the requested number of linear or center-out phase-encoding
lines. `cartesian_line_order` exposes the corresponding line grouping without
constructing a sequence.

`build_pc` constructs a two-point cine PC-GRE acquisition with beat-interleaved
reference and velocity-encoded acquisition windows. It supports physiological
triggering or continuous retrospective acquisition. The requested RR interval
and cardiac-bin count determine the nearest whole-TR bin duration. Incomplete
final line groups are padded with ADC-disabled dummy TRs, while `SET`, `PHS`,
and `LIN` label velocity encoding, cardiac phase, and Cartesian line. Triggered
mode ramps after every trigger; retrospective mode ramps only once at sequence
start. Quadratic RF/receiver spoiling defaults to a 117-degree increment.

`build_pc_encoding_scans` returns separate full-matrix `REF`, `X`, `Y`, and
`Z` velocity-encoding scans. File identity supplies the encoding state, so
these sequences omit `SET` labels. A single-phase scan also omits `PHS`, leaving
only the Cartesian `LIN` labels.

`bssfp_cardiac/generate.jl` generates 128×128 center-out and linear examples
with eight lines per physiological trigger. It writes
`bssfp_cardiac/cardiac.seq` and `bssfp_cardiac/linear.seq`.

`gre_readout_kernel(...; rewind_m0=true)` appends a common-duration decoder so
each Cartesian line returns all three gradient moments to zero. `build_tse`
places those lines on a minimum-time spin-echo train and repeats the train until
the two-dimensional matrix is covered. Linear ordering places `ky=0` on a
selected echo; center-out ordering places it on the first echo. Every train is
either preceded by a physiological input trigger or padded to a requested TR.

`build_grase` places equal odd-length EPI echo groups on a minimum-time train
of spin echoes. The EPI kernel owns line grouping and packet timing; the GRASE
builder owns only the repeated refocusing pulses, crushers, and placement of
each packet center. Effective TE and spin-echo spacing are derived outputs.
The assembled sequence records the three-dimensional imaging box and
`Nx`, `Ny`, and `Nz=1` for KomaMRI raw-data conversion.

## Cartesian indices and labels

Callers pass centered Cartesian indices. Even dimensions use
`-N/2:N/2-1`. Pulseq labels are zero-based array indices:

- GRE: `LIN` and, for 3D, `PAR`;
- EPI: `LIN` plus polarity/segment state;
- SPI: `LIN`, `PAR`, and `SLC` for `(kx, ky, kz)`.

SPI's multi-sample FID remains one ADC event. Pulseq labels apply to that whole
event, not to individual ADC dwell samples.

## Important boundaries

- EPI temporarily represents overlapping blips outside their nominal blocks;
  `epi_readout_kernel` materializes them before returning an M0-refocused
  readout.
- `build_cartesian_bssfp` repeats only the first block and first RF coil of its
  excitation input. Other excitation moments must be included in the kernel's
  fixed-area balance.
- The strict SPI center-out order completes every spherical shell before moving
  outward. Cube-corner shell fragments can therefore produce larger jumps.
- Pulseq files should always be read back and checked with the target scanner.
  `build_spi` compensates for the different scalar block-pulse sample-center
  conventions used by Koma and Pulseq during serialization.
