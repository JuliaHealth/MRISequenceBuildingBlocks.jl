# Migrated sequence builders

This directory is a staging area for high-level sequence construction code. It
is intentionally not a Julia package yet. Numeric design calculations use SI
units internally; Unitful quantities are added only when calling KomaMRI pulse
constructors.

## Suggested include order

```julia
include("gradient_design.jl")
include("utils.jl")
include("preprocess.jl")
include("excitation.jl")
include("gre.jl")
include("epi.jl")
include("bSSFP.jl")
include("spinEcho.jl")
include("spi.jl")
```

`gre.jl`, `epi.jl`, `bSSFP.jl`, and `spi.jl` use `_lobe_timing` from
`gradient_design.jl`. `bSSFP.jl` also uses the GRE kernel and `rf_excitation`.
The public EPI kernel uses `materialize_gradients` from `preprocess.jl`.

## Public surfaces

| File | Main functions |
| --- | --- |
| `excitation.jl` | `build_centered_sinc_pulse`, `slice_selective_sinc`, `rf_sinc` |
| `gre.jl` | `gre_readout_kernel`, `sgre_base` |
| `epi.jl` | `epi_readout_kernel` |
| `bSSFP.jl` | `bssfp_readout_kernel`, `build_cartesian_bssfp` |
| `spinEcho.jl` | `build_refocusing_block`, `build_spin_echo` |
| `spi.jl` | `spi_linear_order`, `spi_center_out_order`, `build_spi` |
| `preprocess.jl` | `materialize_gradients` |
| `utils.jl` | `rf_excitation` |

GRE, EPI, and bSSFP use a small kernel contract: `readout(args...)` builds one
concrete readout and `center_time(args...)` reports the k-space-center time
relative to its beginning. `build_spin_echo` consumes this contract without
branching on the readout family.

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
  `epi_readout_kernel` materializes them before returning a readout.
- `build_cartesian_bssfp` repeats only the first block and first RF coil of its
  excitation input. Other excitation moments must be included in the kernel's
  fixed-area balance.
- The strict SPI center-out order completes every spherical shell before moving
  outward. Cube-corner shell fragments can therefore produce larger jumps.
- Pulseq files should always be read back and checked with the target scanner.
  `build_spi` compensates for the different scalar block-pulse sample-center
  conventions used by Koma and Pulseq during serialization.
