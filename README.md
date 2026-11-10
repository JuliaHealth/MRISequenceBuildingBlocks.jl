<p align="center">
  <img alt="MRISequenceBuildingBlocks.jl" src="./docs/assets/logo-light.png#gh-light-mode-only" width="850">
  <img alt="MRISequenceBuildingBlocks.jl" src="./docs/assets/logo-dark.png#gh-dark-mode-only" width="850">
</p>

<div align="center">

[![Build Status](https://github.com/JuliaHealth/MRISequenceBuildingBlocks.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaHealth/MRISequenceBuildingBlocks.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaHealth/MRISequenceBuildingBlocks.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaHealth/MRISequenceBuildingBlocks.jl)
[![Julia 1.10](https://img.shields.io/badge/julia-v1.10-9558B2?logo=julia)](https://julialang.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Composable MRI pulse-sequence builders with Pulseq support.

</div>

## Installation

```julia-repl
] add MRISequenceBuildingBlocks
```

## Example

```julia
import MRISequenceBuildingBlocks as SBB
using KomaMRI
import KomaMRI.PulseDesigner as PD

sys = Scanner(limits=HardwareLimits(Gmax=30e-3))
FOV = (0.24, 0.24)
matrix = (64, 64)
flip_angle = π / 2
rf_duration = 3e-3
slice_thickness = 5e-3
bandwidth_per_pixel = 500.0
echo_time = 10e-3
phase_encode = 16

excitation = PD.build_sinc_pulse(flip_angle; sys, duration=rf_duration, slice_thickness, use=Excitation())
readout = SBB.gre_readout_kernel(FOV, matrix, sys, bandwidth_per_pixel; rewind_m0=true)
sequence = SBB.build_spin_echo(excitation, readout, sys, phase_encode; slice_thickness, TE=echo_time)

write_seq(sequence, "spin_echo.seq"; sys)
```

Plot with KomaMRI:

```julia
plot_seq(sequence)
plot_kspace(sequence; view_2d=true)
```

<table>
  <tr>
    <td width="63%"><img src="docs/assets/spin_echo_sequence.png" alt="Spin-echo sequence waveforms"></td>
    <td width="37%"><img src="docs/assets/spin_echo_kspace.png" alt="Spin-echo k-space trajectory"></td>
  </tr>
</table>

One spin-echo shot.

Physical inputs use SI units. Unitful quantities are also accepted when
`Unitful` is loaded, for example `slice_thickness = 5u"mm"`.

## Functions

- Spin Echo: `se_epi`, `build_spin_echo`, `build_se_epi`, `build_grase`, `build_tse`
- Excitation: `load_rf_excitation`, `build_hard_pulse`
- Gradient Moments: `build_moment_compensation`
- Labels: `remove_labels!`
- Spoiled Gradient Echo: `sgre_base`
- Balanced SSFP: `build_cartesian_bssfp`, `build_cine_bssfp`
- Phase Contrast: `build_pc`, `build_pc_encoding_scans`
- Single Point Imaging: `build_spi`
- Readout Kernels: `gre_readout_kernel`, `sgre_base`, `bssfp_readout_kernel`, `epi_readout_kernel`
- Navigators: `build_epi_navigator`
- Encoding Order: `cartesian_line_order`, `spi_linear_order`, `spi_center_out_order`


Use Julia help, for example `?build_tse`, for arguments and return values.

## Custom-RF spin-echo EPI

Load a custom excitation and assemble one shot from its gradient prewinder,
RF pulse, crushed slice-selective 180°, and EPI readout:

```julia
using KomaMRIBase
import MRISequenceBuildingBlocks as SBB

sys = Scanner(limits=HardwareLimits(B1=20e-6, Gmax=30e-3, Smax=180,
    ADC_Δt=1e-6, RF_dead_time=120e-6, RF_ring_down_time=50e-6))
fov = (0.256, 0.256, 0.010)
resolution = (0.002, 0.002)
matrix = round.(Int, fov[1:2] ./ resolution)

custom_rf = SBB.load_rf_excitation("pulse.jld2", sys; moment_compensation=true)
epi = SBB.epi_readout_kernel(fov[1:2], matrix, sys;
    BWpp=500.0, n_shots=5, spoiler_phase=6π)
spin_echo_shot = SBB.se_epi(custom_rf, epi, sys;
    slice_thickness=fov[3], crusher_phase=10π)

spin_echo = spin_echo_shot(1)
sequence = custom_rf + spin_echo.delay_before_180 + spin_echo.crushed_180 +
    spin_echo.delay_before_epi + spin_echo.epi_readout
check_timing(sequence, sys)
check_hw_limits(sequence, sys)
```

The JLD2 file supplies `B1` [T], `Gx`/`Gy` and optional `Gz` [T/m], and scalar
duration `T` [s], with uniform waveform knots over `[0,T]`.
`moment_compensation=true` returns two blocks: a prewinder that cancels each
gradient's complete M0, followed by the supplied RF/gradients.

`se_epi` chooses a common minimum TE at the EPI train center and returns a
callable for each shot. Its returned blocks are copies, so disabling ADCs for
a dummy shot does not affect later shots. The EPI kernel includes identical
X/Z spoilers set by the X voxel width when `spoiler_phase` is supplied.

Build NAV separately with `build_hard_pulse` and `build_epi_navigator`, using
the same EPI kernel for matching readout timing and spoilers. EPI/NAV omit AVG
labels. Dummy shots, TR delays, and file export belong in the calling notebook
or script; no external recipe file is required.
