# Migration TODO

This inventory compares `migrated/` with the live project source, including its
uncommitted and untracked Julia files. Gropt code, tests, generated sequence
files, and historical output directories are intentionally excluded.

The items below describe missing behavior and useful ideas. They are not a
requirement to reproduce the original scripts or their APIs exactly.

## Preparation modules

- [ ] Add a frequency-selective fat-saturation preparation based on
  `gaussian_fat_sat_prep` from `fat_sat.jl`:
  - Gaussian 90-degree RF pulse;
  - fat-water offset derived from `B0` and approximately `-3.4 ppm`;
  - B1-limit-aware duration;
  - configurable dephasing across the slice;
  - balanced pre/post z lobes where required by the preparation design.
- [ ] Add support for externally supplied sampled preparation pulses. The
  original code loads `B1`, `Gx`, `Gy`, and `T` from JLD2 files. The reusable API
  should accept waveform data or a pulse constructor; file loading belongs in a
  thin adapter.
- [ ] Add optional x/y moment compensation for custom spatially selective RF
  pulses. This is the useful part of `optimizer_excitation(...; prephase=true)`.
- [ ] Consider a reusable B1-limited nonselective hard-pulse constructor when a
  second sequence needs it. SPI currently performs this design internally, and
  the old spin-echo scripts contain another local implementation; do not extract
  it solely for hypothetical reuse.
- [ ] Add spatial saturation bands or slab-suppression preparation. The
  `logo_gre_single.jl` example creates positive and negative off-slice
  saturation slabs followed by spoilers.
- [ ] Define a simple way for high-level builders to accept an ordered list of
  preparation constructors between a trigger and excitation. This should cover
  fat saturation, custom preparation RF, and spoilers without creating a broad
  preparation framework prematurely.

## EPI navigator and complete EPI acquisitions

- [x] Migrate the EPI navigator from `epi_nav.jl`:
  - independent slice-selective navigator excitation;
  - configurable number of alternating center-k-space lines;
  - readout rewind to zero net x moment;
  - `LIN`, `REV`, `SEG`, and `NAV` labels;
  - navigator timing derived from the same EPI readout design used by the main
    acquisition rather than duplicated formulas.
- [ ] Add a complete gradient-echo EPI builder. `logo_epi.jl` combines trigger,
  fat saturation, navigator, custom excitation, and multishot EPI, while
  `migrated/epi.jl` currently supplies only the readout kernel.
- [ ] Add a complete spin-echo EPI acquisition builder around
  `build_spin_echo` with:
  - all-shot assembly;
  - optional acquisition of selected shots;
  - repeated acquisitions;
  - sequential `LIN` remapping when only selected shots are exported;
  - an explicit delay between navigator and main EPI;
  - coordinated minimum-TE and bandwidth selection that is feasible for every
    acquired shot.
- [x] Define the navigator `AVG` contract: use `AVG=0` on all navigator lines
  except the final line, which uses `AVG=1`; reset it before imaging.

## Spoiled GRE acquisition builders

- [ ] Add a complete Cartesian spoiled-GRE builder on top of `sgre_base`:
  - excitation and ADC phase matched with quadratic RF spoiling;
  - configurable RF-spoiling increment, normally 117 degrees;
  - optional linear ramp preparation with ADC disabled;
  - supplied encoding order or subset;
  - optional trigger/retrigger grouping;
  - minimum or requested TR.
- [ ] Add a full spin-echo GRE acquisition train. `build_spin_echo` can build
  one encoding, but `logo_se_gre_sequence` also handles all phase-encoding
  lines, trigger and preparation insertion, and common echo timing.
- [ ] Consider a generic Cartesian encoding-order helper only when the GRE,
  bSSFP, and flow builders share a concrete need. Do not abstract this merely
  because each currently contains a small loop.

## bSSFP extensions

- [ ] Add radial bSSFP, based on `radial_bSSFP.jl`. The reusable pieces are a
  balanced one-dimensional profile, rotation by each spoke angle, golden-angle
  or supplied angle ordering, ADC-disabled ramp shots, and 180-degree RF/ADC
  phase alternation.
- [ ] Allow the Cartesian bSSFP builder to take a supplied encoding order or
  subset. The original builder also supported random ordering and retriggering
  after a configurable number of lines.
- [ ] Allow optional preparation constructors after each trigger and before the
  ramp. `preppulse_bSSFP.jl` demonstrates a custom RF/gradient preparation and
  spoiler, but its file I/O and top-level generation should not be migrated.
- [ ] If explicit TR control is added, preserve the current symmetry: full
  readout k-space center must remain halfway between RF centers and all gradient
  zeroth moments must remain balanced over the TR.

## Phase contrast and flow encoding

- [ ] Migrate reusable bipolar velocity-encoding design from `pc_gre.jl`:
  - Venc-to-first-moment conversion;
  - minimum-time bipolar timing under gradient and slew limits;
  - common-duration REF/RO/PE/PAR modules;
  - explicit encoding-state representation.
- [ ] Add a triggered two-dimensional cine PC-GRE builder with the useful
  behavior from `make_pcgre`:
  - package or beat interleaving;
  - cardiac phase intervals quantized to complete TR units;
  - configurable phases per trigger;
  - dummy lines where a final phase package is incomplete;
  - `SET` labels for velocity encoding and `PHS` labels for cardiac phase;
  - continuous RF-spoiling state through triggers and packages.
- [ ] Add a retrospective Cartesian 4D-flow spoiled-GRE builder based on
  `make_4dflow_gre`:
  - four REF/RO/PE/PAR encodings;
  - 3D Cartesian sampling;
  - readout partial Fourier;
  - requested or minimum TE and TR;
  - repetitions, ramp shots, RF spoiling, and axis-wise spoilers;
  - externally supplied exact `(ky, kz)` order.
- [ ] Migrate the useful 3D `(ky, kz)` ordering options only with the 4D-flow
  builder: linear, radial/pseudo-spiral center-out, random, exact external order,
  and controlled undersampling.
- [ ] Keep stenotic-flow settings from `stenotic_4dflow_presets.jl` as examples
  or recipes in the future package, not hard-coded library functionality.

## SPI follow-up

- [ ] Improve `spi_center_out_order` at shell transitions. Strict spherical
  shell completion can leave disconnected cube-corner patches and produce large
  inter-shell jumps.
- [ ] Define and verify the even-matrix center convention explicitly.
- [ ] Keep trajectory diagnostics as acceptance criteria: unique coverage,
  monotonic shell membership, maximum local jump, cumulative mean by axis and
  norm, and final mean at k-space center.
- [ ] Consider additional Cartesian orderings only when they have a clear scan
  use case, for example randomized ordering or an externally supplied order with
  validation. `build_spi` already accepts arbitrary valid point lists.

## Other non-Cartesian ideas

- [ ] Add a reusable radial GRE readout and full radial spoiled-GRE builder.
  There is no clean implementation in the source tree, but radial bSSFP already
  demonstrates the rotation and spoke-order concepts.
- [ ] Explore UTE as a separate readout family. Shared minimum-time gradient
  timing may be reusable, but UTE needs explicit RF-to-readout dead-time and
  center-sample semantics rather than being forced through the Cartesian GRE
  kernel.
- [ ] Consider stack-of-stars or 3D radial acquisition only after the 2D radial
  kernel contract is settled.

## Acquisition orchestration and metadata

- [ ] Decide which orchestration belongs in the library: triggers, retrigger
  groups, ramp/dummy shots, repetitions, shot subsets, and preparation insertion
  are recurring across the old bSSFP, GRE, EPI, and PC-GRE scripts.
- [ ] Establish a minimal label policy for `LIN`, `PAR`, `SLC`, `SEG`, `REV`,
  `NAV`, `SET`, `PHS`, `REP`, and `AVG`, including which builder owns each label.
- [ ] Establish a minimal Pulseq-definition policy. Keep timing and construction
  state in constructors or returned kernel data; write only metadata needed by
  scanners, reconstruction, or artifact identification to `seq.DEF`.
- [ ] Keep `write_seq`, filename selection, scanner presets, batch generation,
  and artifact copying outside core constructors. The `sequence_share/` scripts
  are useful recipe examples, not package APIs.

## Deferred API decisions

- [ ] Define the module, exports, and exact KomaMRI dependency boundary when
  these files move into their destination package. Do not add package scaffolding
  to this staging directory.
- [ ] Consider replacing closure-based kernels with callable structs only if a
  real need emerges for kernel inspection, serialization, or multiple dispatch.
  The current closures are otherwise simpler.
- [ ] Normalize public naming when these files become a package. In particular,
  decide whether `sgre_base` should become `sgre_readout_kernel` and define a
  consistent distinction between event makers, one-block builders, readout
  kernels, and full acquisition builders.
- [ ] Decide whether custom RF and preparation waveform adapters live in this
  package or a small I/O companion. Core builders should not depend directly on
  JLD2.

## Reviewed but intentionally not migrated

- The original `make_trigger` wrapper is superseded by KomaMRI's
  `make_trigger`/`build_trigger` functions.
- The original `gre_base`, `bssfp_base`, and `spi_base` are superseded by the
  migrated kernels and builders; only missing acquisition-level behavior is
  listed above.
- `preprocess_seq` is represented by `materialize_gradients`; its old demo and
  stress code do not belong in the library.
- `build_phantom_static.jl` and `build_phantom_motion.jl` are useful Koma phantom
  conversion examples, but they belong in simulation tooling or documentation,
  not this sequence-construction layer.
- Gropt examples and `pc_bssfp_gropt2.jl` remain out of scope.
- Top-level scripts such as `logo_gre_single.jl`, `logo_se_epi_single.jl`,
  `logo_se_readout.jl`, and `stenotic_4dflow_presets.jl` are sources of sequence
  ideas and acceptance cases, not APIs to copy into `migrated/`.
