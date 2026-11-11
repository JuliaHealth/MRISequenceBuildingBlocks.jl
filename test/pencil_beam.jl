using KomaMRIBase
using Unitful
import MRISequenceBuildingBlocks as SBB

function pencil_beam_profile(seq, positions)
    excitation = seq[2]
    rf = excitation.RF[1, 1]
    rf_times = times(rf)[2:end-1]
    grid = sort!(unique([rf_times; times(excitation.GR[1, 1]); times(excitation.GR[2, 1])]))
    gx, gy, _ = KomaMRIBase.get_grads(excitation, grid)
    kx = [0.0; cumsum(diff(grid) .* (gx[1:end-1] + gx[2:end]) .* (γ / 2))]
    ky = [0.0; cumsum(diff(grid) .* (gy[1:end-1] + gy[2:end]) .* (γ / 2))]
    indices = searchsortedfirst.(Ref(grid), rf_times)
    kx, ky = kx[end] .- kx[indices], ky[end] .- ky[indices]
    weights = ComplexF64.(cis(rf.ϕ) .* rf.A) .* dwell(rf)
    weights[begin] /= 2
    weights[end] /= 2
    return [sum(weights .* cis.(-2π .* (kx .* x + ky .* y))) / sum(weights)
        for (x, y) in positions]
end

@testset "Analytic pencil-beam excitation" begin
    sys = Scanner(limits=HardwareLimits(B1=20e-6, Gmax=30e-3, Smax=150.0,
        RF_dead_time=120e-6, RF_ring_down_time=50e-6))
    seq = SBB.build_pencil_beam(π/2, 0.080, 0.256, 40, sys)
    @test length(seq) == 2
    @test !any(is_RF_on, seq[1].RF)
    @test count(is_RF_on, seq.RF) == 1
    @test !any(is_ADC_on, seq.ADC)
    @test !any(is_GR_on, seq.GR[3, :])
    @test is_GR_on(seq.GR[1, 2]) && is_GR_on(seq.GR[2, 2])
    @test all(axis -> abs(sum(area, seq.GR[axis, :])) < 1e-14, 1:3)
    @test get_flip_angles(seq)[2] ≈ 90.0
    @test minimum(seq.RF[1, 2].A) < 0 < maximum(seq.RF[1, 2].A)
    @test maximum(abs, seq.RF[1, 2].A) <= sys.limits.B1
    @test delay(seq.RF[1, 2], sys) >= sys.limits.RF_dead_time
    @test seq.DUR[2] >= dur(seq.RF[1, 2], sys) - 1e-12
    @test isnothing(check_timing(seq, sys))
    @test isnothing(check_hw_limits(seq, sys))

    units = SBB.build_pencil_beam(90u"°", 80u"mm", 256u"mm", 40, sys)
    @test units ≈ seq
    capped = SBB.build_pencil_beam(π/2, 0.080, 0.256, 40, sys; B1=12u"μT")
    @test maximum(abs, capped.RF[1, 2].A) <= 12e-6
    @test dur(capped) < 15e-3
    @test get_flip_angles(capped)[2] ≈ 90.0
    @test SBB.build_pencil_beam(π/2, 0.080, 0.256, 40, sys; B1=30e-6) ≈ seq
    reversed = SBB.build_pencil_beam(-π/2, 0.080, 0.256, 40, sys; B1=12e-6)
    @test area(reversed.RF[1, 2]) ≈ -area(capped.RF[1, 2])
    @test maximum(abs, reversed.RF[1, 2].A) <= 12e-6
    low_B1_sys = Scanner(limits=HardwareLimits(B1=1e-6, Gmax=30e-3, Smax=150.0,
        RF_dead_time=120e-6, RF_ring_down_time=50e-6))
    low_B1 = SBB.build_pencil_beam(π/2, 0.080, 0.256, 40, low_B1_sys)
    @test dur(low_B1) > dur(seq)
    @test maximum(abs, low_B1.RF[1, 2].A) <= low_B1_sys.limits.B1
    @test get_flip_angles(low_B1)[2] ≈ 90.0
    @test isnothing(check_timing(low_B1, low_B1_sys))
    @test isnothing(check_hw_limits(low_B1, low_B1_sys))
    positions = [(x, 0.0) for x in (0.0, 0.020, 0.040, 0.060)]
    @test all(isapprox.(abs.(pencil_beam_profile(low_B1, positions)),
        abs.(pencil_beam_profile(seq, positions)); atol=0.02))

    # Independently integrate the small-tip response from the emitted events.
    for diameter in (0.060, 0.080)
        small_tip = SBB.build_pencil_beam(5u"°", diameter, 0.256, 64, sys)
        radius = diameter / 2
        inside = [(radius * fraction, 0.0) for fraction in (0.0, 0.25, 0.5, 0.7)]
        outside = [(radius * fraction, 0.0) for fraction in (1.4, 1.6, 2.0)]
        diagonal = [(x / sqrt(2), x / sqrt(2)) for (x, _) in inside]
        @test all(value -> isapprox(abs(value), 1.0; atol=0.05),
            pencil_beam_profile(small_tip, inside))
        @test maximum(abs, pencil_beam_profile(small_tip, outside)) < 0.04
        @test abs(pencil_beam_profile(small_tip, [(radius, 0.0)])[1]) ≈ 0.5 atol=0.08
        @test abs.(pencil_beam_profile(small_tip, diagonal)) ≈
            abs.(pencil_beam_profile(small_tip, inside)) atol=0.03
    end

    mktempdir() do directory
        path = joinpath(directory, "pencil_beam.seq")
        write_seq(capped, path; sys, verbose=false)
        serialized = read_seq(path; verify_signature=true, verbose=false)
        @test length(serialized) == 2
        @test !any(is_ADC_on, serialized.ADC)
        @test get_flip_angles(serialized)[2] ≈ 90.0
        @test maximum(abs, serialized.RF[1, 2].A) <= 12e-6
        @test isnothing(check_timing(serialized, sys))
        @test isnothing(check_hw_limits(serialized, sys))
        @test all(axis -> abs(sum(area, serialized.GR[axis, :])) < 1e-10, 1:3)
    end

    epi = SBB.epi_readout_kernel((0.032, 0.032), (16, 8), sys; BWpp=500.0)
    spin_echo = SBB.se_epi(seq, epi, sys; slice_thickness=0.010, crusher_phase=10π)(1)
    acquisition = seq + spin_echo.delay_before_180 + spin_echo.crushed_180 +
        spin_echo.delay_before_epi + spin_echo.epi_readout
    @test isnothing(check_timing(acquisition, sys))
    @test isnothing(check_hw_limits(acquisition, sys))

    @test_throws ErrorException SBB.build_pencil_beam(0.0, 0.080, 0.256, 40, sys)
    @test_throws ErrorException SBB.build_pencil_beam(π/2, -0.080, 0.256, 40, sys)
    @test_throws ErrorException SBB.build_pencil_beam(π/2, 0.300, 0.256, 40, sys)
    @test_throws ErrorException SBB.build_pencil_beam(π/2, 0.080, 0.256, 1, sys)
    @test_throws ErrorException SBB.build_pencil_beam(π/2, 0.080, 0.256, 40, sys; apodization=1.1)
    @test_throws ErrorException SBB.build_pencil_beam(π/2, 0.080, 0.256, 40, sys; B1=0.0)
end
