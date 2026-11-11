using JLD2
using KomaMRIBase
using KomaMRIBase.PulseDesigner: make_delay
import MRISequenceBuildingBlocks as SBB

@testset "Custom RF and SE-EPI composition" begin
    mktempdir() do directory
        rf_file = joinpath(directory, "excitation.jld2")
        waveform = sinpi.(range(0, 1; length=101))
        jldsave(rf_file; B1=2e-6 .* waveform, Gx=1e-3 .* waveform,
            Gy=-0.5e-3 .* waveform, T=1e-3)
        sys = Scanner(limits=HardwareLimits(B1=20e-6, Gmax=30e-3, Smax=180, ADC_Δt=1e-6,
            RF_dead_time=120e-6, RF_ring_down_time=50e-6))
        excitation = SBB.load_rf_excitation(rf_file, sys)
        compensated = SBB.load_rf_excitation(rf_file, sys; moment_compensation=true)
        @test length(excitation) == 1
        @test length(compensated) == 2
        @test compensated[2] ≈ excitation
        @test all(axis -> abs(sum(area, compensated.GR[axis, :])) < 1e-15, 1:3)
        @test isnothing(check_timing(compensated, sys))
        @test isnothing(check_hw_limits(compensated, sys))

        for (fov, resolution, n_shots, n_adc) in (
            ((0.256, 0.256, 0.010), (0.002, 0.002), 5, 131),
            ((0.192, 0.160, 0.008), (0.003, 0.002), 4, 83),
        )
            matrix = round.(Int, fov[1:2] ./ resolution)
            epi = SBB.epi_readout_kernel(fov[1:2], matrix, sys;
                BWpp=500.0, n_shots, spoiler_phase=6π)
            spin_echo_shot = SBB.se_epi(compensated, epi, sys;
                slice_thickness=fov[3], crusher_phase=10π)
            nav_rf = SBB.build_hard_pulse(π/2, sys; B1=5e-6, use=Excitation())
            nav = SBB.build_epi_navigator(nav_rf, epi; n_lines=3, post_delay=0.0)
            TR = 1.0
            seq = Sequence(sys)
            @addblock seq += nav + make_delay(TR - dur(nav))
            for i in -4:n_shots
                spin_echo = spin_echo_shot(max(i, 1))
                if i <= 0
                    spin_echo.epi_readout.ADC .= Ref(ADC(0, 0))
                end
                shot = compensated + spin_echo.delay_before_180 + spin_echo.crushed_180 +
                    spin_echo.delay_before_epi + spin_echo.epi_readout
                @addblock seq += shot
                if i < n_shots
                    @addblock seq += make_delay(TR - dur(shot))
                end
            end
            @test count(is_ADC_on, seq.ADC) == n_adc
            @test count(is_RF_on, seq.RF) == 1 + 2 * (5 + n_shots)
            @test !any(event isa Union{LabelSet,LabelInc} && event.labelstring == "AVG"
                for block in seq.EXT for event in block)
            @test !any(event isa Trigger for block in seq.EXT for event in block)
            @test isnothing(check_timing(seq))
            @test isnothing(check_hw_limits(seq))
            gaps = findall(>(0.5), seq.DUR)
            packets = [seq[first:last] for (first, last) in
                zip([1; gaps .+ 1], [gaps .- 1; length(seq)])]
            @test length(packets) == 1 + 5 + n_shots
            expected_dummy = deepcopy(packets[7])
            expected_dummy.ADC .= Ref(ADC(0, 0))
            @test all(dummy -> dummy ≈ expected_dummy, packets[2:6])
            starts = get_block_start_times(seq)
            rf_blocks = findall(is_RF_on, vec(seq.RF))
            @test all(isapprox.(diff(starts[rf_blocks[2:2:end]]), TR; atol=1e-12, rtol=0))
            @test all(rf_blocks) do block
                rf = seq.RF[1, block]
                delay(rf, sys) >= 120e-6 - 1e-12 &&
                    seq.DUR[block] - (dur(rf, sys) - 50e-6) >= 50e-6 - 1e-12
            end
            path = joinpath(directory, "$(n_shots)shot.seq")
            write_seq(seq, path; sys, verbose=false)
            serialized = read_seq(path; verify_signature=true, verbose=false)
            @test isnothing(check_timing(serialized))
            @test isnothing(check_hw_limits(serialized))
            @test count(is_ADC_on, serialized.ADC) == n_adc
        end
    end
end
