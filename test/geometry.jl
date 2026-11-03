@testset "Cartesian ADC samples and encoding labels" begin
    sys = Scanner()
    readout_fov = 0.24
    phase_fov = 0.24
    readout_samples = 16
    phase_lines = 8
    bandwidth_per_pixel = 500.0
    readout_axis = 1
    phase_axis = 2

    FOV = (readout_fov, phase_fov)
    matrix = (readout_samples, phase_lines)
    readout_indices = (-(readout_samples ÷ 2)):(readout_samples - 1 - readout_samples ÷ 2)
    phase_indices = (-(phase_lines ÷ 2)):(phase_lines - 1 - phase_lines ÷ 2)
    kx = collect(readout_indices) ./ readout_fov

    gre = gre_readout_kernel(FOV, matrix, sys, bandwidth_per_pixel)
    for line in phase_indices
        seq = gre.readout(line)
        _, k = get_kspace(seq)
        @test k[:, readout_axis] ≈ kx
        @test k[:, phase_axis] ≈ fill(line / phase_fov, readout_samples)
        context = adc_label_context(seq)
        @test [context.labels[b].LIN for b in context.adc_blocks] == [line - first(phase_indices)]
    end

    epi_shots = 2
    epi = epi_readout_kernel(FOV, matrix, sys;
        BWpp=bandwidth_per_pixel, n_shots=epi_shots)
    for shot in 1:epi_shots
        seq = epi.readout(shot)
        _, k = get_kspace(seq)
        lines = collect(phase_indices[shot:epi_shots:end])
        expected_kx = reduce(vcat, [isodd(echo) ? kx : reverse(kx) for echo in eachindex(lines)])
        expected_ky = repeat(lines ./ phase_fov; inner=readout_samples)
        @test k[:, readout_axis] ≈ expected_kx
        @test k[:, phase_axis] ≈ expected_ky
        context = adc_label_context(seq)
        @test [context.labels[b].LIN for b in context.adc_blocks] == lines .- first(phase_indices)
        @test [context.labels[b].REV for b in context.adc_blocks] == Int.(iseven.(eachindex(lines)))
    end
end

@testset "bSSFP balances all three gradients including slice selection" begin
    sys = Scanner()
    readout_fov = 0.24
    phase_fov = 0.24
    readout_samples = 16
    phase_lines = 8
    bandwidth_per_pixel = 500.0
    flip_angle = π / 6
    rf_duration = 1e-3
    slice_thickness = 5e-3
    apodization = 0.46
    slice_axis = 3
    excitation_block = 1
    rephaser_block = 2

    FOV = (readout_fov, phase_fov)
    matrix = (readout_samples, phase_lines)
    phase_indices = (-(phase_lines ÷ 2)):(phase_lines - 1 - phase_lines ÷ 2)

    excitation = build_sinc_pulse(flip_angle; sys, use=Excitation(),
        duration=rf_duration, slice_thickness, apodization)
    slice_rephaser = area(excitation.GR[slice_axis, rephaser_block])
    zero_area = zero(slice_rephaser)
    fixed_area = (zero_area, zero_area, slice_rephaser)
    balanced = bssfp_readout_kernel(FOV, matrix, sys, bandwidth_per_pixel; fixed_area)
    for line in phase_indices
        cycle = excitation[excitation_block] + balanced.readout(line)
        for gradients in eachrow(cycle.GR)
            areas = area.(gradients)
            positive_area = sum(max.(areas, zero_area))
            negative_area = -sum(min.(areas, zero_area))
            @test positive_area ≈ negative_area
        end
    end
end

@testset "Pulseq preserves EPI ADC samples, labels, and k-space" begin
    sys = Scanner()
    readout_fov = 0.24
    phase_fov = 0.24
    readout_samples = 16
    phase_lines = 8
    bandwidth_per_pixel = 500.0
    epi_shots = 2

    FOV = (readout_fov, phase_fov)
    matrix = (readout_samples, phase_lines)

    epi = epi_readout_kernel(FOV, matrix, sys;
        BWpp=bandwidth_per_pixel, n_shots=epi_shots)
    original = reduce(+, epi.readout(shot) for shot in 1:epi_shots)
    expected_adc_samples = prod(matrix)
    mktempdir() do dir
        path = joinpath(dir, "epi.seq")
        write_seq(original, path; sys, verbose=false)
        restored = read_seq(path; verify_signature=true, verbose=false)
        @test sum(adc.N for adc in restored.ADC) == expected_adc_samples
        before, after = adc_label_context(original), adc_label_context(restored)
        for name in (:LIN, :REV, :SEG)
            @test [getproperty(after.labels[b], name) for b in after.adc_blocks] ==
                [getproperty(before.labels[b], name) for b in before.adc_blocks]
        end
        @test last(get_kspace(restored)) ≈ last(get_kspace(original))
    end
end
