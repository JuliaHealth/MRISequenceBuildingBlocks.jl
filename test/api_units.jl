using Unitful

@testset "Cartesian kernels accept numeric and Unitful inputs" begin
    sys = Scanner()
    FOV = (0.24, 0.24)
    FOV_units = (24u"cm", 24u"cm")
    matrix = (16, 8)
    bandwidth = 500.0
    bandwidth_units = 500u"Hz"
    fixed_area = (0.0, 0.0, 1e-6)
    fixed_area_units = (0u"T*s/m", 0u"T*s/m", 1e-6u"T*s/m")
    spoil_area = (0.0, 0.0, 1e-5)
    spoil_area_units = (0u"T*s/m", 0u"T*s/m", 1e-5u"T*s/m")
    phase_dimension = 2
    phase_lines = matrix[phase_dimension]
    phase_indices = (-(phase_lines ÷ 2)):(phase_lines - 1 - phase_lines ÷ 2)

    gre = gre_readout_kernel(FOV, matrix, sys, bandwidth; fixed_area)
    gre_units = gre_readout_kernel(FOV_units, matrix, sys, bandwidth_units;
        fixed_area=fixed_area_units)
    spoiled = sgre_base(FOV, matrix, sys, bandwidth; fixed_area, spoil_area)
    spoiled_units = sgre_base(FOV_units, matrix, sys, bandwidth_units;
        fixed_area=fixed_area_units, spoil_area=spoil_area_units)
    balanced = bssfp_readout_kernel(FOV, matrix, sys, bandwidth; fixed_area)
    balanced_units = bssfp_readout_kernel(FOV_units, matrix, sys, bandwidth_units;
        fixed_area=fixed_area_units)
    for line in phase_indices
        @test gre.readout(line) ≈ gre_units.readout(line)
        @test gre.center_time(line) ≈ gre_units.center_time(line)
        @test spoiled(line) ≈ spoiled_units(line)
        @test balanced.readout(line) ≈ balanced_units.readout(line)
    end
end

@testset "EPI kernel and navigator accept numeric and Unitful inputs" begin
    sys = Scanner()
    FOV = (0.24, 0.24)
    FOV_units = (24u"cm", 24u"cm")
    matrix = (16, 8)
    bandwidth = 500.0
    bandwidth_units = 500u"Hz"
    shots = 2
    flip_angle = π / 6
    flip_angle_units = 30u"°"
    rf_duration = 3e-3
    rf_duration_units = 3u"ms"
    slice_thickness = 5e-3
    slice_thickness_units = 5u"mm"
    navigator_lines = 3
    post_delay = 1e-3
    post_delay_units = 1u"ms"

    epi = epi_readout_kernel(FOV, matrix, sys; BWpp=bandwidth, n_shots=shots)
    epi_units = epi_readout_kernel(FOV_units, matrix, sys;
        BWpp=bandwidth_units, n_shots=shots)
    for shot in 1:shots
        @test epi.readout(shot) ≈ epi_units.readout(shot)
        @test epi.center_time(shot) ≈ epi_units.center_time(shot)
    end
    excitation = build_sinc_pulse(flip_angle; sys, use=Excitation(),
        duration=rf_duration, slice_thickness)
    excitation_units = build_sinc_pulse(flip_angle_units; sys, use=Excitation(),
        duration=rf_duration_units, slice_thickness=slice_thickness_units)
    navigator = build_epi_navigator(excitation, epi; n_lines=navigator_lines, post_delay)
    navigator_units = build_epi_navigator(excitation_units, epi_units;
        n_lines=navigator_lines, post_delay=post_delay_units)
    @test navigator ≈ navigator_units
end

@testset "bSSFP acquisitions accept numeric and Unitful inputs" begin
    sys = Scanner()
    FOV = (0.24, 0.24)
    FOV_units = (24u"cm", 24u"cm")
    volume_FOV = (0.24, 0.24, 5e-3)
    volume_FOV_units = (24u"cm", 24u"cm", 5u"mm")
    matrix = (16, 8)
    bandwidth = 500.0
    bandwidth_units = 500u"Hz"
    flip_angle = π / 6
    flip_angle_units = 30u"°"
    rf_duration = 3e-3
    rf_duration_units = 3u"ms"
    slice_thickness = 5e-3
    slice_thickness_units = 5u"mm"
    slice_axis = 3
    rephaser_block = 2
    ramp_shots = 1
    cardiac_bins = 2
    RR = 0.08
    RR_units = 80u"ms"
    post_trigger_delay = 0.0
    post_trigger_delay_units = 0u"ms"
    steady_state_duration = 0.0
    steady_state_duration_units = 0u"ms"
    slice_normal = (1.0, 1.0, 1.0)
    slice_shift = 1e-3
    slice_shift_units = 1u"mm"

    excitation = build_sinc_pulse(flip_angle; sys, use=Excitation(),
        duration=rf_duration, slice_thickness)
    excitation_units = build_sinc_pulse(flip_angle_units; sys, use=Excitation(),
        duration=rf_duration_units, slice_thickness=slice_thickness_units)
    rephaser_area = area(excitation.GR[slice_axis, rephaser_block])
    zero_area = zero(rephaser_area)
    fixed_area = (zero_area, zero_area, rephaser_area)
    balanced = bssfp_readout_kernel(FOV, matrix, sys, bandwidth; fixed_area)
    balanced_units = bssfp_readout_kernel(FOV_units, matrix, sys, bandwidth_units; fixed_area)

    cartesian = build_cartesian_bssfp(excitation, balanced, sys;
        trigger=nothing, n_ramp_shots=ramp_shots, post_trigger_delay)
    cartesian_units = build_cartesian_bssfp(excitation_units, balanced_units, sys;
        trigger=nothing, n_ramp_shots=ramp_shots, post_trigger_delay=post_trigger_delay_units)
    @test cartesian ≈ cartesian_units

    cine = build_cine_bssfp(excitation, balanced, sys;
        cardiac_bins, RR, n_ramp_shots=ramp_shots, steady_state_duration, post_trigger_delay)
    cine_units = build_cine_bssfp(excitation_units, balanced_units, sys;
        cardiac_bins, RR=RR_units, n_ramp_shots=ramp_shots,
        steady_state_duration=steady_state_duration_units, post_trigger_delay=post_trigger_delay_units)
    @test cine ≈ cine_units

    oriented = build_cine_bssfp(volume_FOV, matrix, sys, bandwidth;
        excitation, slice_normal, slice_shift, cardiac_bins, RR,
        n_ramp_shots=ramp_shots, steady_state_duration)
    oriented_units = build_cine_bssfp(volume_FOV_units, matrix, sys, bandwidth_units;
        excitation=excitation_units, slice_normal, slice_shift=slice_shift_units,
        cardiac_bins, RR=RR_units, n_ramp_shots=ramp_shots,
        steady_state_duration=steady_state_duration_units)
    @test oriented ≈ oriented_units
end

@testset "Spin echo API accepts numeric and Unitful inputs" begin
    sys = Scanner()
    FOV = (0.24, 0.24)
    FOV_units = (24u"cm", 24u"cm")
    matrix = (16, 8)
    grase_matrix = (16, 9)
    bandwidth = 500.0
    bandwidth_units = 500u"Hz"
    shots = 2
    grase_shots = 3
    flip_angle = π / 6
    flip_angle_units = 30u"°"
    rf_duration = 3e-3
    rf_duration_units = 3u"ms"
    slice_thickness = 5e-3
    slice_thickness_units = 5u"mm"
    refocusing_bandwidth = 500.0
    refocusing_bandwidth_units = 500u"Hz"
    refocusing_phase = π / 2
    refocusing_phase_units = 90u"°"
    crusher_phase = 4π
    crusher_phase_units = 720u"°"
    echo_time = 30e-3
    echo_time_units = 30u"ms"
    repetition_time = 0.2
    repetition_time_units = 200u"ms"
    post_shot_delay = 0.01
    post_shot_delay_units = 10u"ms"
    center_line = 0
    echo_train_length = 4

    excitation = build_sinc_pulse(flip_angle; sys, use=Excitation(),
        duration=rf_duration, slice_thickness)
    excitation_units = build_sinc_pulse(flip_angle_units; sys, use=Excitation(),
        duration=rf_duration_units, slice_thickness=slice_thickness_units)
    gre = gre_readout_kernel(FOV, matrix, sys, bandwidth; rewind_m0=true)
    gre_units = gre_readout_kernel(FOV_units, matrix, sys, bandwidth_units; rewind_m0=true)
    epi = epi_readout_kernel(FOV, matrix, sys; BWpp=bandwidth, n_shots=shots)
    epi_units = epi_readout_kernel(FOV_units, matrix, sys; BWpp=bandwidth_units, n_shots=shots)
    grase_epi = epi_readout_kernel(FOV, grase_matrix, sys;
        BWpp=bandwidth, n_shots=grase_shots, shot_partition=:contiguous)
    grase_epi_units = epi_readout_kernel(FOV_units, grase_matrix, sys;
        BWpp=bandwidth_units, n_shots=grase_shots, shot_partition=:contiguous)

    spin_echo = build_spin_echo(excitation, gre, sys, center_line;
        slice_thickness, TE=echo_time, crusher_phase, refocusing_phase, refocusing_bandwidth)
    spin_echo_units = build_spin_echo(excitation_units, gre_units, sys, center_line;
        slice_thickness=slice_thickness_units, TE=echo_time_units,
        crusher_phase=crusher_phase_units, refocusing_phase=refocusing_phase_units,
        refocusing_bandwidth=refocusing_bandwidth_units)
    @test spin_echo ≈ spin_echo_units

    se_epi = build_se_epi(excitation, epi, sys;
        slice_thickness, post_shot_delay, crusher_phase, refocusing_phase, refocusing_bandwidth)
    se_epi_units = build_se_epi(excitation_units, epi_units, sys;
        slice_thickness=slice_thickness_units, post_shot_delay=post_shot_delay_units,
        crusher_phase=crusher_phase_units, refocusing_phase=refocusing_phase_units,
        refocusing_bandwidth=refocusing_bandwidth_units)
    @test se_epi ≈ se_epi_units

    grase = build_grase(excitation, grase_epi, sys;
        slice_thickness, crusher_phase, refocusing_phase, refocusing_bandwidth)
    grase_units = build_grase(excitation_units, grase_epi_units, sys;
        slice_thickness=slice_thickness_units, crusher_phase=crusher_phase_units,
        refocusing_phase=refocusing_phase_units, refocusing_bandwidth=refocusing_bandwidth_units)
    @test grase ≈ grase_units

    tse = build_tse(excitation, gre, sys; slice_thickness, echo_train_length,
        TR=repetition_time, crusher_phase, refocusing_phase, refocusing_bandwidth)
    tse_units = build_tse(excitation_units, gre_units, sys;
        slice_thickness=slice_thickness_units, echo_train_length, TR=repetition_time_units,
        crusher_phase=crusher_phase_units, refocusing_phase=refocusing_phase_units,
        refocusing_bandwidth=refocusing_bandwidth_units)
    @test tse ≈ tse_units
end

@testset "Phase contrast API accepts numeric and Unitful inputs" begin
    sys = Scanner()
    FOV = (0.24, 0.24, 5e-3)
    FOV_units = (24u"cm", 24u"cm", 5u"mm")
    matrix = (16, 8)
    velocity = 1.0
    velocity_units = 1u"m/s"
    bandwidth = 500.0
    bandwidth_units = 500u"Hz"
    flip_angle = π / 6
    flip_angle_units = 30u"°"
    rf_bandwidth = 500.0
    rf_bandwidth_units = 500u"Hz"
    spoil_phase = 4π
    spoil_phase_units = 720u"°"
    rf_spoil_increment_degrees = 117
    rf_spoil_increment = deg2rad(rf_spoil_increment_degrees)
    rf_spoil_increment_units = rf_spoil_increment_degrees * u"°"
    steady_state_duration = 0.0
    steady_state_duration_units = 0u"ms"
    cardiac_bins = 1
    ramp_shots = 0

    pc = build_pc(FOV, matrix, sys; venc=velocity, cardiac_bins, BWpp=bandwidth,
        trigger_channel=nothing, n_ramp_shots=ramp_shots, steady_state_duration,
        flip_angle, rf_bandwidth, spoil_phase, rf_spoil_increment)
    pc_units = build_pc(FOV_units, matrix, sys;
        venc=velocity_units, cardiac_bins, BWpp=bandwidth_units,
        trigger_channel=nothing, n_ramp_shots=ramp_shots,
        steady_state_duration=steady_state_duration_units, flip_angle=flip_angle_units,
        rf_bandwidth=rf_bandwidth_units, spoil_phase=spoil_phase_units,
        rf_spoil_increment=rf_spoil_increment_units)
    @test pc ≈ pc_units
    @test pc_units.DEF["Venc"] ≈ velocity
    @test pc_units.DEF["RfSpoilIncrement"] ≈ rf_spoil_increment

    scans = build_pc_encoding_scans(FOV, matrix, sys;
        venc=velocity, cardiac_bins, BWpp=bandwidth, trigger_channel=nothing,
        n_ramp_shots=ramp_shots, steady_state_duration, flip_angle)
    scans_units = build_pc_encoding_scans(FOV_units, matrix, sys;
        venc=velocity_units, cardiac_bins, BWpp=bandwidth_units, trigger_channel=nothing,
        n_ramp_shots=ramp_shots, steady_state_duration=steady_state_duration_units,
        flip_angle=flip_angle_units)
    for encoding in keys(scans)
        @test scans[encoding] ≈ scans_units[encoding]
    end
end

@testset "Ordering API uses integer indices" begin
    phase_indices = -4:3
    lines_per_shot = 3
    expected_center_out = [0, -1, 1, -2, 2, -3, 3, -4]
    groups = cartesian_line_order(phase_indices, lines_per_shot; view_order=:center_out)
    @test reduce(vcat, groups) == expected_center_out
    @test all(group -> length(group) <= lines_per_shot, groups)

    matrix = (4, 4, 4)
    shell_thickness = 1.0
    spiral_spacing = 1.0
    expected_points = prod(matrix)
    linear = spi_linear_order(matrix)
    center_out = spi_center_out_order(matrix; shell_thickness, spiral_spacing)
    @test length(unique(linear)) == expected_points
    @test Set(center_out) == Set(linear)
    @test length(center_out) == expected_points
    shell_indices = [floor(Int, sqrt(sum(abs2, point)) / shell_thickness) for point in center_out]
    @test issorted(shell_indices)
end

@testset "SPI and Pulseq API preserve numeric/Unitful parity" begin
    sys = Scanner()
    FOV = (0.24, 0.24, 0.24)
    FOV_units = (24u"cm", 24u"cm", 24u"cm")
    matrix = (4, 4, 4)
    flip_angle_degrees = 13
    flip_angle = deg2rad(flip_angle_degrees)
    flip_angle_units = flip_angle_degrees * u"°"
    adc_samples = 48
    adc_duration = 2e-3
    adc_duration_units = 2u"ms"
    repetition_time = 0.01
    repetition_time_units = 10u"ms"
    spoil_phase = 4π
    spoil_phase_units = 720u"°"
    order = spi_linear_order(matrix)

    spi = build_spi(FOV, matrix, order, sys;
        flip_angle, adc_samples, adc_duration, TR=repetition_time, spoil_phase)
    spi_units = build_spi(FOV_units, matrix, order, sys;
        flip_angle=flip_angle_units, adc_samples, adc_duration=adc_duration_units,
        TR=repetition_time_units, spoil_phase=spoil_phase_units)
    @test spi ≈ spi_units
    expected_adc_samples = prod(matrix) * adc_samples
    mktempdir() do dir
        numeric_path = joinpath(dir, "numeric.seq")
        unitful_path = joinpath(dir, "unitful.seq")
        write_seq(spi, numeric_path; sys, verbose=false)
        write_seq(spi_units, unitful_path; sys, verbose=false)
        numeric_readback = read_seq(numeric_path; verify_signature=true, verbose=false)
        unitful_readback = read_seq(unitful_path; verify_signature=true, verbose=false)
        @test numeric_readback ≈ unitful_readback
        @test sum(adc.N for adc in unitful_readback.ADC) == expected_adc_samples
    end
end
