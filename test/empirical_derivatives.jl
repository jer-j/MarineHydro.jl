using MarineHydro
using Test

@testset "Empirical derivative formulas" begin
    # KVLCC2 full-scale particulars, as in validation/kvlcc2_maneuvering.
    kvlcc2 = HullParticulars(; length_pp = 320.0, beam = 58.0, draft = 20.8,
        block_coefficient = 0.8098)

    @testset "Reproduce the published KVLCC2 column" begin
        # Chame et al. (2025), Table 9, KVLCC2 column. These are the numbers the
        # notebook compares against, so if a formula is mistyped it shows here
        # rather than silently shifting a plotted curve.
        published = Dict(:low_aspect_ratio => 0.204, :clarke => 0.389,
            :inoue => 0.410, :norrbin => 0.357, :smitt => 0.325,
            :ho_young_lee => 0.308)
        for (method, expected) in published
            value = empirical_sway_derivatives(kvlcc2; method).Y_v
            @test value≈expected rtol=0.03
        end
    end

    @testset "The low-aspect-ratio value is pi T / L" begin
        base = empirical_sway_derivatives(kvlcc2; method = :low_aspect_ratio)
        @test base.Y_v≈π * kvlcc2.draft / kvlcc2.length_pp rtol=1e-12
        @test base.N_v≈base.Y_v / 2 rtol=1e-12
        # Every regression is a correction on that value, so all of them should
        # sit within a factor of a few of it.
        for method in empirical_methods()
            value = empirical_sway_derivatives(kvlcc2; method).Y_v
            @test 0.5 * base.Y_v < value < 3 * base.Y_v
        end
    end

    @testset "Clarke's full set has the expected signs" begin
        clarke = clarke_full_derivatives(kvlcc2)
        # Sway force opposes sway velocity; yaw damping opposes yaw rate.
        @test clarke.Y_v < 0
        @test clarke.N_r < 0
        @test clarke.Y_vdot < 0
        @test clarke.N_rdot < 0
        # And the velocity derivative agrees with the tabulated form, up to the
        # sign convention the two are quoted in.
        @test abs(clarke.Y_v)≈empirical_sway_derivatives(kvlcc2;
            method = :clarke).Y_v * kvlcc2.draft / kvlcc2.length_pp rtol=1e-10
    end

    @testset "Hirano and Takashina" begin
        hirano = hirano_takashina_derivatives(kvlcc2)
        aspect = 2 * kvlcc2.draft / kvlcc2.length_pp
        @test hirano.N_v≈-aspect rtol=1e-12
        @test hirano.Y_r≈0.25π * aspect rtol=1e-12
        @test hirano.N_r < 0

        # On a hull as full as KVLCC2 the +1.4 C_B B/L correction very nearly
        # cancels the -0.5 pi Lambda term, and Y_v comes out marginally
        # POSITIVE, which no real hull does. Recorded rather than asserted away:
        # it marks the formula as outside its useful range at this block
        # coefficient, and is the sort of thing a notebook plotting it without
        # comment would hide.
        @test hirano.Y_v > -0.01
        @test abs(hirano.Y_v) < 0.01
        slender = HullParticulars(; length_pp = 320.0, beam = 40.0, draft = 20.8,
            block_coefficient = 0.6)
        @test hirano_takashina_derivatives(slender).Y_v < 0
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError HullParticulars(; length_pp = -1.0, beam = 1.0,
            draft = 1.0, block_coefficient = 0.8)
        @test_throws ArgumentError empirical_sway_derivatives(kvlcc2;
            method = :not_a_method)
    end
end
