using Test


@testset "barrier_index" begin

    delta_x = 0.1

    x_b = 0.05
    index = MosquitoArbovirus.time_dependent_porosity.barrier_index(x_b, delta_x)
    @test index == 1

    x_b = 0.1
    index = MosquitoArbovirus.time_dependent_porosity.barrier_index(x_b, delta_x)
    @test index == 2

    x_b = 0.10001
    index = MosquitoArbovirus.time_dependent_porosity.barrier_index(x_b, delta_x)
    @test index == 2
end;
