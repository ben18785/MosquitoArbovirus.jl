using Test
using LinearAlgebra
using Statistics

@testset "neumann" begin
    xmax = 0.5
    delta_x = 0.1
    A = MosquitoArbovirus.finite_differences.second_derivative(delta_x, xmax)
    nrows = Int(xmax / delta_x + 1)
    sz = size(A)
    @test nrows == sz[1] == sz[2]

    cons = 1 / delta_x^2

    # diagonal
    diag = LinearAlgebra.diag(A, 0)
    diag_mean = Statistics.mean(diag)
    @test diag_mean ≈ -2 * cons

    # upper diagonal
    offdiag = LinearAlgebra.diag(A, 1)
    @test offdiag[1] == 2 * cons
    offdiag_mean = Statistics.mean(offdiag[2:(nrows - 1)])
    @test offdiag_mean ≈ cons

    # lower diagonal
    offdiag = LinearAlgebra.diag(A, -1)
    @test offdiag[(nrows - 1)] == 2 * cons
    offdiag_mean = Statistics.mean(offdiag[1:(nrows - 2)])
    @test offdiag_mean ≈ cons

end;

@testset "dirichlet" begin
    xmax = 1.0
    delta_x = 0.01
    A = MosquitoArbovirus.finite_differences.second_derivative(delta_x, xmax, "dirichlet")
    nrows = Int(xmax / delta_x + 1)
    sz = size(A)
    @test nrows == sz[1] == sz[2]

    cons = 1 / delta_x^2

    # diagonal
    diag = LinearAlgebra.diag(A, 0)
    diag_mean = Statistics.mean(diag[1:(nrows - 1)])
    @test diag_mean ≈ -2 * cons

    # upper diagonal
    offdiag = LinearAlgebra.diag(A, 1)
    @test offdiag[1] == 2 * cons
    @test offdiag[nrows - 1] == 0.00
    offdiag_mean = Statistics.mean(offdiag[2:(nrows - 2)])
    @test offdiag_mean ≈ cons

    # lower diagonal
    offdiag = LinearAlgebra.diag(A, -1)
    @test offdiag[nrows - 1] == 0.0
    offdiag_mean = Statistics.mean(offdiag[1:(nrows - 2)])
    @test offdiag_mean ≈ cons
end;
