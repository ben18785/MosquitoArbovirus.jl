module time_dependent_porosity

using DifferentialEquations
include("ode.jl")
using .ode
include("finite_differences.jl")
using .finite_differences
using Interpolations
using QuadGK
using DataFrames
using LinearAlgebra

"""
Computes a time-varying function specifying the point in the x plane where the
BL becomes porous to the virus.
"""
function porosity_position(
    volume_solution, mechanical_parameters,
    midgut_parameters)

    first = 3.0 / 4.0 * 1.0 / pi
    power = 1.0 / 3.0
    v_star = mechanical_parameters.v_star
    zeta = midgut_parameters.zeta
    r = t -> first * volume_solution(t)^power
    r_star = first * v_star^power
    function x_b(t)
        pos = 1.0 - zeta * (r(t) - r_star)^2
        if pos < 0
            pos = 0
        end
        pos
    end
    x_b
end

"""
Returns the last index of the discretised x domain which is not in the porous
layer.
"""
function barrier_index(x_b, delta_x)
    # need +1 since starts at 0; need extra little bit to avoid underflow
    Int(floor((x_b + 1e-8) / delta_x) + 1)
end

"""
Computes time-derivative of x-discretised PDE describing midgut virus
dynamics. The PDE has the form:

``\\partial m(x, t)/dt = D * \\partial^2 m(x, t)/\\partial x^2 + growth_rate * m(x, t) * (1 - m(x, t)) - 1(x > x_b)\\gamma x``

where ``1(x > x_b)`` is an indicator function equal to 1 if true and ``0\\leq x_b \\ 1`` is a dynamic
porosity barrier position.

"""
function dmdt(m, p, t)
    laplacian_matrix = p[1]
    vm_parameters = p[2]
    midgut_rate_func = p[3]
    spatial_parameters = p[4]
    barrier_func = p[5]
    viral_growth_rate = vm_parameters.virus_replication_rate
    D = vm_parameters.diffusion_rate
    gamma = vm_parameters.gamma

    delta_x = spatial_parameters.delta_x
    x_b = barrier_func(t)
    last_non_porous_index = barrier_index(x_b, delta_x)

    dmdt = D * laplacian_matrix * m + viral_growth_rate * m .* (1 .- m)
    dmdt[1] += midgut_rate_func(t) / delta_x
    dmdt[(last_non_porous_index + 1):end] -= gamma * m[(last_non_porous_index + 1):end]
    dmdt
end

function simulate_virus_in_midgut(tspan, initial_conditions,
                                  virus_lumen_parameters,
                                  virus_midgut_parameters,
                                  mechanical_parameters,
                                  spatial_domain_parameters,
                                  virus_lumen_solution,
                                  volume_solution,
                                  barrier_func)

   # determine viral flux at x=0
   midgut_rate_func = midgut_initial_rate(virus_lumen_solution,
                                          virus_lumen_parameters)

    # solve discretised spatial PDE
    m0 = initial_conditions.virus_midgut
    delta_x = spatial_domain_parameters.delta_x
    x_max = spatial_domain_parameters.x_max
    laplacian_matrix = second_derivative(delta_x, x_max)
    params = [laplacian_matrix, virus_midgut_parameters,
              midgut_rate_func, spatial_domain_parameters,
              barrier_func]
    m_sol = solve_ode_generic(tspan, m0, dmdt, params)
    m_sol
end

function total_virus_in_tissue(t, virus_tissue_solution, xrange)
    m = virus_tissue_solution(t)
    f = LinearInterpolation(xrange, m)
    quadgk(f, 0, 1)[1]
end

function total_virus_leaving_midgut(t, virus_midgut_solution, xrange,
                                    barrier_func)
    m = virus_midgut_solution(t)
    f = LinearInterpolation(xrange, m)
    x_lower = barrier_func(t)
    quadgk(f, x_lower, 1)[1]
end

function dhdt(h, p, t)
    laplacian_matrix = p[1]
    vh_parameters = p[2]
    bl_rate_func = p[3]
    spatial_parameters = p[4]

    viral_growth_rate = vh_parameters.virus_replication_rate
    D = vh_parameters.diffusion_rate
    delta_y = spatial_parameters.delta_y

    dhdt = D * laplacian_matrix * h + viral_growth_rate * h .* (1 .- h)
    dhdt[1] += bl_rate_func(t) / delta_y
    dhdt
end

# implements theta rule from: http://hplgit.github.io/num-methods-for-PDEs/doc/pub/diffu/html/._diffu001.html
function solve_theta_rule(tspan, initial_conditions,
                          virus_hemolymph_parameters,
                          virus_midgut_solution,
                          spatial_domain_parameters,
                          barrier_func)



end

function simulate_virus_in_hemolymph_mol(tspan, initial_conditions,
                                  virus_hemolymph_parameters,
                                  virus_midgut_solution,
                                  spatial_domain_parameters,
                                  barrier_func)

    # rate at which virus leaves midgut, which determines flux at x=0
    x_range = spatial_domain_parameters.x_range
    function bl_rate_func(t)
        total_virus_leaving_midgut(t, virus_midgut_solution, x_range,
                                   barrier_func)
    end

    # solve discretised spatial PDE
    h0 = initial_conditions.virus_hemolymph
    delta_y = spatial_domain_parameters.delta_y
    y_max = spatial_domain_parameters.y_max
    laplacian_matrix = second_derivative(delta_y, y_max)
    params = [laplacian_matrix, virus_hemolymph_parameters,
              bl_rate_func, spatial_domain_parameters]
    h_sol = solve_ode_generic(tspan, h0, dhdt, params)
    h_sol
end

function simulate_virus_in_hemolymph_implicit(tspan, initial_conditions,
                                  virus_hemolymph_parameters,
                                  virus_midgut_solution,
                                  spatial_domain_parameters,
                                  barrier_func)
  # rate at which virus leaves midgut, which determines flux at x=0
  x_range = spatial_domain_parameters.x_range
  function bl_rate_func(t)
      total_virus_leaving_midgut(t, virus_midgut_solution, x_range,
                                 barrier_func)
  end

  # solve discretised spatial PDE implicitly
  h0 = initial_conditions.virus_hemolymph
  delta_y = spatial_domain_parameters.delta_y
  y_max = spatial_domain_parameters.y_max
  delta_t = spatial_domain_parameters.time_delta
  theta = spatial_domain_parameters.theta
  viral_growth_rate = virus_hemolymph_parameters.virus_replication_rate
  D = virus_hemolymph_parameters.diffusion_rate

  A = implicit_A(delta_y, y_max, delta_t, D, theta)
  C = implicit_C(delta_y, y_max, delta_t, D, theta)

  # solve system
  N = length(h0)
  ts = delta_t:delta_t:tspan[2]
  nt = length(ts)
  H = zeros(nt + 1, N)
  h = h0
  H[1, :] = h
  for i in 1:nt
      flux = bl_rate_func(ts[i])
      d = implicit_d(h, delta_y, y_max, delta_t, D, theta, viral_growth_rate, flux)
      b = C * h + d
      h = A \ b
      H[i + 1, :] = h
  end

  # create an array of interpolation functions: one for each x
  H_dims = size(H)
  nx = H_dims[2]
  ts1 = Vector(0.0:delta_t:tspan[2])
  fvec = [LinearInterpolation(ts1, H[:, i]) for i in 1:nx]
  function h_sol(t)
     [fvec[i](t) for i in 1:nx]
  end
  h_sol
end

function simulate_virus_in_hemolymph(tspan, initial_conditions,
                                  virus_hemolymph_parameters,
                                  virus_midgut_solution,
                                  spatial_domain_parameters,
                                  barrier_func)

    solver_type = spatial_domain_parameters.solver_type
    if solver_type == "mol"
        h_sol = simulate_virus_in_hemolymph_mol(tspan, initial_conditions,
                                          virus_hemolymph_parameters,
                                          virus_midgut_solution,
                                          spatial_domain_parameters,
                                          barrier_func)
    elseif solver_type == "implicit"
        h_sol = simulate_virus_in_hemolymph_implicit(tspan, initial_conditions,
                                          virus_hemolymph_parameters,
                                          virus_midgut_solution,
                                          spatial_domain_parameters,
                                          barrier_func)
    end
    h_sol
end

function dhdy_flux(virus_hemolymph_solution, virus_hemolymph_parameters,
    spatial_domain_parameters)
    delta_y = spatial_domain_parameters.delta_y
    y_boundary_index = spatial_domain_parameters.y_boundary_index
    function dhdy_flux_func(t)
        hs = virus_hemolymph_solution(t)
        h_n = hs[y_boundary_index]
        h_n_minus_1 = hs[y_boundary_index - 1]
        D = virus_hemolymph_parameters.diffusion_rate
        - D * (h_n - h_n_minus_1) / delta_y
    end
    dhdy_flux_func
end

function dsdt(s, p, t)
    vs_parameters = p[1]
    dhdy_flux_func = p[2]
    viral_growth_rate = vs_parameters.virus_replication_rate
    viral_growth_rate * s * (1 - s) + dhdy_flux_func(t)
end

function simulate_virus_in_salivary_glands(tspan, initial_conditions,
    virus_salivary_glands_parameters,
    virus_hemolymph_solution,
    virus_hemolymph_parameters,
    spatial_domain_parameters)

    # flux into salivary glands = -D dhdy|y=1
    dhdy_flux_func = dhdy_flux(virus_hemolymph_solution, virus_hemolymph_parameters,
        spatial_domain_parameters)

    # solve ode
    s0 = initial_conditions.virus_salivary_glands
    params = [virus_salivary_glands_parameters, dhdy_flux_func]
    s_sol = solve_ode_generic(tspan, s0, dsdt, params)
    s_sol
end

function simulate_porosity(simulation_time_parameters,
    virus_lumen_parameters,
    virus_midgut_parameters,
    virus_hemolymph_parameters,
    virus_salivary_glands_parameters,
    mechanical_parameters,
    initial_conditions,
    spatial_domain_parameters,
    refeeding_parameters
    )
    tspan = simulation_time_parameters.tspan
    l_solution = simulate_virus_in_lumen(tspan, initial_conditions,
                                         virus_lumen_parameters,
                                         refeeding_parameters)

    volume_solution = simulate_lumen_volume(initial_conditions,
                                            mechanical_parameters,
                                            refeeding_parameters)

    # determine porosity barrier position
    barrier_func = porosity_position(
        volume_solution, mechanical_parameters,
        virus_midgut_parameters)

    m_solution = simulate_virus_in_midgut(tspan, initial_conditions,
        virus_lumen_parameters, virus_midgut_parameters,
        mechanical_parameters, spatial_domain_parameters,
        l_solution, volume_solution, barrier_func)

    h_solution = simulate_virus_in_hemolymph(tspan, initial_conditions,
        virus_hemolymph_parameters, m_solution,
        spatial_domain_parameters, barrier_func)

    s_solution = simulate_virus_in_salivary_glands(tspan,
            initial_conditions,
            virus_salivary_glands_parameters,
            h_solution,
            virus_hemolymph_parameters,
            spatial_domain_parameters)

    # calculate totals in midgut and hemolymph
    ts = simulation_time_parameters.ts

    x_range = spatial_domain_parameters.x_range
    total_midgut = map(t -> total_virus_in_tissue(t, m_solution, x_range), ts)

    y_range = spatial_domain_parameters.y_range
    total_hemolymph = map(t -> total_virus_in_tissue(t, h_solution, y_range), ts)

    (virus_lumen=l_solution, virus_midgut=m_solution,
     virus_hemolymph=h_solution, virus_salivary_glands=s_solution,
     total_midgut=total_midgut, total_hemolymph=total_hemolymph,
     lumen_volume=volume_solution,
     barrier_func=barrier_func)
end

function dataframe_aggregates(simulation_time_parameters, system_solution)
    ts = simulation_time_parameters.ts
    sims = system_solution
    lumen = map(t -> sims.virus_lumen(t), ts)
    midgut = sims.total_midgut
    hemolymph = sims.total_hemolymph
    salivary_glands = map(t -> sims.virus_salivary_glands(t), ts)
    volume = map(sims.lumen_volume, ts)
    barrier = map(sims.barrier_func, ts)
    DataFrame(time=ts, lumen=lumen, midgut=midgut,
              hemolymph=hemolymph, salivary_glands=salivary_glands,
              volume=volume, barrier=barrier)
end

function dataframe_spatial(simulation_time_parameters, tissue_solution,
    x_range)
    ts = simulation_time_parameters.ts
    for i in 1:length(ts)
        temp = tissue_solution(ts[i])
        df_tmp = DataFrame(time=ts[i], value=temp, x=x_range)
        if i == 1
            global df_big = df_tmp
        else
            global df_big = [df_big; df_tmp]
        end
    end
    df_big
end

end
