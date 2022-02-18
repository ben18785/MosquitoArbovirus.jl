module time_dependent_diffusion

using DifferentialEquations
include("ode.jl")
using .ode
include("finite_differences.jl")
using .finite_differences

"""
Computes a time-varying function specifying the diffusion rate of
virus in the midgut layer, where the virus moves towards
penetrating the basal lamina layer.
"""
function diffusion_rate_function(
    volume_solution, mechanical_parameters,
    midgut_parameters)

    first = 3.0 / 4.0 * 1.0 / pi
    power = 1.0 / 3.0
    v_star = mechanical_parameters.v_star
    zeta = midgut_parameters.zeta
    r = t -> first * volume_solution(t)^power
    r_star = first * v_star^power
    t -> zeta * (r(t) - r_star)^2
end

"""
Computes time-derivative of x-discretised PDE describing midgut virus
dynamics. The derivatives take the form:

``dm(x_i, t)/dt = D * central_diff_laplacian(i) + growth_rate * m(x_i, t) * (1 - m(x_i, t))``

"""
function dmdt(m, p, t)
    D_func = p[1]
    laplacian_matrix = p[2]
    vm_parameters = p[3]
    midgut_rate_func = p[4]
    spatial_parameters = p[5]
    viral_growth_rate = vm_parameters.virus_replication_rate
    delta_x = spatial_parameters.delta_x

    D = D_func(t)
    dmdt = D * laplacian_matrix * m + viral_growth_rate * m .* (1 .- m)
    dmdt[1] += midgut_rate_func(t) / (delta_x)
    dmdt
end


function simulate_virus_in_midgut(tspan, initial_conditions,
                                  virus_lumen_parameters,
                                  virus_midgut_parameters,
                                  mechanical_parameters,
                                  spatial_domain_parameters,
                                  virus_lumen_solution,
                                  volume_solution)

   # determine viral flux at x=0
   midgut_rate_func = midgut_initial_rate(virus_lumen_solution, virus_lumen_parameters)

    # determine diffusion rate of virus through BL
    D_func = diffusion_rate_function(
        volume_solution,
        mechanical_parameters,
        virus_midgut_parameters)

    # solve discretised spatial PDE
    m0 = initial_conditions.virus_midgut
    delta_x = spatial_domain_parameters.delta_x
    x_max = spatial_domain_parameters.x_max
    laplacian_matrix = second_derivative(delta_x, x_max)
    params = [D_func, laplacian_matrix, virus_midgut_parameters,
              midgut_rate_func, spatial_domain_parameters]
    m_sol = solve_ode_generic(tspan, m0, dmdt, params)
    m_sol
end

function simulate_diffusion_rate(tspan,
    virus_lumen_parameters,
    virus_midgut_parameters,
    mechanical_parameters,
    initial_conditions,
    spatial_domain_parameters,
    refeeding_parameters
    )

    l_solution = simulate_virus_in_lumen(tspan, initial_conditions,
                                         virus_lumen_parameters,
                                         refeeding_parameters)

    volume_solution = simulate_lumen_volume(tspan, initial_conditions,
                                            mechanical_parameters,
                                            refeeding_parameters)

    m_solution = simulate_virus_in_midgut(tspan, initial_conditions,
        virus_lumen_parameters, virus_midgut_parameters,
        mechanical_parameters, spatial_domain_parameters,
        l_solution, volume_solution)

    (virus_lumen=l_solution, virus_midgut=m_solution, lumen_volume=volume_solution)
end


end # module
