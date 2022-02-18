module ode

using DifferentialEquations
export solve_ode_generic, midgut_initial_rate, simulate_virus_in_lumen, simulate_lumen_volume


function solve_ode_generic(tspan, inits, derivative_function, parameter_list)
    prob = ODEProblem(derivative_function, inits, tspan, parameter_list)
    sol = solve(prob)
    sol
end

"""
Derivative of ODE specifying the lumen virus quantity.
"""
function dldt(l, p, t)
    lumen_params = p[1]
    gamma = lumen_params.gamma
    kmax = lumen_params.kmax
    a = lumen_params.a

    - gamma * l - kmax * l^2 / (a^2 + l^2)
end

"""
Computes the rate at which virus enters the midgut (i.e. at ``x=0``). This function
returns a function of time which specifies this quantity.
"""
function midgut_initial_rate(virus_lumen_solution, lumen_parameters)
    kmax = lumen_parameters.kmax
    a = lumen_parameters.a
    l = virus_lumen_solution
    t -> kmax * l(t)^2 / (a^2 + l(t)^2)
end

function simulate_virus_in_lumen(simulation_parameters, initial_conditions,
    virus_lumen_parameters, refeeding_parameters)
    time_refeed = refeeding_parameters.time_refeed
    tspan = simulation_parameters.tspan
    time_end = tspan[2]
    time_diff = time_end - time_refeed
    l0 = initial_conditions.virus_lumen
    if time_diff <= 0
        l_sol = solve_ode_generic(tspan, l0, dldt, [virus_lumen_parameters])
    else
        # solve from 1st feed - 2nd
        tspan_1 = (tspan[1], time_refeed)
        l_sol_1 = solve_ode_generic(tspan_1, l0, dldt, [virus_lumen_parameters])

        # solve from 2nd feed to end
        tspan_2 = (tspan[1], time_diff)
        l0_after = l_sol_1(time_refeed) + refeeding_parameters.virus_added
        l_sol_2 = solve_ode_generic(tspan_2, l0_after, dldt, [virus_lumen_parameters])

        function l_sol(t)
            if t < time_refeed
                l_sol_1(t)
            else
                l_sol_2(t - time_refeed)
            end
        end
    end
    l_sol
end

"""
Calculates lumen volume over time via an analytical solution of:

``\\frac{dv}{dt} = -\\eta (v - v^*)``

where refeeding at time ``t_1`` causes volume to discontinuously change to
a new volume, after which the above ODE dictates the dynamics.
"""
function simulate_lumen_volume(initial_conditions, mechanical_parameters,
                               refeeding_parameters)
    v0 = initial_conditions.volume
    time_refeed = refeeding_parameters.time_refeed
    v_after = refeeding_parameters.volume_after_refeed
    eta = mechanical_parameters.eta
    v_star = mechanical_parameters.v_star

    function vol_sol(t)
        if t < time_refeed
            v_star + exp(-eta * t) * (v0 - v_star)
        else
            v_star + exp(-eta * (t - time_refeed)) * (v_after - v_star)
        end
    end
    vol_sol
end

end
