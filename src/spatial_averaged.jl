module spatial_averaged

using DifferentialEquations
include("ode.jl")
using .ode


"""
Derivative of ODE specifying the lumen virus quantity.
"""
function dldt(l, p, t)
    lumen_params = p[1]
    gamma = lumen_params.gamma
    kmax = lumen_params.kmax
    a = lumen_params.a
    X = p[2]

    - gamma * l - kmax * X(t) * l^2 / (a^2 + l^2)
end


# contains Helen Byrne's spatially averaged equations
function dmdt(m, p, t)
    lumen_params = p[1]
    kmax = lumen_params.kmax
    a = lumen_params.a
    X = p[2]
    vm_parameters = p[3]
    alpha = vm_parameters.virus_replication_rate
    kappa = vm_pararameters.carrying_capacity
    keh = virus_midgut_parameters.gamma
    mech_parameters = p[4]
    Xstar = mech_parameters.v_star

    dmdt = (
        kmax * l(t)^2 * X(t) / (a^2 + l(t)^2) + alpha * m * (1 - m / kappa) -
        keh * m * (X(t) - Xstar)
    )
    dmdt
end

function dhdt(h, p, t)
    vm_parameters = p[1]
    vh_parameters = p[2]
    alpha = vh_parameters.virus_replication_rate
    kappa = vh_pararameters.carrying_capacity
    keh = virus_midgut_parameters.gamma
    
    keh * m * (X(t) - Xstar) + alpha * h (1 - h / kappa)
end


function dsdt(s, p, t)
    vs_parameters = p[1]
    alpha = vs_parameters.virus_replication_rate
    kappa = vs_pararameters.carrying_capacity
    alpha * s * (1 - s / kappa)
end

end
