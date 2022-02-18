module parameter_structs

Base.@kwdef mutable struct SimulationTimeParameters
    time_min::Float64
    time_max::Float64
    time_delta::Float64

    tspan = (time_min, time_max) # for ODE solutions
    ts = time_min:time_delta:time_max
end

Base.@kwdef struct VirusLumenParameters
    gamma::Float64
    kmax::Float64
    a::Float64
end

Base.@kwdef struct VirusMidgutParameters
    zeta::Float64
    virus_replication_rate::Float64
    diffusion_rate::Float64
    gamma::Float64 # rate at which virus passes out through BL
end

Base.@kwdef struct VirusHemolymphParameters
    virus_replication_rate::Float64
    diffusion_rate::Float64
end

Base.@kwdef struct VirusSalivaryGlandsParameters
    virus_replication_rate::Float64
end

Base.@kwdef struct MechanicalParameters
    eta::Float64
    v_star::Float64
end

Base.@kwdef struct InitialConditions
    virus_lumen::Float64
    virus_midgut::Array{Float64}
    virus_hemolymph::Array{Float64}
    virus_salivary_glands::Float64
    volume::Float64
end

Base.@kwdef mutable struct SpatialDomainParameters
    delta_x::Float64
    x_max::Float64
    delta_y::Float64
    y_max::Float64

    x_range = 0:delta_x:x_max
    y_range = 0:delta_y:y_max

    # determine y index closest to y=1
    y_boundary_index = argmin(map(y -> abs(y - 1.0), y_range))
end

Base.@kwdef struct RefeedingParameters
    time_refeed::Float64
    volume_after_refeed::Float64
    virus_added::Float64
end

end
