# MosquitoArbovirus.jl
A Julia package for simulating viral infection dynamics in mosquitoes.

Can be installed via the Julia REPL by first pressing `]` then:

`add https://github.com/ben18785/MosquitoArbovirus.jl`


# Example: single infectious bloodmeal
The following code simulates viral dynamics following a single infectious bloodmeal on a mosquito at `t=0`. See the paper for parameter definitions.

```
using Pkg
using MosquitoArbovirus
using Plots
using CSV

ma = MosquitoArbovirus.time_dependent_porosity
st = MosquitoArbovirus.parameter_structs

# define parameters and initial conditions
x_max = 1.0
delta_x = 0.01
y_max = 2.0
delta_y = 0.02
simulation_parameters = st.SimulationParameters(
    time_min=0.0, time_max=10.0, time_delta_pde=0.05,
    delta_x=delta_x, x_max=x_max,
    delta_y=delta_y, y_max=y_max,
    time_delta=0.05, theta=0.5,
    solver_type="implicit")
refeeding_parameters = st.RefeedingParameters(
    time_refeed=15.0,
    volume_after_refeed=1.0,
    virus_added=0.0)
virus_lumen_parameters = st.VirusLumenParameters(
    gamma=1.5, kmax=0.0001, a=1.0)
virus_growth_rate = 2.7
virus_midgut_parameters = st.VirusMidgutParameters(
    zeta=50.0, virus_replication_rate=virus_growth_rate,
    diffusion_rate=0.0055,
    gamma=1.0)
virus_hemolymph_growth_rate = 3 * virus_growth_rate
virus_hemolymph_parameters = st.VirusHemolymphParameters(
    virus_replication_rate=virus_hemolymph_growth_rate,
    diffusion_rate=0.033)
virus_salivary_glands_parameters = st.VirusSalivaryGlandsParameters(
    virus_replication_rate=1.6 * virus_growth_rate)
mechanical_parameters = st.MechanicalParameters(
    eta=2.0, v_star=0.1)
m0 = zeros(Int(x_max / delta_x) + 1)
h0 = zeros(Int(y_max / delta_y) + 1)
initial_conditions = st.InitialConditions(
    virus_lumen=1.0, virus_midgut=m0,
    virus_hemolymph=h0, virus_salivary_glands=0.0,
    volume=1.0)

# simulate
sims = ma.simulate(simulation_parameters,
    virus_lumen_parameters,
    virus_midgut_parameters,
    virus_hemolymph_parameters,
    virus_salivary_glands_parameters,
    mechanical_parameters,
    initial_conditions,
    refeeding_parameters
    )

# plot aggregate quantities on a single graph
ts = simulation_parameters.ts
lumen = map(t -> sims.virus_lumen(t), ts)
salivary_gland = map(t -> sims.virus_salivary_glands(t), ts)
midgut = sims.total_midgut
hemolymph = sims.total_hemolymph
plot(ts, lumen, label="lumen")
plot!(ts, midgut, label="midgut")
plot!(ts, hemolymph, label="hemolymph")
plot!(ts, salivary_gland, label="salivary gland")
```
