# Ben: this is a convenience file that makes it quicker to run tests

include("../src/MosquitoArbovirus.jl")
using .MosquitoArbovirus

include("runtests_finite_differences.jl")
include("runtests_porosity.jl")
