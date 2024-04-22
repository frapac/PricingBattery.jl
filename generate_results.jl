include("config.jl")
include("markov_sddp.jl")

RESULTS_DIR = joinpath(@__DIR__, "results")

if !isdir(RESULTS_DIR)
    mkpath(RESULTS_DIR)
end

compute_figure1(config)
compute_figure2(config)
compute_figure3(config)

