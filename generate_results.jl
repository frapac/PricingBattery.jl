include("markov_sddp.jl")

RESULTS_DIR = joinpath(@__DIR__, "results")

if !isdir(RESULTS_DIR)
    mkpath(RESULTS_DIR)
end

maxit = 50

compute_figure1(; total_it=maxit)
compute_figure2(; total_it=maxit)
compute_figure3(; total_it=maxit)

