using Comonicon

include("config.jl")
include("markov_sddp.jl")

RESULTS_DIR = joinpath(@__DIR__, "results")

if !isdir(RESULTS_DIR)
    mkpath(RESULTS_DIR)
end

@main function main(
    arg;
    max_iter=300,
    risk_aversion=1e-4,
    Nd=8,
)
    if arg == "fig2"
        config = Config(
            ARParams(0.0, 0.54, 0.39),
            ModelParams(364, 0.2, 1e5, 0.0, 5.0),
            SDDPParams(500, 1000),
            risk_aversion,
        )
        compute_figure2(config, [2^i for i in 0:6])
    elseif arg == "fig3"
        config = Config(
            ARParams(0.0, 0.54, 0.39),
            ModelParams(364, 0.2, 1e5, 0.0, 5.0),
            SDDPParams(max_iter, 1000),
            risk_aversion,
        )
        sigmas = 0.4:0.05:0.8
        compute_figure3(config, Nd, sigmas)
    elseif arg == "fig4"
        config = Config(
            ARParams(0.0, 0.54, 0.39),
            ModelParams(364, 0.2, 1e5, 0.0, 5.0),
            SDDPParams(max_iter, 1000),
            0.0001,
        )
        rhos = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
        compute_figure4(config, Nd, rhos)
    elseif arg == "fig5"
        config = Config(
            ARParams(0.0, 0.54, 0.39),
            ModelParams(364, 0.2, 1e5, 0.0, 5.0),
            SDDPParams(max_iter, 1000),
            risk_aversion,
        )
        alphas = 0.2:0.2:1.0
        compute_figure4(config, Nd, alphas)
    end
end
