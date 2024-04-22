
using DelimitedFiles
using Statistics
using Plots
using StatsPlots

include("config.jl")

function plot_figure1(config)
    nbins = config.Nd
    files = ["results/lower_bound_nd$(n).txt" for n in nbins]
    labels = ["#nodes=$(n)" for n in nbins]
    lbs = readdlm.(files)

    fig = plot()
    for (lb, lab) in zip(lbs, labels)
        plot!(lb, label=lab, lw=2.0)
    end
    xlabel!("SDDP iteration")
    ylabel!("SDDP lower-bound")
    title!("Convergence of MC-SDDP")
    savefig("convergence_sddp.pdf")
    return fig
end

function plot_figure2(config)
    nbins = config.Nd
    insamp_files = ["results/simulation_$(n)_insamp.txt" for n in nbins]
    outsamp_files = ["results/simulation_$(n)_outsamp.txt" for n in nbins]
    lb_files = ["results/lower_bound_nd$(n).txt" for n in nbins]

    colors = palette(:darktest, 3)
    fig = plot(layout=(length(nbins), 1), link=:both, legend=nothing, size=(500, 500), bottom_margin=0*Plots.mm, yaxis=nothing,
               guidefontsize=7,
               legendfontsize=6,
    )
    k = 0
    for (f1, f2, f3, nb) in zip(insamp_files, outsamp_files, lb_files, nbins)
        k += 1
        insamp_gain = readdlm(f1)
        outsamp_gain = readdlm(f2)
        lb = readdlm(f3)

        leg = (k == 1)
        if k != length(nbins)
            plot!(xaxis=nothing, subplot=k)
        end
        density!(insamp_gain, label="in-sample model", subplot=k, lw=1.0, color=colors[1], legend=leg, ls=:dot)
        density!(outsamp_gain, label="out-sample model", subplot=k, lw=1.0, color=colors[2], ls=:dot)

        vline!([lb[end]], subplot=k, label="SDDP LB", lw=3.0, color=colors[3])
        vline!([mean(insamp_gain)], subplot=k, label="avg in-sample", lw=3.0, color=colors[1])
        vline!([mean(outsamp_gain)], subplot=k, label="avg out-sample", lw=3.0, color=colors[2])

        # annotate!(4.92, 60.00006, text("#nodes=$(nb)", 8), subplot=k)
        ylabel!("#nodes=$(nb)", subplot=k)
    end
    # title!("Density", subplot=1)
    xlabel!("Utility", subplot=6)
    savefig("density_sddp.pdf")
    return fig
end

function plot_figure3(config)
    nbins = config.Nd
    results = readdlm("results/sensitivity_std.txt")

    res = results[:, 2] .== 2.0
    display(results[res, :])

    colors = palette(:berlin10, 6)
    fig = plot(layout=(2, 1), link=:all)
    for (k, nb) in enumerate(nbins)
        scan = results[:, 2] .== nb
        σ = results[scan, 1]
        plot!(σ, results[scan, 3], subplot=1, color=colors[k], marker=:d, label="#nodes=$(nb)")
        plot!(xaxis=nothing, subplot=1)
    end
    for (k, nb) in enumerate(nbins)
        scan = results[:, 2] .== nb
        σ = results[scan, 1]
        plot!(σ, results[scan, 4], subplot=2, color=colors[k], marker=:d, label="#nodes=$(nb)", legend=false)
    end
    xlabel!("Standard deviation σ", subplot=2)
    ylabel!("Out-of-sample", subplot=2)
    ylabel!("SDDP LB", subplot=1)
    title!("Expected utility", subplot=1)
    savefig("sensitivity_std.pdf")
    return fig
end

f1 = plot_figure1(config)
f2 = plot_figure2(config)
f3 = plot_figure3(config)
