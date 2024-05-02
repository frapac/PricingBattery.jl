
using DelimitedFiles
using Statistics
using Plots

DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

function detrend(input)
    @assert length(input) == 8760
    ndays = 365
    nhours = ndays * 24
    nweeks = 52 + 1

    detrended = copy(input)
    # Remove annual variation
    Δ_w = 7 * 24
    weekly_price = zeros(nweeks)
    for i in 1:nweeks
        k = (i-1)*Δ_w + 1
        l = min(i * Δ_w, nhours)
        weekly_price[i] = mean(input[k:l])
    end

    for i in 1:nhours
        k = div(i, Δ_w) + 1
        detrended[i] -= weekly_price[k]
    end

    # Remove weekly variation
    days_per_week = 7
    daily_price = zeros(ndays)
    for i in 1:ndays
        k = (i-1) * 24 + 1
        l = i * 24
        daily_price[i] = mean(detrended[k:l])
    end
    avg_daily_price = zeros(days_per_week)
    for i in 1:days_per_week
        avg_daily_price[i] = mean(daily_price[i:days_per_week:end])
    end


    Δ_d = 24
    for i in 1:nhours-1
        k = (div(i, Δ_d) % days_per_week) + 1
        detrended[i] -= avg_daily_price[k]
    end

    # Remove daily variation
    hours_per_day = 24
    hourly_price = zeros(hours_per_day)
    for i in 1:hours_per_day
        hourly_price[i] = mean(detrended[i:Δ_d:end])
    end

    for i in 1:nhours
        k = i % hours_per_day + 1
        detrended[i] -= hourly_price[k]
    end

    mean_reversion = zeros(nhours)
    for i in 1:nhours
        i_w = div(i, Δ_w) + 1
        i_d = (div(i, Δ_d) % days_per_week) + 1
        i_h = i % hours_per_day + 1
        mean_reversion[i] = weekly_price[i_w] + avg_daily_price[i_d] + hourly_price[i_h]
    end
    return (mean=mean_reversion, week=weekly_price, day=avg_daily_price, hour=hourly_price, residual=detrended)
end

function plot_mean_reversion(input, results)
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    pos = 15*24:30*24:365*24

    fig = plot(layout=(2, 1), link=:both, xticks=(pos, months), ylims=(2, 5), xlims=(1, 365*24))
    plot!(input, lw=0.5, label="", subplot=1)
    plot!(results.mean, label="", lw=0.5, subplot=2)

    title!("Original", subplot=1)
    title!("Mean reversion", subplot=2)
    for s in [1,2]
        vline!(24 .* cumsum([0; DAYS_PER_MONTH]), color="black", label="", subplot=s)
    end
    ylabel!("Price [€/MWh]")
    savefig("time-series.pdf")
    return fig
end

function plot_mean_reversion_decomposition(results)
    fig = plot(layout=(3, 1), link=:y)
    ylabel!("Price")
    bar!(results.week, subplot=1, label="", xlims=(0.5, 52.5))
    xlabel!("Week", subplot=1)
    title!("Annual variation", subplot=1)
    plot!(subplot=2, xticks=(1:7, ["mon", "tue", "wed", "thu", "fri", "sat", "sun"]))
    bar!(results.day[[5, 6, 7, 1, 2, 3, 4]], subplot=2, label="", xlims=(0.5, 7.5))
    xlabel!("Day", subplot=2)
    title!("Weekly variation", subplot=2)
    bar!(results.hour, subplot=3, label="", xlims=(0.5, 24.5))
    xlabel!("Hour", subplot=3)
    title!("Daily variation", subplot=3)
    savefig("mean-reversion-decomposition.pdf")
    return fig
end

epex_price = readdlm("epexprice.txt")
results = detrend(log.(epex_price))
f1 = plot_stats(epex_price, results)
f2 = plot_mean_reversion_decomposition(results)
writedlm("logmeanreversion.txt", results.mean)

