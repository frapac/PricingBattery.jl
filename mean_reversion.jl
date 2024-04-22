
using DelimitedFiles
using Statistics

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

epex_price = readdlm("epexprice.txt")
results = detrend(epex_price)
writedlm("meanreversion.txt", results.mean)

