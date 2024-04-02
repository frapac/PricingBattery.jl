#=
    Implement the Markov chain discretization using the Tauchen-Hussey algorithm.
=#

using Random
using Distributions
using FastGaussQuadrature

abstract type StochasticProcess end

# z_{t+1} = (1 - phi) * mu + phi z_t + sigma εₜ
struct AR1{T} <: StochasticProcess
    mu::T
    phi::T
    sigma::T
end

function simulate(process::AR1, x0, T, nscen)
    noise = Normal(0.0, process.sigma)
    scenarios = zeros(T+1, nscen)
    for n in 1:nscen
        scenarios[1, n] = x0
        for t in 1:T
            epsilon = rand(noise)
            scenarios[t+1, n] = (1.0 - process.phi) * process.mu + process.phi * scenarios[t, n] + epsilon
        end
    end
    return scenarios
end

struct StationaryMarkovChain{T} <: StochasticProcess
    x::Vector{T}
    proba::Matrix{T}
end

function _sample(probs)
    n = length(probs)
    u = rand()
    cumsum = 0.0
    j = -1
    for k in 1:n
        cumsum += probs[k]
        if cumsum >= u
            j = k
            break
        end
    end
    println(probs)
    return j
end

function project(process::StationaryMarkovChain, x0)
    nx = length(process.x)
    ind0, val_min = -1, Inf
    for i in 1:nx
        current_val = abs(process.x[i] - x0)
        if current_val < val_min
            ind0 = i
            val_min = current_val
        end
    end
    return ind0, process.x[ind0]
end

function simulate(process::StationaryMarkovChain, x0, T, nscen)
    scenarios = zeros(T+1, nscen)
    # Project initial point onto the Markov Chain
    ind0, _ = project(process, x0)
    for n in 1:nscen
        scenarios[1, n] = process.x[ind0]
        ind = ind0
        for t in 1:T
            probs = view(process.proba, ind, :)
            ind = _sample(probs)
            scenarios[t+1, n] = process.x[ind]
        end
    end
    return scenarios
end


abstract type DiscretizationAlgorithm end

struct TauchenHussey <: DiscretizationAlgorithm
    N::Int
end

function _gaussian_quadrature(law::Normal, N::Int)
    x0, w0 = FastGaussQuadrature.gausshermite(N)
    x = (sqrt(2.0) * law.σ) .* x0 .+ law.μ
    w = w0 ./ sqrt(pi)
    return (x, w)
end

function fit_markov(process::AR1, algorithm::TauchenHussey)
    N = algorithm.N
    μ, σ = process.mu, process.sigma
    law = Normal(μ, σ)
    x, w = _gaussian_quadrature(law, N)
    proba = zeros(N, N)
    for i in 1:N, j in 1:N
        prior = (1 - process.phi) * μ + process.phi * x[i]
        zj = x[j]
        proba[i, j] = w[j] * pdf(Normal(prior, σ), zj) / pdf(Normal(μ, σ), zj)
    end
    # Renormalization
    for i in 1:N
        proba[i, :] ./= sum(proba[i, :])
    end
    return StationaryMarkovChain(x, proba)
end

