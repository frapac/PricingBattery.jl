
using SDDP
using HiGHS
using Statistics
using DelimitedFiles

include("markov_chain.jl")

function epex_data()
    epex_price = readdlm("logmeanreversion.txt")
    Δ = 24
    ndays = 365
    @assert length(epex_price) == Δ * ndays
    epex_ln = zeros(ndays)
    for i in 1:ndays
        k = (i-1)*Δ + 1
        l = i*Δ
        epex_ln[i] = mean(epex_price[k:l])
    end
    return epex_ln
end

function markov_sddp(
    markov_model,
    x0,
    horizon,
    average_price;
    iteration_limit = 25,
    linearize=true,
    rho=0.2,
    C_max=5.0,
    alpha=0.2,
    Q=0.0,
)
    r = 0.0001
    l = 0.05

    # create the Markovian graph
    ind0, x0d = project(markov_model, x0)
    graph = SDDP.Graph((0, ind0))
    # Add nodes
    nb_nodes = length(markov_model.x)
    for t in 1:horizon
        for k in 1:nb_nodes
            SDDP.add_node(graph, (t, k))
        end
    end

    # Initial stage
    for i in 1:nb_nodes
        xi = markov_model.x[i]
        pij = markov_model.proba[ind0, i]
        SDDP.add_edge(graph, (0, ind0) => (1, i), pij)
    end
    for t in 2:horizon
        for i in 1:nb_nodes, j in 1:nb_nodes
            xi, xj = markov_model.x[i], markov_model.x[j]
            pij = markov_model.proba[i, j]
            SDDP.add_edge(graph, (t-1, i) => (t, j), pij)
        end
    end

    # create the model
    model = SDDP.PolicyGraph(
        graph,
        sense = :Min,
        lower_bound = -1e4,
        optimizer = linearize ? HiGHS.Optimizer : Ipopt.Optimizer,
    ) do sp, node
        t, k = node
        price = 1.0
        # State
        @variable(sp, 0 <= C <= C_max, SDDP.State, initial_value = 0.0)
        @variable(sp, Z, SDDP.State, initial_value = -Q)
        # Control
        @variable(sp, -alpha * C_max <= U <= alpha * C_max)
        if t == horizon
            @variable(sp, 0 <= utility)
        end

        # Dynamics
        @constraint(sp, gain, Z.out == (1+r)*Z.in - price*U)
        @constraint(sp, charge, C.out == (1 - l)*C.in + U)

        if t == horizon
            if linearize
                zmin, zmax = 0, 10000
                # vmin, vmax = (zmin, zmax) .* rho
                for v in range(zmin, zmax, 5000)
                    exp_v = exp(-rho * v)
                    zk = v
                    @constraint(sp, utility >= exp_v / rho - exp_v *(Z.out - zk))
                end
                @stageobjective(sp, utility) # maximize the final amount of money of the investor
            else
                @NLconstraint(sp, utility == exp(-rho*(Z_plus.out - Z_moins.out))/rho)
                @stageobjective(sp, utility) # maximize the final amount of money of the investor
            end
        end

        SDDP.parameterize(sp, [exp(average_price[t+1] + markov_model.x[k])]) do ω
            JuMP.set_normalized_coefficient(gain, U,  ω)
        end

    end

    # train the model
    SDDP.train(model; iteration_limit)

    # print final wealth
    price = SDDP.calculate_bound(model)

    return model, price, graph
end

function _utility(z)
    rho = 0.2
    return 1.0 / rho * exp(- rho * z)
end

function simulate(model, nsimu)
    sim = SDDP.simulate(model, nsimu, [:Z])
    final_wealth = Float64[]
    for s in sim
        z = s[end][:Z]
        push!(final_wealth, z.out)
    end
    return final_wealth
end

function _next!(x, pb::SDDP.Node, price)
    m = pb.subproblem
    # Set current price in JuMP model
    JuMP.set_normalized_coefficient(m[:gain], m[:U],  price)

    for (k, state) in enumerate(keys(pb.states))
        JuMP.fix(pb.states[state].in, x[k])
    end

    JuMP.optimize!(m)

    for (k, state) in enumerate(keys(pb.states))
        x[k] = JuMP.value(pb.states[state].out)
    end
    return
end

function _out_of_sample(model, train_proc::StationaryMarkovChain, test_proc::AR1)
    T = 364 # number of timesteps
    ε = Normal(0.0, test_proc.sigma)
    moy_ln = epex_data()

    nx = length(model.initial_root_state)
    x0 = zeros(nx)

    Δp = 0.0
    for t in 1:T
        # Simulate AR
        Δp = (1.0 - test_proc.phi) * test_proc.mu + test_proc.phi * Δp + rand(ε)
        # Get current selling price
        price = exp(moy_ln[t+1] + Δp)
        # Get nearest point in Markov chain
        k, val = project(train_proc, Δp)
        # Find nearest model
        pb = model.nodes[(t, k)]
        _next!(x0, pb, price)
    end
    return x0[1]
end

function simulate(model, train_proc, test_proc, nsimu)
    final_gain = Float64[]
    for k in 1:nsimu
        push!(final_gain, _out_of_sample(model, train_proc, test_proc))
    end
    return final_gain
end

# Figure 1: SDDP's Convergence against size of
function compute_figure1(config, Nds)
    # Stochastic model
    price_process = AR1(config.stats.mu, config.stats.phi, config.stats.sigma)
    moy_ln = epex_data()
    horizon = config.model.horizon

    x0 = 0.0
    for Nd in Nds
        @info "Nd=$(Nd) rho=$(config.rho)"
        markov_model = fit_markov(price_process, TauchenHussey(Nd))
        model, price_lb, graph = markov_sddp(markov_model, x0, horizon, moy_ln; iteration_limit=config.sddp.sddp_it, rho=config.rho)
        lb = [log.bound for log in model.most_recent_training_results.log]
        writedlm(joinpath("results", "lower_bound_nd$(Nd).txt"), lb)
    end
end

function compute_figure2(config, Nds)
    x0 = 0.0
    # Stochastic model
    price_process = AR1(config.stats.mu, config.stats.phi, config.stats.sigma)
    moy_ln = epex_data()
    horizon = config.model.horizon

    for Nd in Nds
        @info "Nd=$(Nd) rho=$(config.rho)"
        markov_model = fit_markov(price_process, TauchenHussey(Nd))
        model, price_lb, graph = markov_sddp(markov_model, x0, horizon, moy_ln; iteration_limit=config.sddp.sddp_it, rho=config.rho)
        lb = [log.bound for log in model.most_recent_training_results.log]
        writedlm(joinpath("results", "lower_bound_nd$(Nd).txt"), lb)

        in_utility = simulate(model, config.sddp.nsimus)
        out_utility = simulate(model, markov_model, price_process, config.sddp.nsimus)

        writedlm(joinpath("results", "simulation_$(Nd)_insamp.txt"), in_utility)
        writedlm(joinpath("results", "simulation_$(Nd)_outsamp.txt"), out_utility)
    end
end

# sensitivity w.r.t. volatility
function compute_figure3(config, Nd, sigmas)
    x0 = 0.0
    horizon = config.model.horizon
    moy_ln = epex_data()
    results = zeros(length(sigmas), 4)

    k = 0
    for sigma in sigmas
        @info "sigma=$(sigma) Nd=$(Nd)"
        k += 1
        price_process = AR1(config.stats.mu, config.stats.phi, sigma)
        markov_model = fit_markov(price_process, TauchenHussey(Nd))
        model, price_lb, graph = markov_sddp(markov_model, x0, horizon, moy_ln; iteration_limit=config.sddp.sddp_it, rho=config.rho)

        out_utility = simulate(model, markov_model, price_process, config.sddp.nsimus)

        results[k, 1] = sigma
        results[k, 2] = Nd
        results[k, 3] = SDDP.calculate_bound(model)
        results[k, 4] = mean(out_utility)
    end

    writedlm(joinpath("results", "sensitivity_sigmas.txt"), results)
end

# sensitivity w.r.t. risk aversion
function compute_figure4(config, Nd, rhos)
    x0 = 0.0
    horizon = config.model.horizon
    moy_ln = epex_data()
    results = zeros(length(rhos), 4)

    price_process = AR1(config.stats.mu, config.stats.phi, config.stats.sigma)
    markov_model = fit_markov(price_process, TauchenHussey(Nd))

    k = 0
    for rho in rhos
        @info "rho=$(rho) Nd=$(Nd)"
        k += 1
        model, price_lb, graph = markov_sddp(markov_model, x0, horizon, moy_ln; iteration_limit=config.sddp.sddp_it, rho=rho)

        out_utility = simulate(model, markov_model, price_process, config.sddp.nsimus)

        results[k, 1] = rho
        results[k, 2] = Nd
        results[k, 3] = SDDP.calculate_bound(model)
        results[k, 4] = mean(out_utility)
    end

    writedlm(joinpath("results", "sensitivity_rhos.txt"), results)
end

# sensitivity w.r.t. charge/discharge rate
function compute_figure5(config, Nd, alphas)
    x0 = 0.0
    horizon = config.model.horizon
    moy_ln = epex_data()
    results = zeros(length(rhos), 4)

    price_process = AR1(config.stats.mu, config.stats.phi, config.stats.sigma)
    markov_model = fit_markov(price_process, TauchenHussey(Nd))

    k = 0
    for alpha in alphas
        @info "rho=$(rho) Nd=$(Nd)"
        k += 1
        model, price_lb, graph = markov_sddp(
            markov_model,
            x0,
            horizon,
            moy_ln;
            iteration_limit=config.sddp.sddp_it,
            rho=config.rho,
            alpha=alpha,
        )

        out_utility = simulate(model, markov_model, price_process, config.sddp.nsimus)

        results[k, 1] = alpha
        results[k, 2] = Nd
        results[k, 3] = SDDP.calculate_bound(model)
        results[k, 4] = mean(out_utility)
    end

    writedlm(joinpath("results", "sensitivity_charging_rate.txt"), results)
end

