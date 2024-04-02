
using SDDP
using HiGHS
using Statistics
using DelimitedFiles

include("markov_chain.jl")

function epex_data()
    data = readdlm("weatherdata.csv", ','; header=true)[1]
    epex_price = data[:, end]
    epex_ln = zeros(365)
    for i in 1:365
        k = (i-1)*96 + 1
        l = i*96
        epex_ln[i] = log(mean(epex_price[k:min(l, size(data, 1))]))
    end
    return epex_ln
end

function markov_sddp(
    markov_model,
    x0,
    average_price;
    iteration_limit = 25,
    linearize=true,
    rho=0.2,
)
    T = 364
    Q = 0 # initial amount of money invested
    C_max = 5 # capacity of the storage
    r = 0.0001
    l = 0.05
    sc = 1e5 #renormalization coefficient

    # create the Markovian graph
    ind0, x0d = project(markov_model, x0)
    graph = SDDP.Graph((0, ind0))
    # Add nodes
    nb_nodes = length(markov_model.x)
    for t in 1:T
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
    for t in 2:T
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
        lower_bound = -1e6,
        optimizer = linearize ? HiGHS.Optimizer : Ipopt.Optimizer,
    ) do sp, node
        t, k = node
        price = 1.0
        # State
        @variable(sp, 0 <= C <= C_max, SDDP.State, initial_value = 0.0)
        @variable(sp, 0 <= Z_plus, SDDP.State, initial_value = 0.0)
        @variable(sp, 0 <= Z_moins, SDDP.State, initial_value = Q/sc)
        # Control
        @variable(sp, 0 <= U_plus <=  C_max)
        @variable(sp, 0 <= U_moins <=   C_max)
        if t == T
            @variable(sp, 0 <= utility)
        end

        # Dynamics
        @constraint(sp, gain, Z_plus.out - Z_moins.out == (1+r)*Z_plus.in - (1+r)*Z_moins.in  - (price)*U_plus/sc + (price)*U_moins/sc)
        @constraint(sp, charge, C.out == (1 - l)*C.in + U_plus - U_moins)

        if t == T
            if linearize
                zmin, zmax = 0, 10000
                vmin, vmax = (zmin, zmax) .* (rho / sc)
                for v in range(vmin, vmax, 5000)
                    exp_v = exp(-v)
                    zk = v / rho
                    @constraint(sp, utility >= exp_v / rho - exp_v *(Z_plus.out - Z_moins.out - zk))
                end
                @stageobjective(sp, utility) # maximize the final amount of money of the investor
            else
                @NLconstraint(sp, utility == exp(-rho*(Z_plus.out - Z_moins.out))/rho)
                @stageobjective(sp, utility) # maximize the final amount of money of the investor
            end
        end

        SDDP.parameterize(sp, [exp(average_price[t+1] + markov_model.x[k])]) do ω
            JuMP.set_normalized_coefficient(gain, U_plus,  ω / sc)
            JuMP.set_normalized_coefficient(gain, U_moins, -ω / sc)
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
    sc = 10e4 #renormalization coefficient
    sim = SDDP.simulate(model, nsimu, [:Z_plus, :Z_moins])
    utilities = Float64[]
    for s in sim
        z_plus, z_moins = s[end][:Z_plus], s[end][:Z_moins]
        push!(utilities, _utility(z_plus.out - z_moins.out))
    end
    return utilities
end

function _next!(x, pb::SDDP.Node, price)
    m = pb.subproblem
    sc = 1e5
    # Set current price in JuMP model
    JuMP.set_normalized_coefficient(m[:gain], m[:U_plus],  price / sc)
    JuMP.set_normalized_coefficient(m[:gain], m[:U_moins], -price / sc)

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
    sc = 1e5
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
    z_moins, z_plus = x0[1], x0[2]
    return _utility(z_plus - z_moins)
end

function simulate(model, train_proc, test_proc, nsimu)
    final_gain = Float64[]
    for k in 1:nsimu
        push!(final_gain, _out_of_sample(model, train_proc, test_proc))
    end
    return final_gain
end

function convert_gain(val)
    sc, rho, r, tf = 1e5, 0.2, 1e-4, 364
    return -sc / rho * log(rho*val) / (1 + r)^tf
end

# Figure 1: SDDP's Convergence against size of
function compute_figure1(; total_it=50)
    phi, sigma = 0.54, 0.39
    # Stochastic model
    price_process = AR1(0.0, phi, sigma)

    moy_ln = epex_data()

    x0 = 0.0
    for Nd in [2, 5, 10, 20]
        @info "Nd=$(Nd)"
        markov_model = fit_markov(price_process, TauchenHussey(Nd))
        model, price_lb, graph = markov_sddp(markov_model, x0, moy_ln; iteration_limit=total_it)
        lb = [log.bound for log in model.most_recent_training_results.log]
        writedlm(joinpath("results", "lower_bound_nd$(Nd).txt"), lb)
    end
end

function compute_figure2(; total_it=50, nsimus=1000)
    phi, sigma = 0.54, 0.39
    x0 = 0.0
    # Stochastic model
    price_process = AR1(0.0, phi, sigma)
    moy_ln = epex_data()

    for Nd in [2, 5, 10, 20]
        @info "Nd=$(Nd)"
        markov_model = fit_markov(price_process, TauchenHussey(Nd))
        model, price_lb, graph = markov_sddp(markov_model, x0, moy_ln; iteration_limit=total_it)

        in_utility = simulate(model, nsimus)
        out_utility = simulate(model, markov_model, price_process, nsimus)

        writedlm(joinpath("results", "simulation_$(Nd)_insamp.txt"), in_utility)
        writedlm(joinpath("results", "simulation_$(Nd)_outsamp.txt"), out_utility)
    end
end

function compute_figure3(; total_it=50, nsimus=1000)
    phi = 0.54
    x0 = 0.0

    sigmas = 0.39:0.06:0.8
    Nds = [1, 2, 5, 10, 20]

    results = zeros(length(sigmas) * length(Nds), 4)

    k = 0
    for sigma in sigmas, Nd in Nds
        @info "sigma=$(sigma) Nd=$(Nd)"
        k += 1
        price_process = AR1(0.0, phi, sigma)
        markov_model = fit_markov(price_process, TauchenHussey(Nd))
        model, price_lb, graph = markov_sddp(markov_model, x0, moy_ln; iteration_limit=total_it)

        out_utility = simulate(model, markov_model, price_process, nsimus)

        results[k, 1] = sigma
        results[k, 2] = Nd
        results[k, 3] = SDDP.calculate_bound(model)
        results[k, 4] = mean(out_utility)
    end

    writedlm(joinpath("results", "sensitivity_std_det.txt"), results)
end

