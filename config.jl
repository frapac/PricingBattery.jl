
struct ARParams
    mu::Float64
    phi::Float64
    sigma::Float64
end

struct ModelParams
    rho::Float64
    sc::Float64
    Q::Float64
    C_max::Float64
end

struct SDDPParams
    sddp_it::Int
    nsimus::Int
end

struct Config
    stats::ARParams
    model::ModelParams
    sddp::SDDPParams
    Nd::Vector{Int}
    sigmas::AbstractRange{Float64}
end

config = Config(
    ARParams(0.0, 0.54, 0.39),
    ModelParams(0.2, 1e5, 0.0, 5.0),
    SDDPParams(50, 1000),
    [1, 2, 4, 8, 16, 32],
    0.39:0.06:0.8,
)

