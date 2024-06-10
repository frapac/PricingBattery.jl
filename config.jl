
struct ARParams
    mu::Float64
    phi::Float64
    sigma::Float64
end

struct ModelParams
    horizon::Int
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
    rho::Float64
end

