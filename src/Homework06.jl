module Homework06

using Plots
using Base: product
using LinearAlgebra: norm
using NLsolve
using DifferentialEquations

export AbstractTVC, AbstractGrowthModel,
    TVC,
    kstart,
    kstop,
    tmax,
    tmin,
    ConstrGrowthModel,
    UnconstrGrowthModel,
    production,
    mglproduction,
    utilty,
    mglutility,
    steady_state_k,
    steady_state_psi,
    steady_state,
    Vconstr,
    dVconstr_dthat,
    dotk,
    dotpsi,
    constr_binds,
    dot_kpsi,
    dot_kpsi!,
    dotψ,
    constr_boundary,
    myquiver!,
    psi_nullcline,
    k_nullcline,
    khat,
    psihat,
    that_root,
    that_root!,
    find_that,
    that_residual,
    ODE_that_to_t0,
    ODE_that_to_T,
    consumption,
    switch_solution,
    solution,
    kψpath,
    ψpath,
    kpath,
    tgrid


abstract type AbstractTVC end
abstract type AbstractGrowthModel end

struct TVC <: AbstractTVC
    kstart::Float64
    kstop::Float64
    tmax::Float64
end

TVC(;kstart=1.15, kstop=0.6, tmax=10.0) = TVC(kstart,kstop,tmax)

kstart(x::AbstractTVC) = x.kstart
kstop(x::AbstractTVC) = x.kstop
tmax(x::AbstractTVC) = x.tmax
tmin(x::AbstractTVC) = zero(Float64)

struct ConstrGrowthModel <: AbstractGrowthModel
    alpha::Float64 # production
    rho::Float64   # discount
    depr::Float64  # depreciation
end

struct UnconstrGrowthModel <: AbstractGrowthModel
    alpha::Float64 # production
    rho::Float64   # discount
    depr::Float64  # depreciation
end

# outer constructor
(::Type{T})(;alpha=0.5, rho=0.1, depr=0.3) where {T<:AbstractGrowthModel} = T(alpha, rho, depr)


const GrowthModel = AbstractGrowthModel

#access fields of struct
alpha(m::GrowthModel) = m.alpha
rho(  m::GrowthModel) = m.rho
depr( m::GrowthModel) = m.depr
rho_depr(m::GrowthModel) = rho(m) + depr(m)

# because I changed the name of the function
@deprecate rho_delta(m) rho_depr(m)

# conversion
ConstrGrowthModel(m::UnconstrGrowthModel) = ConstrGrowthModel(alpha(m), rho(m), depr(m))
UnconstrGrowthModel(m::ConstrGrowthModel) = UnconstrGrowthModel(alpha(m), rho(m), depr(m))

# -------------------------
# pdxn + mgl utility
# -------------------------

production(m::GrowthModel,k) = k^alpha(m)
mglproduction(m::GrowthModel,k) = alpha(m)*k^(alpha(m)-1)

utility(m::GrowthModel, c) = log(c)
mglutility(m::GrowthModel, c) = 1/c

invmglutility(m::GrowthModel, psi) = 1/psi

# -------------------------
# steady states
# -------------------------

# steady states
function steady_state_k(m::GrowthModel)
    frac = alpha(m)/rho_depr(m)
    return frac^( 1/(1-alpha(m)) )
end
function steady_state_psi(m::GrowthModel)
    khat = steady_state_k(m)
    fk = production(m, khat)
    return 1/(fk - depr(m)*khat)
end
steady_state(m::GrowthModel) = (steady_state_k(m), steady_state_psi(m),)

# -------------------------
# scrap value
# -------------------------

# Value of being in constrained part
function Vconstr(x::ConstrGrowthModel, tvc::AbstractTVC, that)
    that <= tmax(tvc) || throw(DomainError(that))
    exp_rho_that = exp(-rho(x)*that)
    exp_rho_T = exp(-rho(x)*tmax(tvc))

    a = log(kstop(tvc)) + depr(x)*tmax(tvc) - depr(x)/rho(x)
    b = (exp_rho_that - exp_rho_T)/rho(x)

    c = depr(x)/rho(x)
    d = exp_rho_that*that - exp_rho_T*tmax(tvc)

    return a*b - c*d
end

# Value of being in constrained part
function dVconstr_dthat(x::ConstrGrowthModel, tvc::AbstractTVC, that::Number)
    that <= tmax(tvc) || throw(DomainError(that))
    exp_rho_that = exp(-rho(x)*that)
    a = log(kstop(tvc)) + depr(x)*tmax(tvc) - depr(x)/rho(x)
    b = depr(x)/rho(x)*(rho(x)*that - 1)
    return exp_rho_that*(-a + b)
end

# -------------------------
# ODEs
# -------------------------

# ODEs
dotk(  m::UnconstrGrowthModel, k, psi) = production(m,k) - invmglutility(m,psi) - depr(m)*k
dotpsi(m::UnconstrGrowthModel, k, psi) = (rho_depr(m) - mglproduction(m,k))*psi

consumption(m::UnconstrGrowthModel, k, psi) = invmglutility(m,psi)

constr_boundary(m::ConstrGrowthModel, k) = mglutility(m,production(m,k))
constr_binds(m::ConstrGrowthModel, k, psi) = (psi <= constr_boundary(m,k))

function dotk(m::ConstrGrowthModel, k, psi)
    if constr_binds(m,k,psi)
        return - depr(m)*k
    else
        m_unconstr = UnconstrGrowthModel(m)
        return dotk(m_unconstr, k, psi)
    end
end
function dotpsi(m::ConstrGrowthModel, k, psi)
    if constr_binds(m,k,psi)
        return rho_depr(m)*psi - mglproduction(m,k) * mglutility(m,production(m,k))
    else
        m_unconstr = UnconstrGrowthModel(m)
        return dotpsi(m_unconstr, k, psi)
    end
end

function consumption(m::ConstrGrowthModel, k, psi)
    if constr_binds(m,k,psi)
        return production(m,k)
    else
        m_unconstr = UnconstrGrowthModel(m)
        return consumption(m_unconstr, k, psi)
    end

end


dotψ(args...) = dotpsi(args...)

dot_kpsi(m::GrowthModel, y...) = (dotk(m,y...), dotpsi(m,y...))

# solvers in DifferentialEquations.jl require functions that
# take the form   `f!(dy, y, params, t)`
function dot_kpsi!(dy, y, m, t)
    dy .= dot_kpsi(m, y...)
end

# -------------------------
# plotting
# -------------------------

function myquiver!(plt, model::AbstractGrowthModel; scalek=1, scalepsi=1, kwargs...)

    ssk, ssψ = steady_state(model)

    # grid of points in K, I space
    kspace = range(ssk*0.5, stop = 1.2*ssk, length=15)
    ψspace = range(ssψ*0.5, stop = 1.2*ssψ, length=15)

    # create a "mesh"
    kψspace = product(kspace, ψspace)
    KK = zeros(length(kψspace))
    ψψ = similar(KK)
    for (n, (k,ψ)) in enumerate(kψspace)
        KK[n], ψψ[n] = k,ψ
    end

    # comptue ODEs for all values in mesh
    dotKK = broadcast(x -> dotk(model, x...)*scalek, zip(KK,ψψ))
    dotψψ = broadcast(x -> dotψ(model, x...)*scalepsi, zip(KK,ψψ))

    quiver!(plt, KK, ψψ, quiver=(dotKK,dotψψ); kwargs...)
end

k_nullcline(  m::UnconstrGrowthModel, k) = mglutility(m, production(m,k) - depr(m)*k)

function psi_nullcline(m::ConstrGrowthModel, k)
    k < steady_state_k(m) && throw(DomainError(k))
    return mglproduction(m,k) / (production(m,k) * rho_depr(m))
end

# -------------------------
# switchtime
# -------------------------

function khat(x::ConstrGrowthModel, tvc::TVC, that::Number)
    kT = kstop(tvc)
    tmax(tvc) <= that && @warn "T = $(tmax(tvc)) <= t hat = $that"
    T = tmax(tvc)
    m = depr(x)
    return kT*exp(m*(T - that))
end

function psihat(x::ConstrGrowthModel, tvc::TVC, that::Number)
    k = khat(x,tvc,that)
    fk = production(x,k)
    return mglutility(x, fk)
end

k_psi_hat(x, tvc, that) = [khat(x,tvc,that), psihat(x,tvc,that)]


function bc!(residual, u, model, tvc, t)
    residual[1] = first(u[end]) - kstart(tvc)  # k(0) = 15
    residual[2] = first(u[1]) - kstop(tvc)  # k(T) = 20
    return residual
end

function makeBVProblem(x::UnconstrGrowthModel, tvc, u0)
    tspan = (tmax(tvc), 0.0)
    bcinner!(res, u, model, t) = bc!(res, u, model, tvc, t)
    return BVProblem(dot_kpsi!, bcinner!, u0, tspan, x)
end

# -------------------------
# Define ODE problems
# -------------------------

function ODE_that_to_t0(x,tvc,that)
    u0 = k_psi_hat(x, tvc, that)
    tspan = (that, 0.0)
    prob = ODEProblem(dot_kpsi!, u0, tspan, x)
end

function ODE_that_to_T(x,tvc,that)
    u0 = k_psi_hat(x, tvc, that)
    tspan = (that,tmax(tvc))
    prob = ODEProblem(dot_kpsi!, u0, tspan, x)
end

# -------------------------
# Boundary conditions for \hat t
# -------------------------

function that_residual(x::ConstrGrowthModel, tvc::TVC, that)
    prob = ODE_that_to_t0(x,tvc,that)
    sol = solve(prob, Tsit5())
    return sol(0.0)[1] - kstart(tvc)
end

function that_residual!(res, x, tvc, that)
    res[1] = that_residual(x,tvc, that[1])
end

function find_that(x,tvc,that0; kwargs...)
    f(res,that) = that_residual!(res, x, tvc, that)
    nlsolve(f, that0; kwargs...)
end

# -------------------------
# Solution concept
# -------------------------

abstract type AbstractGrowthModelSolution end

struct ConstrainedGrowthModelSolution{M<:ConstrGrowthModel, T<:ODESolution} <: AbstractGrowthModelSolution
    model::M
    sol0::T
    sol1::T
    function ConstrainedGrowthModelSolution(model::M, sol0::T, sol1::T) where {M,T}
        maximum(sol0.t) == minimum(sol1.t) || throw(error("solutions in wrong order"))
        return new{M,T}(model,sol0,sol1)
    end
end

struct UnconstrainedGrowthModelSolution{M<:UnconstrGrowthModel, T<:ODESolution} <: AbstractGrowthModelSolution
    model::M
    sol0::T
    function UnconstrainedGrowthModelSolution(model::M, sol0::T) where {M,T}
        return new{M,T}(model,sol0)
    end
end

start(x::AbstractGrowthModelSolution) = minimum(x.sol0.t)
stop(x::ConstrainedGrowthModelSolution) = maximum(x.sol1.t)
stop(x::UnconstrainedGrowthModelSolution) = maximum(x.sol0.t)

function solution(x::ConstrGrowthModel,tvc,that)
    sol0 = solve(ODE_that_to_t0(x,tvc,that), Tsit5())
    sol1 = solve(ODE_that_to_T( x,tvc,that), Tsit5())
    return ConstrainedGrowthModelSolution(x, sol0, sol1)
end

function solution(x::UnconstrGrowthModel, tvc, uend; method=Shooting(Tsit5()), kwargs...)
    bvp = makeBVProblem(x, tvc, uend)
    sol = solve(bvp, method; kwargs...)
    return UnconstrainedGrowthModelSolution(x, sol)
end

t_in_sol(x::ODESolution, t) = (minimum(x.t) <= t <= maximum(x.t))
pick_solution(x::ConstrainedGrowthModelSolution, t) = t_in_sol(x.sol0, t) ? x.sol0 : x.sol1
pick_solution(x::UnconstrainedGrowthModelSolution, t) = x.sol0

kψpath(x::AbstractGrowthModelSolution, t) = pick_solution(x,t)(t)
kpath( x::AbstractGrowthModelSolution, t) = first(kψpath(x,t))
ψpath( x::AbstractGrowthModelSolution, t) = last(kψpath(x,t))
dotk(  x::AbstractGrowthModelSolution, t) = dotk(x.model, kψpath(x,t)...)
dotψ(  x::AbstractGrowthModelSolution, t) = dotψ(x.model, kψpath(x,t)...)
consumption(x::AbstractGrowthModelSolution, t) = consumption(x.model, kψpath(x,t)...)
get_that(x::AbstractGrowthModelSolution) = maximum(x.sol0.t)

kψpath(     x) = map(t -> kψpath(     x, t), tgrid(x))
kpath(      x) = map(t -> kpath(      x, t), tgrid(x))
ψpath(      x) = map(t -> ψpath(      x, t), tgrid(x))
dotk(       x) = map(t -> dotk(       x, t), tgrid(x))
dotψ(       x) = map(t -> dotψ(       x, t), tgrid(x))
consumption(x) = map(t -> consumption(x, t), tgrid(x))


function tgrid(x::AbstractGrowthModelSolution, n=101)
    return range(start(x); stop=stop(x), length=n)
end

function solutionplots(solutions::AbstractGrowthModelSolution)
    ts = tgrid(solutions) #  = range(0; stop=tmax(tvc), length=500)

    plt1 = plot(title="optimal path", legend=:left)
    plot!(ts, t -> kpath(solutions, t); label="\$k(t)\$")
    plot!(ts, t -> ψpath(solutions, t); label="\$\\psi(t)\$")
    plot!(ts, t -> consumption(solutions, t); label="\$c(t)\$")

    plt2 = plot(title="kdot, psidot", legend=:bottomleft)
    plot!(ts, t -> dotk(solutions,t); label="\$\\dot k(t)\$")
    plot!(ts, t -> dotψ(solutions,t); label="\$\\dot \\psi(t)\$")

    return plt1, plt2
end


# -------------------------
# Optimal switchtime
# -------------------------

function that_root(x, tvc, that)
    k = khat(x,tvc,that)
    ψ = psihat(x,tvc,that)
    m = depr(x)

    fk = production(x,k)
    u = utility(x, fk)

    dVdt = dVconstr_dthat(x, tvc, that)

    return u - ψ*m*k + dVdt
end


function that_root!(F, x, tvc, that)
    F[1] = that_root(x, tvc, that[1])
end


# find_that(x,tvc, that0::Vector=[1.0,]) = nlsolve((F,t) -> that_root!(F,x,tvc,t), that0)

#
# function bc!(residual, sol, p, t)
#     residual[1] = sol(pi/4)[1] + pi/2 # use the interpolation here, since indexing will be wrong for adaptive methods
#     residual[2] = sol(pi/2)[1] - pi/2
# end


end # module
