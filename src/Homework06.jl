module Homework06

using Plots
using Base: product
using LinearAlgebra: norm
using NLsolve

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
    find_that

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

dotψ(args...) = dotpsi(args...)

# solvers in DifferentialEquations.jl require functions that
# take the form   `f!(dy, y, params, t)`
function dot_kpsi!(dy, y, m, t)
    k,i = y
    dy[1] = dotk(m, k, i)
    dy[2] = dotpsi(m, k, i)
    return dy
end

# -------------------------
# plotting
# -------------------------

function myquiver!(plt, model::AbstractGrowthModel; scalek=1, scalepsi=1, kwargs...)

    ssk, ssψ = steady_state(model)

    # grid of points in K, I space
    kspace = range(ssk*0.8, stop = 1.2*ssk, length=15)
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
    T = tmax(tvc)
    m = depr(x)
    return kT*exp(m*(T - that))
end

function psihat(x::ConstrGrowthModel, tvc::TVC, that::Number)
    k = khat(x,tvc,that)
    fk = production(x,k)
    return mglutility(x, fk)
end

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


find_that(x,tvc, that0::Vector=[1.0,]) = nlsolve((F,t) -> that_root!(F,x,tvc,t), that0)






end # module
