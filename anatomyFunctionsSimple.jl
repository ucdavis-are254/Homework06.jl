module anatomyFunctionsSimple

using NamedTuples, DifferentialEquations
using Roots

export search_drill_f!, VuTresid, steadystate, condition_terminal, elapsed_time, elapsed_callback

#----------------- ODEs ---------------


function lambda(U::Real,dA_U::Real,da::Real,s::Real,ψw::Real,θ::NamedTuple)
    (ϕ, m, τ, κ, α, β, ubar, n, ρ, p) = (θ.ϕ, θ.m, θ.τ, θ.κ, θ.α, θ.β, θ.ubar, θ.n, θ.ρ, θ.p)

    tmp = 2α*ubar/U + dA_U*(τ*β - 2α) + ρ*β
    λ = s / ( (2-β)*s + ϕ*da*β ) * ( ϕ*da*tmp - (2-β)*ρ*ψw )

    return λ
end


function search_drill_f!(dx::AbstractArray, x::AbstractArray, regime::Integer, t::Real, θ::NamedTuple)

    (ϕ, m, τ, κ, α, β, ubar, n, ρ, p) = (θ.ϕ, θ.m, θ.τ, θ.κ, θ.α, θ.β, θ.ubar, θ.n, θ.ρ, θ.p)

    (a, w, ψa, ψw, U, Vu)  = (x...)

    # match surpluses
    s = ψa - Vu

    # a couple things we'll re-use
    tmβ = 2.0-β
    mτ_κ = m * (1-τ) / κ

    @assert U >= 0
    @assert s >= 0

    dA =   m * U^(2α/tmβ  ) * (mτ_κ * n * s)^(β/tmβ)
    dA_U = m * U^(2α/tmβ-1) * (mτ_κ * n * s)^(β/tmβ)

    # ̇a and ̇w
    dx[1] = (da = dA / n)
    dx[2] = (dw = (p + ψw) / ϕ)

    #costates
    λ = regime == 0 ? lambda(U,dA_U,da,s,ψw,θ) : 0.0
    dx[3] = -λ + ρ*ψa
    dx[4] =  λ + ρ*ψw

    # ̇U and ̇Vᵘ
    dx[5] = ubar - dA             # dU
    dx[6] = -τ * dA_U * s + ρ*Vu  # dVu

    return dx
end

search_drill_f(x, regime, t, θ) = search_drill_f!(similar(x), x, regime, t, θ)

#--------------- steady state functions ----------------


function VuTresid(Vu::Real, θ::NamedTuple)
    (ϕ, m, τ, κ, α, β, ubar, n, ρ, p) = (θ.ϕ, θ.m, θ.τ, θ.κ, θ.α, θ.β, θ.ubar, θ.n, θ.ρ, θ.p)

    ψaT = p - ϕ*ubar/n
    s = ψaT - Vu

    mτ_κ = m * (1-τ) / κ
    tα = 2*α
    βtα = β / tα
    tmβ = 2-β

    return -τ * ubar * (m/ubar) ^ (tmβ/tα) * (mτ_κ*n)^βtα * s^(1+βtα) + ρ*Vu

end

function steadystate(θ::NamedTuple; aT::Real=200.0)

    (ϕ, m, τ, κ, α, β, ubar, n, ρ, p) = (θ.ϕ, θ.m, θ.τ, θ.κ, θ.α, θ.β, θ.ubar, θ.n, θ.ρ, θ.p)
    @assert ubar > 0.0

    ψaT = p - ϕ*ubar/n
    ψwT = -ψaT

    VuT = find_zero((z) -> VuTresid(z,θ), (0.0, ψaT-1e-10))

    @assert VuT < ψaT

    mτ_κ = m * (1-τ) / κ
    tα = 2*α
    βtα = β / tα
    tmβ = 2-β

    sT = ψaT - VuT
    UT = (ubar/m)^( tmβ/tα ) * (mτ_κ*n*sT)^( -βtα )

    @assert UT > 0.0
    @assert lambda(UT,ubar/UT,ubar,sT,ψwT,θ) ≈ ρ*ψaT

    return [aT, aT, ψaT, ψwT, UT, VuT]

end

#--------------- callback functions ----------------

function condition_terminal(u,t,i)
    minacres = minimum(u[[1,2,5]])
    # s = u[3] - u[end]
    # return min(minacres, s)
    return minacres
end


elapsed_time(t::Real,tstart::Real,tdir::Real) = (t - tstart)*tdir
elapsed_time(t::Real,tspan::NTuple{2,<:Real}) = elapsed_time(t, tspan[1], tspan[1] > tspan[2] ? -1 : 1 )
elapsed_time(t::Real,i::AbstractODEIntegrator) = elapsed_time(t,i.sol.prob.tspan[1],i.tdir)

function elapsed_callback(τ::Real)
    condition(u,t,i) = elapsed_time(t,i) - τ
    affect!(i) = i.p += 1
    return ContinuousCallback(condition, affect!; save_positions=(true,true) )
end


#--------------- solve wrapper ----------------

function solve_search_drill(τ::Real, x0::AbstractArray, tspan::NTuple{2,<:Real}, θ::NamedTuple)
    cbs = CallbackSet(
        ContinuousCallback(condition_terminal, terminate!; save_positions=(true,true)),
        elapsed_callback(τ)
    )

    prob = ODEProblem((dx,x,p,t) -> search_drill_f!(dx,x,p,t,θ), x0, tspan, 0)
    sol = solve(prob, Tsit5(); callback=cbs, reltol=1e-8, abstol=1e-8)

    return sol
end

#----------------- dot -----------------------------

function dotx_search_drill(sol::AbstractODESolution, τ::Real, θ::NamedTuple)

    tspan = sol.prob.tspan
    dx = zeros(size(sol))

    for ti in 1:size(dx,2)
        t = sol.t[ti]
        regime = elapsed_time(t, tspan) > τ ? 1 : 0
        @views search_drill_f!(dx[:,ti], sol[:,ti], regime, t, θ)
    end

    return dx
end



#--------------- additional fct ----------------


function dsu(u,s,θ::NamedTuple)

    (ϕ, m, τ, κ, α, β, ubar, n, ρ, p) = (θ.ϕ, θ.m, θ.τ, θ.κ, θ.α, θ.β, θ.ubar, θ.n, θ.ρ, θ.p)

    mτ_κ = m * (1-τ) / κ
    tα = 2*α
    βtα = β / tα
    tmβ = 2-β


    da   = m * u^(α/tmβ)   * (mτ_κ * n * s)^(β/tmβ)
    da_u = m * u^(α/tmβ-1) * (mτ_κ * n * s)^(β/tmβ)

    du = ubar - da
    ds = s * (τ * da_u + ρ)

    return (du,ds)

end















end
