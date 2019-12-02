using DifferentialEquations, Plots, NamedTuples, Optim

# push!(LOAD_PATH, joinpath(pwd(), "julia"))
# using anatomyFunctionsSimple

gr()

#-------------- search / drilling ------------------------

θsl = @NT(
    J = 1, dF = [1.0,], ϕ = 2.0,
    m = 0.5, τ = 0.3, κ = 0.5, α = 0.75, β = 0.25, ubar = 0.1,
    n = 1.0, ρ = 0.15, p = 5.0
)


tspan = (100.0, 0.0)
x0 = steadystate(θsl; aT = 77.0)
x0[3] -= 1e-6
# x0[4] -= 1e-6
x0[end-1] += 1e-6
x0[end] -= 1e-6

τ0 = 18.4
τlims = (17.0, 19.0)

search = true

if search
    obj(xx) = norm(solve_search_drill(xx,x0,tspan,θsl)[1:2,end])
    @show res = Optim.optimize(obj, τlims...)
    τ1 = res.minimizer
else
    τ1 = τ0
end

sol = solve_search_drill(τ1, x0, tspan, θsl)
dx = dotx_search_drill(sol, τ1, θsl)

#----------- plots -----------

tt = sol.t .- minimum(sol.t)

plot(
    plot(tt, sol[[1,2],:]', label=["a","w"]),
    plot(tt, sol[5,:], label=["U"]),
    plot(tt, sol[1,:] .- sol[2,:], label=["a-w"]),
    plot(tt, sol[[3,4,6],:]', label=["\\psi_a","\\psi_w","Vu"]),
    # plot(tt, sol[3,:] .- sol[end,:], label="s"),
    # plot(tt, v, label="v", legend=:left),
    layout= grid(4,1)
)



plot(
    plot(tt, dx[1:2,:]', label=["da", "dw"], ylims=(-0.1,Inf)),
    plot(tt, dx[[3,6],:]', label=["d\\psi_a","Vu"], ylims=(-0.1,Inf)),
    layout = grid(2,1)
)



# tmβ = 2 - θsl.β
# mτ_κ = θsl.m * (1-θsl.τ)/θsl.κ
# v = θsl.m .* sol[end-1,:] .^ (θsl.α/tmβ) .* (mτ_κ .* θsl.n .* (sol[3,:] .- sol[end,:]) ) .^ (1/tmβ)











#
