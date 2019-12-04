# so that when we change our module, it reloads
using Revise

using Test
using NLsolve
using Plots
using Base: product
using Homework06
using DifferentialEquations
using BoundaryValueDiffEq

cmodel = ConstrGrowthModel()
umodel = UnconstrGrowthModel()
tvc = TVC()

@test UnconstrGrowthModel(ConstrGrowthModel(umodel)) == umodel
@test ConstrGrowthModel(UnconstrGrowthModel(cmodel)) == cmodel

# ----------------------------
# make plot
# ----------------------------

# create a new plot
plt_quiv = plot(;legend=:topright, xlabel="\$k\$", ylabel="\$\\psi\$")

# add quivers
myquiver!(plt_quiv, umodel; scalek=0.2, scalepsi=1.0, arrow=arrow(:closed), label="Unconstrained", alpha=0.1)
myquiver!(plt_quiv, cmodel; scalek=0.2, scalepsi=1.0, arrow=arrow(:closed), label="Constrained", alpha=0.1)

# steady state
ssk, ssψ = steady_state(umodel)
scatter!(plt_quiv, [ssk],[ssψ], label="") # label="SS")

# nullclines
plot!(kstop(TVC()) : 0.1 : 2, k -> k_nullcline(umodel,k); label="") # "\$\\dot k = 0\$")
vline!([ssk,]; label="") # "\$\\dot \\psi = 0\$")
plot!(ssk : 0.1 : 2, k -> psi_nullcline(cmodel,k); label="") # "\$\\dot \\psi = 0\$", linestyle=:dot)

# constraint boundary
plot!(kstop(TVC()) : 0.1 : 2, k -> constr_boundary(cmodel,k); label="") # "\$\\psi = u'(f(k))\$", linestyle=:dash)

# ----------------------------
# solve unconstrained bvp
# ----------------------------

# solve for \hat t
tvc1 = TVC(;tmax=2.5)
tvc2 = TVC(;tmax=7.0)
tvc3 = TVC(;tmax=20.0)

sol1u = solution(umodel, tvc1, [0.6, 0.88])
sol2u = solution(umodel, tvc2, [0.6, 0.88])
sol3u = solution(umodel, tvc3, [0.6, 0.8831])
plot(Homework06.solutionplots(sol1u)...)
plot(Homework06.solutionplots(sol2u)...)
plot(Homework06.solutionplots(sol3u)...)

# ----------------------------
# solve constrained bvp
# ----------------------------

# -------
# TVC 1
# -------

# solve constrained problem for \hat t
tvc = tvc1
sol = sol1u
plot(0 : 0.1 : tmax(tvc)-0.1, t -> that_residual(cmodel,tvc,t), legend=:none, title = "\$k(0) - k_0\$")

that_solution = find_that(cmodel, tvc, [1.0,]; show_trace=true)
that = that_solution.zero[1]
sol1c = solution(cmodel, tvc, that)
plt1, plt2 = Homework06.solutionplots(sol1c)

plot!(plt1, tgrid(sol), t -> kpath(sol, t); label="\$k^u(t)\$", linestyle=:dot, color=:blue)
plot!(plt1, tgrid(sol), t -> ψpath(sol, t); label="\$\\psi^u(t)\$", linestyle=:dot, color=:orange)
plot!(plt1, tgrid(sol), t -> consumption(sol, t); label="\$c^u(t)\$", linestyle=:dot, color=:green)
plot!(plt2, tgrid(sol), t -> dotk(sol, t); label="\$\\dot k^u(t)\$", linestyle=:dot, color=:blue)
plot!(plt2, tgrid(sol), t -> dotψ(sol, t); label="\$\\dot \\psi^u(t)\$", linestyle=:dot, color=:orange)
# vline!(plt1, [that,]; label="", linestyle = :dash)
# vline!(plt2, [that,]; label="", linestyle = :dash)
plot(plt1, plt2)

# -------
# TVC 2
# -------

# solve constrained problem for \hat t
tvc = tvc2
sol = sol2u
plot(0 : 0.1 : tmax(tvc)-2, t -> that_residual(cmodel,tvc,t), legend=:none, title = "\$k(0) - k_0\$")

that_solution = find_that(cmodel, tvc, [tmax(tvc)-2,]; show_trace=true)
that = that_solution.zero[1]
sol2c = solution(cmodel, tvc, that)

plt1, plt2 = Homework06.solutionplots(sol2c)

plot!(plt1, tgrid(sol), t -> kpath(sol, t); label="\$k^u(t)\$", linestyle=:dot, color=:blue)
plot!(plt1, tgrid(sol), t -> ψpath(sol, t); label="\$\\psi^u(t)\$", linestyle=:dot, color=:orange)
plot!(plt1, tgrid(sol), t -> consumption(sol, t); label="\$c^u(t)\$", linestyle=:dot, color=:green)
plot!(plt2, tgrid(sol), t -> dotk(sol, t); label="\$\\dot k^u(t)\$", linestyle=:dot, color=:blue)
plot!(plt2, tgrid(sol), t -> dotψ(sol, t); label="\$\\dot \\psi^u(t)\$", linestyle=:dot, color=:orange)
# vline!(plt1, [that,]; label="", linestyle = :dash)
# vline!(plt2, [that,]; label="", linestyle = :dash)
plot(plt1, plt2)

# -------
# TVC 3
# -------

# solve constrained problem for \hat t
tvc = tvc3
sol = sol3u
plot(18 : 0.05 : tmax(tvc)-1.501, t -> that_residual(cmodel,tvc,t), legend=:none, title = "\$k(0) - k_0\$", ylims=(0,Inf))

that_solution = find_that(cmodel, tvc, [18.44,]; show_trace=true)
that = that_solution.zero[1]
sol3c = solution(cmodel, tvc, that)

plt1, plt2 = Homework06.solutionplots(sol3c)

plot!(plt1, tgrid(sol), t -> kpath(sol, t); label="\$k^u(t)\$", linestyle=:dot, color=:blue)
plot!(plt1, tgrid(sol), t -> ψpath(sol, t); label="\$\\psi^u(t)\$", linestyle=:dot, color=:orange)
plot!(plt1, tgrid(sol), t -> consumption(sol, t); label="\$c^u(t)\$", linestyle=:dot, color=:green)
plot!(plt2, tgrid(sol), t -> dotk(sol, t); label="\$\\dot k^u(t)\$", linestyle=:dot, color=:blue)
plot!(plt2, tgrid(sol), t -> dotψ(sol, t); label="\$\\dot \\psi^u(t)\$", linestyle=:dot, color=:orange)
# vline!(plt1, [that,]; label="", linestyle = :dash)
# vline!(plt2, [that,]; label="", linestyle = :dash)
plot(plt1, plt2)

# ----------------------------
# solve constrained bvp
# ----------------------------

p = deepcopy(plt_quiv)
plot!(p, kpath(sol1c), ψpath(sol1c), label = "", linestyle=:dash) # "\$T=$(tmax(tvc1))\$")
plot!(p, kpath(sol1u), ψpath(sol1u), label = "") # "\$T=$(tmax(tvc1))\$")

p = deepcopy(plt_quiv)
plot!(p, kpath(sol2c), ψpath(sol2c), label = "", linestyle=:dash) # "\$T=$(tmax(tvc1))\$")
plot!(p, kpath(sol2u), ψpath(sol2u), label = "") # "\$T=$(tmax(tvc1))\$")

p = deepcopy(plt_quiv)
plot!(p, kpath(sol3c), ψpath(sol3c), label = "", linestyle=:dash) # "\$T=$(tmax(tvc1))\$")
plot!(p, kpath(sol3u), ψpath(sol3u), label = "") # "\$T=$(tmax(tvc1))\$")

# ----------------------------
# About the free end-time TVC
# ----------------------------

# create a new plot
plt_quiv2 = plot(;legend=:topright, xlabel="\$k\$", ylabel="\$\\psi\$")

# add quivers
myquiver!(plt_quiv2, umodel; arrow=arrow(:closed), label="Unconstrained", alpha=0.1)

kTs = 0.5 : 0.05 : 1.5
psiTs = [
    nlsolve( (res, psiT) -> Homework06.hamiltonian!(res, umodel, TVC(;kstop=k), psiT), [1.0,]).zero[1]
    for k in kTs
]

# steady state
ssk, ssψ = steady_state(umodel)
scatter!([ssk],[ssψ], label="") # label="SS")
plot!(kstop(TVC()) : 0.1 : 2, k -> k_nullcline(umodel,k); label="") # "\$\\dot k = 0\$")
vline!([ssk,]; label="") # "\$\\dot \\psi = 0\$")
plot!(kpath(sol3u), ψpath(sol3u), label = "") # "\$T=$(tmax(tvc1))\$")
scatter!(kTs, psiTs, label="\$\\mathcal H = 0\$")

##
