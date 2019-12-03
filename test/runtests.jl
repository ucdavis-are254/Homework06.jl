# so that when we change our module, it reloads
using Revise

using Test
using Calculus
using NLsolve
using Plots
using Base: product
using Homework06


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
myquiver!(plt_quiv, umodel; scalek=0.2, arrow=arrow(:closed), label="Unconstrained", alpha=0.1)
myquiver!(plt_quiv, cmodel; scalek=0.2, arrow=arrow(:closed), label="Constrained", alpha=0.1)

# steady state
ssk, ssψ = steady_state(umodel)
scatter!(plt_quiv, [ssk],[ssψ], label="SS")

# nullclines
plot!(1.2 : 0.1 : 2, k -> k_nullcline(umodel,k); label="\$\\dot k = 0\$")
vline!([ssk,]; label="\$\\dot \\psi = 0\$")
plot!(ssk : 0.1 : 2, k -> psi_nullcline(cmodel,k); label="\$\\dot \\psi = 0\$", linestyle=:dot)

# constraint boundary
plot!(1.2 : 0.1 : 2, k -> constr_boundary(cmodel,k); label="\$\\psi = u'(f(k))\$", linestyle=:dash)

# ----------------------------
# check derivatives for scrap value
# ----------------------------

@test dVconstr_dthat(cmodel, tvc, 5.0) ≈ Calculus.derivative(t -> Vconstr(cmodel, tvc, t), 5.0)

@show thatsol = find_that(cmodel, TVC(;tmax=40.0), [25.0,])
t_hat = thatsol.zero[1]

plot(20 : 0.1 : 40, t -> Homework06.that_root(cmodel, TVC(;tmax=40.0), t))
