using JuMP
using CPLEX

include("lecture.jl")

println(Data.Q)
println(Data.S)

# Data = donnees(3, 3, 2, 1)
# Data.Q = [[1, 0],
#           [0, 1],
#           [1, 0]]
# Data.S = [[1, 0, 0],
#           [0, 1, 0],
#           [0, 0, 1]]

P = 1
capa = Vector{Int}([12, 12, 12, 12, 12])
FO = 1:25
SO = 26:50

Data.P = P
Data.Capa = capa
Data.FO = FO
Data.SO = SO
sort!(Data.Capa, rev=true)


m = Model(CPLEX.Optimizer)

@variable(m, alpha >= 0)
@variable(m, x[1:Data.O, 1:Data.R], Bin)
@variable(m, z[1:Data.R, 1:Data.R], Bin)
@variable(m, y[1:Data.O, 1:Data.N, 1:Data.R] >= 0)
@variable(m, u[1:Data.O, 1:Data.R], Bin)
@variable(m, v[1:Data.O, 1:Data.R], Bin)

@constraint(m, [o in 1:Data.O], sum(x[o, t] for t in 1:Data.R) >= 1)

@constraint(m, [o in 1:Data.O, t in 1:Data.R], x[o, t] <= u[o, t])
@constraint(m, [o in 1:Data.O, t in 1:Data.R], x[o, t] <= 1 - v[o, t])
@constraint(m, [o in 1:Data.O, t in 1:Data.R], x[o, t] >= u[o, t] - v[o, t])

@constraint(m, [o in 1:Data.O, t in 2:Data.R], v[o, t-1] <= v[o, t])
@constraint(m, [o in 1:Data.O, t in 2:Data.R], u[o, t-1] <= u[o, t])

@constraint(m, [o in 1:Data.O, t in 1:Data.R, i in 1:Data.N], y[o, i, t] <= x[o, t] * Data.Q[i][o])
@constraint(m, [o in 1:Data.O, i in 1:Data.N], sum(y[o, i, t] for t in 1:Data.R) >= Data.Q[i][o])
@constraint(m, [t in 1:Data.R, i in 1:Data.N], sum(y[o, i, t] for o in 1:Data.O) <= sum(z[r, t] * Data.S[i][r] for r in 1:Data.R))

@constraint(m, [r in 1:Data.R], sum(z[r, t] for t in 1:Data.R) == 1)
@constraint(m, [t in 1:Data.R], sum(z[r, t] for r in 1:Data.R) == 1)

@constraint(m, [t in 1:Data.R], alpha >= sum(x[o, t] for o in 1:Data.O))

@objective(m, Min, alpha)

optimize!(m)

feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
println(feasible_solution_found)
if feasible_solution_found
    v_obj = JuMP.objective_value(m)
    println("Valeur fonction objective : ", v_obj)
end
