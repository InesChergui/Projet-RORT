using JuMP
using CPLEX

include("lecture.jl")





Data.P = 5
Data.Capa = Vector{Int}([12,12,12,12,12])

Data.FO = [1,2,3,4,5]
Data.SO = [6,7,8,9,10]

m = Model(CPLEX.Optimizer)

@variable(m, x[1:Data.P, 1:Data.O], Bin)
@variable(m, y[1:Data.P, 1:Data.R], Bin)

@constraint(m, [r in 1:Data.R], sum( y[p, r] for p in 1:Data.P ) <= 1 )
@constraint(m, [o in Data.FO], sum( x[p, o] for p in 1:Data.P) == 1 )
@constraint(m, [o in Data.SO], sum( x[p, o] for p in 1:Data.P) <= 1 )
@constraint(m, [p in 1:Data.P], sum(x[p,o] for o in 1:O) <= Data.Capa[p])
@constraint(m, [i in 1:Data.N, p in 1:Data.P], sum(Data.Q[i][o]*x[p,o] for o in 1:O) <= sum(Data.S[i][r]*y[p,r] for r in 1:R))

# @objective(m, Min, sum( y[p,r] for p in 1:Data.P for r in 1:Data.R))
@objective(m, Min, sum( y[p,r] for p in 1:Data.P for r in 1:Data.R)*(size(Data.SO)[1]+1) - sum( x[p,o] for p in 1:Data.P for o in Data.SO))


optimize!(m)

global feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
println(feasible_solution_found)
global isOpt = termination_status(m) == MOI.OPTIMAL
if feasible_solution_found
    global v_obj = JuMP.objective_value(m)
    global vy = value.(y)
    global vx = value.(x)

    println("Valeur fonction objective : ", v_obj)
    println("valeur x : ", vx)
    println("valeur y : ", vy)
    
end
