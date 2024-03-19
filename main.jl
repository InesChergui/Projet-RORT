using JuMP
using CPLEX
using Plots

include("lecture.jl")
include("greedy.jl")


function PLNE()
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, x[1:Data.P, 1:Data.O], Bin)
    @variable(m, y[1:Data.P, 1:Data.R], Bin)

    @constraint(m, [r in 1:Data.R], sum( y[p, r] for p in 1:Data.P ) <= 1 )
    @constraint(m, [o in Data.FO], sum( x[p, o] for p in 1:Data.P) == 1 )
    @constraint(m, [o in Data.SO], sum( x[p, o] for p in 1:Data.P) <= 1 )
    @constraint(m, [p in 1:Data.P], sum(x[p, o] for o in 1:O) <= Data.Capa[p])
    @constraint(m, [i in 1:Data.N, p in 1:Data.P], sum(Data.Q[i][o]*x[p, o] for o in 1:O) <= sum(Data.S[i][r]*y[p, r] for r in 1:R))

    @objective(m, Min, sum(y[p, r] for p in 1:Data.P, r in 1:Data.R)*(length(Data.SO) + 1) - sum(x[p, o] for p in 1:Data.P, o in Data.SO))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    println(feasible_solution_found)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        println("Valeur fonction objective : ", v_obj)
    end
end


function slaveProblemOrders(alpha)
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, x[1:Data.P, 1:Data.O], Bin)
    
    @constraint(m, [o in Data.FO], sum(x[p, o] for p in 1:Data.P) == 1 )
    @constraint(m, [o in Data.SO], sum(x[p, o] for p in 1:Data.P) <= 1 )
    @constraint(m, [p in 1:Data.P], sum(x[p, o] for o in 1:Data.O) <= Data.Capa[p])

    @objective(m, Min, sum(x[p, o] * alpha[p, i] * Data.Q[i][o] for o in 1:Data.O, p in 1:Data.P, i in 1:Data.N) - sum(x[p, o] for p in 1:Data.P, o in Data.SO))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        vX = JuMP.value.(x)
        #println("Modèle orders ")
        #println(m)
        return v_obj, vX
    else
        println("No solution for order subproblem")
    end
end

function slaveProblemRacks(alpha)
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, y[1:Data.P, 1:Data.R], Bin)
    @constraint(m, [r in 1:Data.R], sum(y[p, r] for p in 1:Data.P) <= 1)

    @objective(m, Min, sum(y[p, r] for p in 1:Data.P, r in 1:Data.R) * (length(Data.SO) + 1) - sum(alpha[p, i] * y[p, r] * Data.S[i][r] for p in 1:Data.P, r in 1:Data.R, i in 1:Data.N))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        vY = JuMP.value.(y)
        #println("Modèle racks ")
        #println(m)
        return v_obj, vY
    else
        println("No solution for rack subproblem")
    end
end

function masterProblem(orders, racks)
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, lo[1:length(orders)] >= 0)
    @variable(m, lr[1:length(racks)] >= 0)
    @variable(m, s >= 0)

    @constraint(m, etao, sum(lo[i] for i in 1:length(orders)) == 1)
    @constraint(m, etar, sum(lr[i] for i in 1:length(racks)) == 1)
    @constraint(m, alpha[p in 1:Data.P, i in 1:Data.N], sum(lo[k] * Data.Q[i][o] * orders[k][p, o] for k in 1:length(orders), o in 1:Data.O) <= s + sum(lr[k] * Data.S[i][r] * racks[k][p, r] for k in 1:length(racks), r in 1:Data.R))
    @objective(m, Min, 10000 * s + (length(Data.SO) + 1) * sum(lr[k] * racks[k][p, r] for k in 1:length(racks), p in 1:Data.P, r in 1:Data.R) - sum(lo[k] * orders[k][p, o] for k in 1:length(orders), p in 1:Data.P, o in Data.SO))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        Alpha = JuMP.dual.(alpha)
        Etar = JuMP.dual(etar)
        Etao = JuMP.dual(etao)
        return v_obj, - Alpha, Etao, Etar
    else
        println("No solution for master problem")
    end

end


function columnGeneration()
    Alpha = zeros(Data.P, Data.N)
    orders = [slaveProblemOrders(Alpha)[2]]
    racks = [slaveProblemRacks(Alpha)[2]]
    optimal = false
    v_obj = -1
    iter = 1
    plotDataMaster = []
    plotDataDual = []
    iters = []
    while !optimal
        v_obj, Alpha, Etao, Etar = masterProblem(orders, racks)
        push!(plotDataMaster, v_obj)
        println("Valeur courrante du problème maître : ", v_obj)
        vr, xr = slaveProblemRacks(Alpha)
        vo, xo = slaveProblemOrders(Alpha)
        push!(plotDataDual, vo + vr)
        push!(iters, iter)
        if vr - Etar < 1e-5
            push!(racks, xr)
        end
        if vo - Etao < 1e-5
            push!(orders, xo)
        end
        if vr - Etar >= -1e-5 && vo - Etao >= -1e-5
            optimal = true
        end
        iter += 1
    end
    return v_obj, iters, plotDataMaster, plotDataDual
end


P = 5
capa = Vector{Int}([12, 12, 12, 12, 12])
FO = 1:5
SO = 6:10

Data.P = P
Data.Capa = capa
Data.FO = FO
Data.SO = SO
sort!(Data.Capa, rev=true)


println("TEST")
#PLNE()
v_obj, iters, plotMaster, plotDual = columnGeneration()
plot(iters[7:end], plotMaster[7:end])
plot!(iters[7:end], plotDual[7:end])
