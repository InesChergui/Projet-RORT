using JuMP
using CPLEX

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
    @constraint(m, [p in 1:Data.P], sum(x[p, o] for o in 1:O) <= Data.Capa[p])

    @objective(m, Min, sum(x[p, o] * alpha[p, i] * Data.Q[i][o] for o in 1:Data.O, p in 1:Data.P, i in 1:Data.N) - sum(x[p, o] for p in 1:Data.P, o in Data.SO))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        vX = JuMP.value.(x)
        return v_obj, vX
    else
        println("No solution for order subproblem")
    end
end

function slaveProblemRacks(alpha)
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, y[1:Data.P, 1:Data.R], Bin)
    @constraint(m, [r in 1:Data.R], sum(y[p, r] for p in 1:Data.P) <= 1)

    @objective(m, Min, sum(y[p, r] for p in 1:Data.P, r in 1:Data.R) * (length(Data.SO) + 1) - sum(alpha[p, i] * y[p, r] * Data.S[i][r] for p in 1:Data.P, r in 1:Data.R, i in 1:N))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        vY = JuMP.value.(y)
        return v_obj, vY
    else
        println("No solution for rack subproblem")
    end
end

function masterProblem(orders, racks)
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, lo[1:length(orders)] >= 0)
    @variable(m, lr[1:length(racks)] >= 0)

    @constraint(m, etao, sum(lo[i] for i in 1:length(orders)) == 1)
    @constraint(m, etar, sum(lr[i] for i in 1:length(racks)) == 1)
    @constraint(m, alpha[p in 1:Data.P, i in 1:Data.N], sum(lo[k] * Data.Q[i][o] * orders[k][p, o] for k in 1:length(orders), o in 1:Data.O) <= sum(lr[k] * Data.S[i][r] * racks[k][p, r] for k in 1:length(racks), r in 1:Data.R))
    @objective(m, Min, (length(Data.SO) + 1) * sum(lr[k] * racks[k][p, r] for k in 1:length(racks), p in 1:Data.P, r in 1:Data.R) - sum(lo[k] * orders[k][p, o] for k in 1:length(orders), p in 1:Data.P, o in Data.SO))

    optimize!(m)

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        Alpha = JuMP.dual.(alpha)
        Etar = JuMP.dual(etar)
        Etao = JuMP.dual(etao)
        return v_obj, Alpha, Etao, Etar
    else
        println("No solution for master problem")
    end

end


function dualCuttingPlanes()
    Alpha = zeros(Data.P, Data.N)
    xo = initialOrders()
    xr = initialRacks(xo)

    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, etao)
    @variable(m, etar)
    @variable(m, alpha[1:Data.P, 1:Data.N] >= 0)

    @constraint(m, etao - sum(alpha[p, i] * Data.Q[i][o] * xo[p, o] for p in 1:Data.P, i in 1:Data.N, o in 1:Data.O) + sum(xo[p, o] for p in 1:Data.P, o in Data.SO) <= 0)
    @constraint(m, etar + sum(alpha[p, i] * Data.S[i][r] * xr[p, r] for p in 1:Data.P, i in 1:Data.N, r in 1:Data.R) - (length(Data.SO) + 1) * sum(xr[p, r] for p in 1:Data.P, r in 1:Data.R) <= 0)
    @objective(m, Max, etao + etar)

    optimal = false

    while !optimal
        optimize!(m)
        v_obj = JuMP.objective_value(m)
        println("Valeur de l'objectif courrant : ", v_obj)
        Alpha = JuMP.value.(alpha)
        Etao = JuMP.value(etao)
        Etar = JuMP.value(etar)

        vr, xr = slaveProblemRacks(Alpha)
        vo, xo = slaveProblemOrders(Alpha)

        if vr - Etar < 0
            @constraint(m, etar + sum(alpha[p, i] * Data.S[i][r] * xr[p, r] for p in 1:Data.P, i in 1:Data.N, r in 1:Data.R) - (length(Data.SO) + 1) * sum(xr[p, r] for p in 1:Data.P, r in 1:Data.R) <= 0)
        end
        if vo - Etao < 0
            @constraint(m, etao - sum(alpha[p, i] * Data.Q[i][o] * xo[p, o] for p in 1:Data.P, i in 1:Data.N, o in 1:Data.O) + sum(xo[p, o] for p in 1:Data.P, o in Data.SO) <= 0)
        end
        if vr - Etar >= 0 && vo - Etao >= 0
            optimal = true
        end
    end
end

function columnGeneration()
    Alpha = zeros(Data.P, Data.N)
    orders = [initialOrders()]
    racks = [initialRacks(orders[1])]
    optimal = false
    v_obj = -1
    iter = 0
    while !optimal
        v_obj, Alpha, Etao, Etar = masterProblem(orders, racks)
        println("Valeur courrante du problème maître : ", v_obj)
        vr, xr = slaveProblemRacks(Alpha)
        vo, xo = slaveProblemOrders(Alpha)
        if vr - Etar < 0
            push!(racks, xr)
        end
        if vo - Etao < 0
            push!(orders, xo)
        end
        if vr - Etar >= 0 && vo - Etao >= 0
            optimal = true
        end
        iter += 1
    end
    return v_obj
end


P = 5
capa = Vector{Int}([12, 12, 12, 12, 12])
FO = 1:5
SO = 6:50

Data.P = P
Data.Capa = capa
Data.FO = FO
Data.SO = SO
sort!(Data.Capa, rev=true)


PLNE()
columnGeneration()
