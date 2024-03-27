using JuMP
using CPLEX
using Plots
using Statistics

include("donnees.jl")


function PLNE(Data)
    m = Model(CPLEX.Optimizer)

    @variable(m, x[1:Data.P, 1:Data.O], Bin)
    @variable(m, y[1:Data.P, 1:Data.R], Bin)

    @constraint(m, [r in 1:Data.R], sum( y[p, r] for p in 1:Data.P ) <= 1 )
    @constraint(m, [o in Data.FO], sum( x[p, o] for p in 1:Data.P) == 1 )
    @constraint(m, [o in Data.SO], sum( x[p, o] for p in 1:Data.P) <= 1 )
    @constraint(m, [p in 1:Data.P], sum(x[p, o] for o in 1:Data.O) <= Data.Capa[p])
    @constraint(m, [i in 1:Data.N, p in 1:Data.P], sum(Data.Q[i][o]*x[p, o] for o in 1:Data.O) <= sum(Data.S[i][r]*y[p, r] for r in 1:Data.R))

    @objective(m, Min, sum(y[p, r] for p in 1:Data.P, r in 1:Data.R)*(length(Data.SO) + 1) - sum(x[p, o] for p in 1:Data.P, o in Data.SO))

    start_time = time()
    optimize!(m)
    end_time = time()

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    println(feasible_solution_found)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        #println("Valeur fonction objective : ", v_obj)
        return v_obj, JuMP.objective_bound(m), end_time - start_time
    end
end


function slaveProblemOrders(Data, alpha)
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
        return v_obj, vX
    else
        println("No solution for order subproblem")
    end
end

function slaveProblemRacks(Data, alpha)
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_ScreenOutput" => 0))

    @variable(m, y[1:Data.P, 1:Data.R], Bin)
    @constraint(m, [r in 1:Data.R], sum(y[p, r] for p in 1:Data.P) <= 1)

    @objective(m, Min, sum(y[p, r] for p in 1:Data.P, r in 1:Data.R) * (length(Data.SO) + 1) - sum(alpha[p, i] * y[p, r] * Data.S[i][r] for p in 1:Data.P, r in 1:Data.R, i in 1:Data.N))

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

function masterProblem(Data, orders, racks)
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


function columnGeneration(Data)
    start_time = time()
    Alpha = zeros(Data.P, Data.N)
    orders = [slaveProblemOrders(Data, Alpha)[2]]
    racks = [slaveProblemRacks(Data, Alpha)[2]]
    optimal = false
    v_obj = -1
    iter = 1
    plotDataMaster = []
    plotDataDual = []
    iters = []
    while !optimal && iter < 1000
        v_obj, Alpha, Etao, Etar = masterProblem(Data, orders, racks)
        push!(plotDataMaster, v_obj)
        if iter % 200 == 0
            println("Valeur courrante du problème maître : ", v_obj)
        end
        vr, xr = slaveProblemRacks(Data, Alpha)
        vo, xo = slaveProblemOrders(Data, Alpha)
        if length(plotDataDual) == 0
            push!(plotDataDual, vo + vr)
        else
            push!(plotDataDual, max(plotDataDual[end], vo + vr))
        end
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
    end_time = time()
    return v_obj, iters, plotDataMaster, plotDataDual, end_time - start_time
end

function plotAndSave(instanceName, FO, iters, plotDataMaster, plotDataDual)
    plot(iters, plotDataMaster, label="Valeur du problème maître")
    plot!(iters, plotDataDual, label="Borne inférieure")
    masterMean = mean(plotDataMaster)
    dualMean = mean(plotDataDual)
    marginMaster = abs(masterMean - plotDataMaster[end])
    marginDual = abs(plotDataDual[end] - dualMean)
    ylims!((plotDataDual[end] - div(marginDual, 3), plotDataMaster[end] + div(marginMaster, 2)))
    savefig(instanceName * "_" * string(floor(FO)) * ".png")
end

function writeLine(instanceName, Data)
    v_plne, LB_plne, t_plne = PLNE(Data)
    v_cg, iters, plotDataMaster, plotDataDual, t_cg = columnGeneration(Data)
    file = open("resultTable.csv", "a")
    write(file, instanceName)
    write(file, ";")
    write(file, string(floor(Data.P)))
    write(file, ";")
    write(file, string(floor(Data.Capa[1])))
    write(file, ";")
    write(file, string(floor(Data.FO[end])))
    write(file, ";")
    write(file, string((round(v_plne, digits=3))))
    write(file, ";")
    write(file, string((round(LB_plne, digits=3))))
    write(file, ";")
    write(file, string((round(t_plne, digits=3))))
    write(file, ";")
    write(file, string((round(v_cg, digits=3))))
    write(file, ";")
    write(file, string(floor(iters[end])))
    write(file, ";")
    write(file, string((round(maximum(plotDataDual), digits=3))))
    write(file, ";")
    write(file, string((round(t_cg, digits=3))))
    write(file, "\n")
    close(file)
    plotAndSave(instanceName, length(Data.FO), iters, plotDataMaster, plotDataDual)
end

function resultTable(instanceName)
    Data = readInstance("instances/" * instanceName * ".txt")
    Data.Capa = Vector{Int}([12, 12, 12, 12, 12])
    if Data.N <= 20 #small instances
        Data.P = 2
        if Data.O % 2 == 0
            FO = div(Data.O, 2)
        else
            FO = div(Data.O, 2) + 1
        end
        Data.FO = 1:FO
        Data.SO = (Data.FO[end] + 1):Data.O
        println(instanceName)
        writeLine(instanceName, Data)
    else
        Data.P = 5
        for l in 1:5
            Data.FO = 1:(5 * l)
            Data.SO = (Data.FO[end] + 1):Data.O
            println(instanceName, " ", length(Data.FO))
            writeLine(instanceName, Data)
        end
    end
end

instanceNames = [#"Data_test_N5_R2_O3_RS2",
                 #"Data_test_N5_R3_O3_RS5",
                 #"Data_test_N5_R4_O3_RS2",
                 #"Data_test_N7_R4_O6_RS7",
                 #"Data_test_N7_R5_O5_RS7",
                 #"Data_test_N7_R5_O6_RS7",
                 #"Data_test_N10_R10_O10_RS7",
                 #"Data_test_N12_R12_O12_RS8",
                 #"instance_N100_R50_O50_RS25",
                 #"instance_N100_R100_O100_RS25",
                 #"instance_N100_R100_O150_RS25",
                 #"instance_N200_R50_O50_RS25",
                 #"instance_N200_R100_O100_RS25",
                 #"instance_N200_R100_O150_RS25",
                 #"instance_N300_R50_O50_RS25",
                 #"instance_N300_R100_O100_RS25",
                 #"instance_N300_R100_O150_RS25"
                 ]


for instance in instanceNames
    resultTable(instance)
end



# Data = readInstance("instances/Data_test_N10_R10_O10_RS7.txt")
# Data.P = 2
# Data.Capa = Vector{Int}([12, 12, 12, 12, 12])
# Data.FO = 1:div(Data.O, 2)
# Data.SO = (Data.FO[end] + 1):Data.O

# println(Data.FO)

# println("TEST")
# #println(PLNE(Data))
# v_obj, iters, plotMaster, plotDual = columnGeneration(Data)
# plot(iters, plotMaster)
# plot!(iters, plotDual)
# masterMean = mean(plotMaster)
# dualMean = mean(plotDual)
# marginMaster = abs(masterMean - plotMaster[end])
# marginDual = abs(plotDual[end] - dualMean)
# ylims!((plotDual[end] - div(marginDual, 2), plotMaster[end] + div(marginMaster, 2)))
# #savefig("output.png")
