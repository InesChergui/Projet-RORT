using JuMP
using CPLEX

include("donnees.jl")


function PLNEOrdo(Data, maxTime = 600)
    m = Model(CPLEX.Optimizer)
    set_attribute(m, "CPXPARAM_TimeLimit", maxTime)
    #set_attribute(m, "CPXPARAM_ScreenOutput", 0)

    @variable(m, alpha >= 0)
    @variable(m, x[1:Data.O, 1:Data.R], Bin)
    @variable(m, z[1:Data.R, 1:Data.R], Bin)
    @variable(m, y[1:Data.O, 1:Data.N, 1:Data.R], Int)
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

    start_time = time()
    optimize!(m)
    end_time = time()

    feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    println(feasible_solution_found)
    if feasible_solution_found
        v_obj = JuMP.objective_value(m)
        println("Valeur fonction objective : ", v_obj)
        return v_obj, objective_bound(m), end_time - start_time
    end
end




function writeLine(instanceName, Data)
    v_plne, LB_plne, t_plne = PLNEOrdo(Data)
    file = open("resultTableOrdo.csv", "a")
    write(file, instanceName)
    write(file, ";")
    write(file, string((round(v_plne, digits=3))))
    write(file, ";")
    write(file, string((round(LB_plne, digits=3))))
    write(file, ";")
    write(file, string((round(t_plne, digits=3))))
    write(file, "\n")
    close(file)
end

function resultTable(instanceName)
    Data = readInstance("instances/" * instanceName * ".txt")
    writeLine(instanceName, Data)
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