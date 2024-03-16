using JuMP
using CPLEX
include("lecture.jl")

sortedRacks = []
for i in 1:N
    push!(sortedRacks, sort(1:Data.R, rev=true, by = r -> Data.S[i][r]))
end

function initialOrders()
    xo = zeros(Data.P, Data.O)
    p = 1
    o = 1
    capa = 0
    while o <= length(Data.FO)
        if capa > Data.Capa[p]
            p += 1
            capa = 0
        end
        xo[p, Data.FO[o]] = 1
        capa += 1
        o += 1
    end
    return xo
end

function initialRacks(xo)
    xr = zeros(Data.P, Data.R)
    assignedRacks = []
    for p in 1:Data.P
        sortedProductsByQuantity = sort(1:Data.N, rev=true, by = i -> sum(Data.Q[i][o] * xo[p, o] for o in 1:Data.O))
        for product in 1:length(sortedProductsByQuantity)
            while sum(Data.Q[sortedProductsByQuantity[product]][o]*xo[p, o] for o in 1:O) > sum(Data.S[sortedProductsByQuantity[product]][r]*xr[p, r] for r in 1:R)
                rack = 1
                while sortedRacks[sortedProductsByQuantity[product]][rack] in assignedRacks
                    rack += 1
                    if rack == length(sortedRacks[sortedProductsByQuantity[product]]) + 1
                        println("Infaisable!")
                        break
                    end
                end
                xr[p, sortedRacks[sortedProductsByQuantity[product]][rack]] = 1
                push!(assignedRacks, sortedRacks[sortedProductsByQuantity[product]][rack])
            end
        end
    end
    for p in 1:Data.P
        for i in 1:Data.N
            @assert sum(Data.Q[i][o]*xo[p, o] for o in 1:O) <= sum(Data.S[i][r]*xr[p, r] for r in 1:R)
        end
    end
    return xr
end


