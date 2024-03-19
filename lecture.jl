# lecture donnees
mutable struct donnees
  N::Int # nbre de produits
  R::Int # nbre racks
  O::Int # nbre d'ordre
  RS::Int # nbre shelves dans un rack
  S::Vector{Vector{Int}} # matrice produits - racks
  Q::Vector{Vector{Int}} # matrice produits - ordres
  P::Int # nbre de pickers
  Capa::Vector{Int} # capacite des pickers
  FO::Vector{Int} # first orders prioritaires
  SO::Vector{Int} # second orders moins prioritaires
  function donnees(N,R,O,RS)
    this=new()
    this.N=N
    this.R=R
    this.O=O
    this.RS=RS
    this.P=0
    # init matrice S
    this.S=[] # vide
    for i in 1:N 
      push!(this.S,[])
    end
    for i in 1:N 
      for j in 1:R
        push!(this.S[i],0)
      end
    end
    # init matrice Q
    this.Q=[] # vide
    for i in 1:N 
      push!(this.Q,[])
    end
    for i in 1:N 
      for j in 1:O
        push!(this.Q[i],0)
      end
    end
    # init matrice Capacite des pickers
    this.Capa=[] # vide
    # init matrice FO vide
    this.FO=[]
    # init matrice SO vide
    this.SO=[]
    
    return this
  end
end # fin de la struct donnees

  println("je suis dans la fonction lecture des donnees 0")

  #lines = readlines("test_data.txt") # le fichier contenant les donnees
  #lines = readlines("test_data_1.txt")
  # lines = readlines("instances/instance_N100_R50_O50_RS25.txt")
  #lines = readlines("instances/instance_N100_R100_O100_RS25.txt")
  #lines = readlines("instances/instance_N100_R100_O150_RS25.txt")
  #lines = readlines("instances/instance_N200_R50_O50_RS25.txt")
  #lines = readlines("instances/instance_N200_R100_O100_RS25.txt")
  #lines = readlines("instances/instance_N300_R50_O50_RS25.txt")

  #lines = readlines("Data_test.txt")
  #lines = readlines("instances/Data_test_N5_R4_O3_RS2.txt")
  #lines = readlines("instances/Data_test_N5_R3_O3_RS5.txt")
  #lines = readlines("instances/Data_test_N5_R2_O3_RS2.txt")
  #lines = readlines("instances/Data_test_N7_R5_O8_RS5.txt") # ici 5/8 des ordres peuvent etre satisfaits
  #lines = readlines("instances/Data_test_N7_R5_O6_RS7.txt")
  lines = readlines("instances/Data_test_N10_R10_O10_RS7.txt") 
  #lines = readlines("instances/Data_test_N12_R12_O12_RS8.txt")
  #lines = readlines("instances/Data_test_N14_R14_O14_RS10.txt")  



println("je suis dans la fonction lecture des donnees 1")

  line=lines[1]
  line_decompose=split(line)
  N=parse(Int64, line_decompose[2])
  #n=parse(Float64, line_decompose[1])
println("nbre produits total N ",N)

  line=lines[2]
  line_decompose=split(line)
  R=parse(Int64, line_decompose[2])
println("nbre racks R ",R)

  line=lines[3]
  line_decompose=split(line)
  O=parse(Int64, line_decompose[2])
println("nbre ordres O ",O)

  line=lines[5]
  line_decompose=split(line)
  RS=parse(Int64, line_decompose[2])
println("nbre shelves par rack ",RS)

Data=donnees(N,R,O,RS)
println(Data.N," ",Data.R," ",Data.O," ",Data.RS)
for r in 1:R # parcours les racks
  num_line=7+r
  global line=lines[num_line]
  global line_decompose=split(line)
  # print("\n rack ",r,"\n")
  for i in 1:RS # parcours les shelves 
    num_prod=parse(Int64,line_decompose[2*i])
    quantite=parse(Int64,line_decompose[1+2*i])

    # print(num_prod," ",quantite," ")
    # Attention les produits vont de 0 - N-1
    Data.S[num_prod+1][r]=quantite
  end
  
end
for o in 1:O # parcours les ordres
  num_line=(7+R+2)+o
  global line=lines[num_line]
  global line_decompose=split(line)
  
  nbre_prod_inside_ordre=parse(Int64,line_decompose[2])
  # print("\n ordre ",o," ",nbre_prod_inside_ordre,"\n")
  for i in 1:nbre_prod_inside_ordre
    num_prod=parse(Int64,line_decompose[2+i])
    # Attention numero de produit vont de 0 - N-1
    Data.Q[num_prod+1][o]+=1
    # println("num produit ", num_prod," ",Data.Q[num_prod+1][o])
  end
end





 