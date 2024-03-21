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

  
  
function readInstance(path::String)
    lines = readlines(path) 
    line=lines[1]
    line_decompose=split(line)
    N=parse(Int64, line_decompose[2])
    #n=parse(Float64, line_decompose[1])
    #println("nbre produits total N ",N)

    line=lines[2]
    line_decompose=split(line)
    R=parse(Int64, line_decompose[2])
    #println("nbre racks R ",R)

    line=lines[3]
    line_decompose=split(line)
    O=parse(Int64, line_decompose[2])
    #println("nbre ordres O ",O)

    line=lines[5]
    line_decompose=split(line)
    RS=parse(Int64, line_decompose[2])
    #println("nbre shelves par rack ",RS)

    Data=donnees(N,R,O,RS)
    #println(Data.N," ",Data.R," ",Data.O," ",Data.RS)
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
    return Data
end