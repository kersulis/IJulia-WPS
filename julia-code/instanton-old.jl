# Julia packages used in this code
using JuMP   # numerical optimization package
using Mosek  # Julia interface to MOSEK solver
using MAT    # Julia interface to MATLAB .mat binary files
using PyCall # used to call python commands
using PyPlot # Julia plotting package; wraps around Matplotlib
PyPlot.svg(true)

function createY(f,t,r,x,b,s)
    # Create an admittance matrix from three vectors: from, to, and adm. value
    # Note: 'f' and 't' must be integer vectors.
    # DC if s == true
    
    if s == true
        y = 1./x
        return sparse([f; t; t; f],[t; f; t; f],[-y; -y; y; y])
    else
        G = 1./r
        G[G.==Inf] = 0
        B = 1./x
        y = complex(G,B)
        return sparse([f; t; t; f],[t; f; t; f],[-y; -y; y + b./2; y + b./2])
    end
end

function zero2NaN(x) # replace 0 with NaN for plotting purposes
    y = deepcopy(x)
    y[y.==0] = NaN
    return y
end

function readRTS96Data();
    # Load RTS96 data
    # uses MAT
    path="/home/jkersulis/Documents/MATPOWERcases/caseRTS96.mat"
    caseRTS96 = matread(path); # import MATLAB workspace

    # Connect relevant MATLAB variables to Julia
    bus_i = caseRTS96["bus_i"][:,1]; # vector of unique bus indices (73)
    bus = caseRTS96["bus"][:,1]; # vector of generator bus indices (99)
    Sb = caseRTS96["Sb"]; # base MVA
    bustype = caseRTS96["type"][:,1]

    Gp_long = caseRTS96["Pg"]; # conventional active power output, divide by Sb later
    Gq_long = caseRTS96["Qg"]; # conventional reactive power output, divide by Sb later

    Rp = caseRTS96["Wind"]./Sb; Rp = Rp[:,1]; # renewable active generation, has zero where no farm exists
    Rq = caseRTS96["Wind"].*(0.484/Sb); Rq = Rq[:,1];

    Dp = caseRTS96["Pd"]./Sb + Rp; # I previously subtracted wind from active demand
    Dp = Dp[:,1];
    Dq = caseRTS96["Qd"]./Sb; # reactive power demand vector
    Dq = Dq[:,1];

    Pmax = caseRTS96["Pmax"]./Sb; Pmax = Pmax[:,1]; # conventional gen. max active output
    Pmin = caseRTS96["Pmin"]./Sb; Pmin = Pmin[:,1]; # conventional gen. min active output

    Qmax = caseRTS96["Qmax"]./Sb; Qmax = Qmax[:,1]; # conventional gen. max reactive output
    Qmin = caseRTS96["Qmin"]./Sb; Qmin = Qmin[:,1]; # conventional gen. min reactive output

    Vm = caseRTS96["Vm"][:,1]; # should be vector of ones
    Va = caseRTS96["Va"][:,1]; # should be vector of zeros

    Plim = caseRTS96["rateA"]./Sb; Plim = Plim[:,1];
    Vceiling = caseRTS96["Vmax"][:,1]; # one voltage mag. ceiling constraint per node
    Vfloor = caseRTS96["Vmin"][:,1]; # one voltage mag. floor constraint per node

    Vceiling[bustype .== 2] = 1; # ceiling is 1 for voltage-controlled nodes
    Vfloor[bustype .== 2] = 1; # floor is 1 for voltage-controlled nodes

    Vg = caseRTS96["Vg"][:,1]; # voltage magnitudes from Jenny's solved ACOPF

    f = int(caseRTS96["fbus"][:,1]); # "from bus" ...
    t = int(caseRTS96["tbus"][:,1]); # ... "to bus"
    r = caseRTS96["r"][:,1]; # resistance, pu
    x = caseRTS96["x"][:,1]; # reactance, pu
    b = caseRTS96["b"][:,1];  # susceptance, pu

    # map bus numbers to 1:73
    for i = 1:length(f)
        if f[i] < 201
            f[i] -= 100
        elseif f[i] < 301
            f[i] -= 176;
        else
            f[i] -= 252;
        end
    end
    for i = 1:length(t)
        if t[i] < 201
            t[i] -= 100;
        elseif t[i] < 301
            t[i] -= 176;
        else
            t[i] -= 252;
        end
    end

    busidx = Int64[1];
    for i = 2:length(bus)
        if bus[i] < 201
            push!(busidx, bus[i]-100);
        elseif bus[i] < 301
            push!(busidx, bus[i]-176);
        else
            push!(busidx, bus[i]-252);
        end
    end

    # area 1: buses 1-24
    # area 2: buses 25-48
    # area 3: buses 49-73

    Gp = zeros(length(bus_i))
    Gq = zeros(length(bus_i))

    for i in unique(busidx)
        Gp[i] = sum(Gp_long.*[busidx .== i])/Sb;
        Gq[i] = sum(Gq_long.*[busidx .== i])/Sb;
    end
    # Now Gp and Gq reflect active and reactive generation at buses 1:73 consecutively.

    N = length(bus_i);
    Nr = length(find(Rp));
    Ng = length(find(Gp));

    Y = createY(f,t,r,x,b,true)
    # spy(Y)
    # title(L"Spy plot of $Y_{bus}$")

    return  Sb, f, t, r, x, b, Y, bustype, 
            Gp, Gq, Dp, Dq, Rp, Rq, 
            Pmax, Pmin, Qmax, Qmin, 
            Plim, Vg, Vceiling, Vfloor, 
            bus_i, N, Nr, Ng
end

function DistDCModel(idx, sense, Rp, Gp, Dp, f, t, x, Y, bustype, Plim)
    # DISTRIBUTED SLACK
    # Create model saturating line 'idx' in direction 'sense' (±1)
    # This function uses JuMP and Mosek
    
    m = Model(solver = MosekSolver()) # Use MOSEK to solve model
    @defVar(m, ρ[1:N] >= 0) # Add decision variables to model (renewable gen)
    @defVar(m, θ[1:N]) # Add bus angles
    @defVar(m, α) # conv. gen. participation factor
    setObjective(m, :Min, 0.5*sum([(ρ[i] - Rp[i])^2 for i in find(Rp)]))

    # add power balance constraints:
    @defConstrRef PowerBalance[1:N]
    for i in setdiff(1:N,find(bustype.==3))
        PowerBalance[i] = @addConstraint(m, sum([Y[i,j]*θ[j] for j in 1:N]) == Gp[i]*α + ρ[i] - Dp[i])
    end
    @addConstraint(m, NonWind, sum([ρ[i] for i in setdiff(collect(1:N),find(Rp))]) == 0) # ρ=0 for non-wind nodes
    @addConstraint(m, Slack, θ[find(bustype.==3)[1]] == 0) # θ=0 for slack bus
    @addConstraint(m, Congestion, θ[f[idx]] - θ[t[idx]] == sense*x[idx]*Plim[idx]) # saturate a line
    @addConstraint(m, Participation, α == (sum(Dp) - sum([ρ[i] for i in find(Rp)]))/sum(Gp))
    return m, ρ, θ, α, PowerBalance, Slack, Congestion, Participation
end

function writeDot(name, busIdx, busInj, renGen, f, t, lineFlow, lineLim, size=(11,14))
    # This function generates a graph that richly expresses the RTS96 system state.
    
    busInj = round(busInj,2)
    lineFlow = round(lineFlow,2)
    
    # Open the dot file, overwriting anything already there:
    dotfile = open("$(name).dot","w")
    
    # Begin writing the dot file:
    write(dotfile, "digraph $(name) {\nnewrank=true;\n")

    # Set graph properties:
    write(dotfile, 
    "graph [fontname=helvetica, tooltip=\" \", overlap=false, size=\"$(size[1]),$(size[2])\", ratio=fill, orientation=\"portrait\",layout=dot];\n")

    # Set default node properties:
    write(dotfile, "node [fontname=helvetica, shape=square, style=filled, fontsize=20, color=\"#bdbdbd\"];\n")

    # Set default edge properties:
    write(dotfile, "edge [fontname=helvetica, style=\"setlinewidth(5)\"];\n")

    # Write bus data to dot file:
    write(dotfile, 
    "subgraph cluster_a1 {\nlabel=\"Area 1: High Wind\";\nfontcolor=\"#5677fc\";\nfontsize=24;\ncolor=\"#ffffff\";\nlabeljust=\"c\";\n")

    for i = 1:24
        write(dotfile, 
        "$(i) [label=$(int(busIdx[i])), tooltip=\"Inj = $(busInj[i])\"") # bus label and tooltip

        # Represent renewable nodes with blue circles:
        if union(find(renGen),i) == find(renGen)
            write(dotfile, ", shape=circle, color=\"#5677fc\"")
        end

        write(dotfile, "];\n")
    end
    write(dotfile, "}\n")

    write(dotfile, 
    "subgraph cluster_a2 {\nlabel=\"Area 2: Moderate Wind\";\nfontcolor=\"#5677fc\";\nfontsize=24;\ncolor=\"#ffffff\";\nlabeljust=\"l\";\n")

    for i = 25:48
        write(dotfile, 
        "$(i) [label=$(int(busIdx[i])), tooltip=\"Inj = $(busInj[i])\"") # bus label and tooltip

        # Represent renewable nodes with blue circles:
        if union(find(renGen),i) == find(renGen)
            write(dotfile, ", shape=circle, color=\"#5677fc\"")
        end

        write(dotfile, "];\n")
    end
    write(dotfile, "}\n")

    write(dotfile, 
    "subgraph cluster_a3 {\nlabel=\"Area 3: Low Wind\";\nfontcolor=\"#5677fc\";\nfontsize=24;\ncolor=\"#ffffff\";\nlabeljust=\"r\";\n")
    for i = 49:length(busIdx)
        write(dotfile, 
        "$(i) [label=$(int(busIdx[i])), tooltip=\"Inj = $(busInj[i])\"") # bus label and tooltip

        # Represent renewable nodes with blue circles:
        if union(find(renGen),i) == find(renGen)
            write(dotfile, ", shape=circle, color=\"#5677fc\"")
        end

        write(dotfile, "];\n")
    end
    write(dotfile, "}\n")

    # Write line data to file:

    for i = 1:length(f)

        normStrain = abs(lineFlow[i])/lineLim[i] # normalized strain on line i

        # Use flow direction to determine arrow direction,
        # label with flow,
        # and color according to strain
        if lineFlow[i] > 0
            write(dotfile, 
            "$(f[i]) -> $(t[i]) [label=$(lineFlow[i])")
        else
            write(dotfile, 
            "$(t[i]) -> $(f[i]) [label=$(-lineFlow[i])")
        end
        write(dotfile,
        ", tooltip=\" \", labeltooltip=\"Flow = $(int(normStrain*100))%\", color=\"$(round((1 - normStrain)/3,3)) 1.000 0.700\"];\n")
    end

    # Cap off the dot file and close it:
    write(dotfile, "}")
    close(dotfile)

    println("$(name).dot generated.\nUse \";dot -Tsvg dotfunctest.dot -o dotfunctest.svg\" to create an SVG.")
end

# Define a type to hold power system data:
type powerSystemData
    Sb
    f
    t
    r
    x
    b
    Y
    bustype
    Gp
    Gq
    Dp
    Dq
    Rp
    Rq
    Pmax
    Pmin
    Qmax
    Qmin
    Plim
    Vg
    Vceiling
    Vfloor
    busIdx
    N
    Nr
    Ng
end

# Define a type to hold instanton results. If there are N constraints, there will be 2N elements in 'score', 2N vectors in 'ρ',
# 2N vectors in 'θ', and 2N elements in 'α'. There will also be 2N values in 'constrIdx', each of which should be a string of the
# form:
#    "activeFlow: i -> k @ x pu"
#    "voltageMagnitude: i @ x pu"
#    "apparentFlow: i -> k @ x pu"

# I may need to modify this class to contain voltage magnitudes later on.
type instantonResults
    score
    ρ
    θ
    α
    constrIdx
end

# Function used to create instance of powerSystemData type from MATLAB workspace:
function psDataLoad()
    (Sb, f, t, r, x, b, Y, bustype, 
     Gp, Gq, Dp, Dq, Rp, Rq, 
     Pmax, Pmin, Qmax, Qmin, 
     Plim, Vg, Vceiling, Vfloor, 
    busIdx, N, Nr, Ng) = readRTS96Data()
    return powerSystemData(Sb,f,t,r,x,b,Y,bustype,Gp,Gq,Dp,Dq,Rp,Rq,Pmax,Pmin,Qmax,Qmin,Plim,Vg,Vceiling,Vfloor,busIdx,N,Nr,Ng)
end

function unpack_psDL(psDL)
# code used to unpack psDL type instance:
return (psDL.Sb,
    psDL.f,
    psDL.t,
    psDL.r,
    psDL.x,
    psDL.b,
    psDL.Y,
    psDL.bustype,
    psDL.Gp,
    psDL.Gq,
    psDL.Dp,
    psDL.Dq,
    psDL.Rp,
    psDL.Rq,
    psDL.Pmax,
    psDL.Pmin,
    psDL.Qmax,
    psDL.Qmin,
    psDL.Plim,
    psDL.Vg,
    psDL.Vceiling,
    psDL.Vfloor,
    psDL.busIdx,
    psDL.N,
    psDL.Nr,
    psDL.Ng)
end

function instantonAnalysis(psData; method="lagrange",constraintType="activeFlow",genResponse="droop")
    # This function accepts the following arguments:
    # psData, a dict containing all system data
    # method, a property that can be set to either "lagrange" or "solver"
    # constraintType, a property that can be set to "activeFlow", "apparentFlow", or "voltageMagnitude"
    # genResponse, a property that can be set to "none", "proportional", or "droop"
    # The function returns the following arguments:
    # rankedList, a dict containing score, renewable generation, phase angle, and response parameter data.

    # Copy data (will be unpacked by functions later):
    psDL = deepcopy(psData)

    # The code takes on an if/elseif/else structure to handle various combinations of constraint types and generator responses.
    # The matrix is far from filled in, but this exhaustive enumeration is probably the best way to structure the code.
    if method=="lagrange"
        if constraintType=="activeFlow"
            if genResponse=="droop"
                # DONE: perform Lagrangian analysis with active power flows and droop response
                return lagrange_activeFlow_droop(psDL)
                
                elseif genResponse=="proportional"
                # function for performing Lagrangian analysis with active power flows and proportional response
                
                elseif genResponse=="none"
                # function for performing Lagrangian analysis with active power flows and no response
                
                else println("genResponse value invalid.")
            end
            
            elseif constraintType=="voltageMagnitude"
            # function for performing Lagrangian analysis with voltage magnitudes
            # Because any model of generator reactive power response is beyond the scope of this project,
            # assume there is no conventional generator response.
            
            elseif constraintType=="apparentFlow"
            if genResponse=="droop"
                # DONE: perform Lagrangian analysis with apparent power flows and droop response
                
                elseif genResponse=="proportional"
                # function for performing Lagrangian analysis with apparent power flows and proportional response
                
                elseif genResponse=="none"
                # function for performing Lagrangian analysis with apparent power flows and no response
                
                else println("genResponse value invalid.")
            end
            
            else println("constraintType value invalid.")
        end
        
        elseif method=="solver"
        if constraintType=="activeFlow"
            if genResponse=="droop"
                # DONE: perform solver analysis with active power flows and droop response
                
                elseif genResponse=="proportional"
                # DONE: perform solver analysis with active power flows and proportional response
                
                elseif genResponse=="none"
                # perform solver analysis with active power flows and no response
                
                else println("genResponse value invalid.")
            end
            
            elseif constraintType=="voltageMagnitude"
            # function for performing solver analysis with voltage magnitudes
            # Because any model of generator reactive power response is beyond the scope of this project,
            # assume there is no conventional generator response.
            
            elseif constraintType=="apparentFlow"
            if genResponse=="droop"
                # perform solver analysis with apparent power flows and droop response
                
                elseif genResponse=="proportional"
                # perform solver analysis with apparent power flows and proportional response
                
                elseif genResponse=="none"
                # perform Lagrangian analysis with apparent power flows and no response
                
                else println("genResponse value invalid.")
            end
            
            else println("constraintType value invalid.")
        end
        else println("method value invalid.")
    end
end

function lagrange_activeFlow_droop(psDL)
    # This block of code handles the case with "lagrange", "activeFlow", and "droop" settings.
    
    # unpack psDL (boilerplate):
    (Sb,f,t,r,x,b,Y,bustype,
    Gp,Gq,Dp,Dq,Rp,Rq,
    Pmax,Pmin,Qmax,Qmin,Plim,
    Vg,Vceiling,Vfloor,
    busIdx,N,Nr,Ng) = unpack_psDL(psData)
    
    # Prepare all variables:
    sorted = [find(Rp); setdiff(1:N,find(Rp))]
    ρ0 = Rp[find(Rp)]
    Λ = speye(Nr)
    Nc = N - Nr

    Yr = Y[sorted[1:Nr],:]
    Yc = Y[sorted[Nr+1:end],:]

    G0 = sparse(Gp)
    k = spzeros(N,1)
    for i = 1:N
        if G0[i] != 0
            k[i] = 1/length(find(Gp))
        end
    end
    k = k[sorted]
    kr = k[1:Nr]
    kc = k[Nr+1:end]

    G0r = G0[sorted[1:Nr]]
    G0c = G0[sorted[Nr+1:end]]

    D = Dp
    Dr = D[sorted[1:Nr]]
    Dc = D[sorted[Nr+1:end]]

    # sref is the column vector establishing angle reference:
    sref = spzeros(N,1)
    sref[find(bustype.==3)] = 1

    ρ = Array(Vector{Float64},0) # array of vectors with Float64 values
    θ = Array(Vector{Float64},0) # array of vectors with Float64 values
    α = Float64[]
    score = Float64[]

    instInvert(A,B,i) = try
        \(A,B)
        catch y
        println("$(y), i = $(i)")
        return zeros(length(A)).*NaN
    end

    # loop through all lines. SingularException: 52, 90
    for i = 1:2length(f)

        sik = spzeros(N,1)
        if i <= length(f)
            sik[f[i]] =  Y[f[i],t[i]]
            sik[t[i]] = -Y[f[i],t[i]]
            Plimik = Plim[i]
        else
            i2 = i - length(f)
            sik[f[i2]] = -Y[f[i2],t[i2]]
            sik[t[i2]] =  Y[f[i2],t[i2]]
            Plimik = Plim[i2]
        end

    #                   Nr cols        N cols     1 col     Nr cols     Nc cols          1 col             1 col      1 col
            A = sparse([   Λ         zeros(Nr,N+1)         -eye(Nr)   zeros(Nr,Nc)  -ones(Nr,1)                zeros(Nr,2); # Nr rows
                               zeros(N,Nr+N+1)               Yr'          Yc'        zeros(N,1)             sref       sik; # N rows
                               zeros(1,Nr+N+1)              -kr'         -kc'            -1                    zeros(1,2) ; # 1 row
                       -eye(Nr)        Yr        -kr                          zeros(Nr,N+3)                               ; # Nr rows
                     zeros(Nc,Nr)      Yc        -kc                          zeros(Nc,N+3)                               ; # Nc rows
                      -ones(1,Nr)   zeros(1,N)   -1                           zeros(1,N+3)                                ; # 1 row
                     zeros(1,Nr)     sref'                                 zeros(1,N+4)                                   ; # 1 row
                     zeros(1,Nr)     sik'                                  zeros(1,N+4)                                  ]) # 1 row

            B = Float64[]
            append!(B,Λ*ρ0)
            append!(B,zeros(N+1))
            append!(B,[G0r-Dr][:,1])
            append!(B,[G0c-Dc][:,1])
            push!(B,sum(G0) - sum(D))
            push!(B,0)
            push!(B,Plimik)

            X = instInvert(A,B,i)
            push!(ρ, X[1:Nr])
            push!(θ, X[Nr+1:Nr+N]) # angles sorted by original bus numbers
            push!(α, X[Nr+N+1])

        if sum(abs(X[1:Nr])) == 0
            push!(score, NaN)
        else
            push!(score, [0.5(X[1:Nr] - ρ0)'*(X[1:Nr] - ρ0)][1])
        end
    end
    
    constrIdx = String[]
    # Generate strings to tell user which constraint was violated:
    for i = 1:2length(f)
        if i <= length(f)
            push!(constrIdx, "activeFlow: $(f[i]) -> $(t[i]) @ $(Plim[i])")
        else
            i2 = i - length(f)
            push!(constrIdx, "activeFlow: $(t[i2]) -> $(f[i2]) @ $(Plim[i2])")
        end
    end
    # Obtain ranked list of events:
    rankedList = sortrows([zero2NaN(score) [ρ θ α constrIdx]]) # sort candidates; place NaNs at bottom
    
    # Pack results into instance of "instantonResults":
    return instantonResults(rankedList[:,1],rankedList[:,2],rankedList[:,3],rankedList[:,4],rankedList[:,5])
end
