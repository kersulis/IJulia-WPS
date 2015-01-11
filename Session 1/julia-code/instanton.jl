# This julia file loads all packages used in the remainder of the instanton analysis code.
# It also defines two custom types used to store power system data and instanton results.

# Julia packages used in this code
using JuMP   # numerical optimization package
using Ipopt  # Julia interface to Ipopt solver
using MAT    # Julia interface to MATLAB .mat binary files
using PyCall # used to call python commands
using PyPlot # Julia plotting package; wraps around Matplotlib
PyPlot.svg(true)

# Other routines needed to perform instanton analysis:
include("instantonAnalysis.jl")                 # Main function
include("plotInstantonResults.jl")              # Plotting function

include("psDataLoad.jl")                        # Loads power system data from MATLAB data
include("readRTS96Data.jl")                     # Called by psDataLoad()
include("writeDot.jl")                          # Produces a dot file showing power flow on the RTS-96 network

#include("lagrange_activeFlow_droop.jl")         # Lagrange with active flows and droop response
#include("solver_activeFlow_proportional.jl")    # Solver with active flows and proportional response
include("solver_activeFlow_droop.jl")           # Solver with active flows and droop response
#include("solver_apparentFlow_droop.jl")         # Solver with active flows and droop response
#include("solver_currentFlow_droop.jl")          # Solver with current magnitude limits and droop response

#include("solver_activeFlow_droop.jl")  # Solver with active flows and droop response

# Define a type to hold power system data:
type powerSystemData
    Sb          # Base complex power
    f           # Lines: "from"
    t           # Lines: "to"
    r           # Lines: resistance
    x           # Lines: reactance
    b           # Lines: susceptance
    Y           # Lines: admittance
    bustype     # Buses: type
    Gp          # Buses: conv. active gen
    Gq          # Buses: conv. reactive gen
    Dp          # Buses: active demand
    Dq          # Buses: reactive demand
    Rp          # Buses: wind active gen
    Rq          # Buses: wind reactive gen
    Pmax        # Buses: max. active gen
    Pmin        # Buses: min. active gen
    Qmax        # Buses: max. reactive gen
    Qmin        # Buses: min. reactive gen
    Plim        # Lines: flow limit
    Vg          # Buses: nominal voltage
    Vceiling    # Buses: max. voltage
    Vfloor      # Buses: min. voltage
    busIdx      # Buses: index
    N           # Buses: total
    Nr          # Buses: renewable
    Ng          # Buses: conventional
    k           # Buses: participation factors
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
    Gp
end

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

function zero2NaN(x)
    # replace 0 with NaN for plotting purposes
    y = deepcopy(x)
    y[y.==0] = NaN
    return y
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
    psDL.Ng,
    psDL.k)
end