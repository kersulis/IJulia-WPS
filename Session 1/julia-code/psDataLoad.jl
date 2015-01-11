# Function used to create instance of powerSystemData type from MATLAB workspace:
function psDataLoad()
    (Sb, f, t, r, x, b, Y, bustype, 
     Gp, Gq, Dp, Dq, Rp, Rq, 
     Pmax, Pmin, Qmax, Qmin, 
     Plim, Vg, Vceiling, Vfloor, 
    busIdx, N, Nr, Ng) = readRTS96Data()

    # Allow each generator to participate equally in droop response.
    # Note: this only applies to analysis types with droop response!
    k = Float64[]
    for i = 1:length(Gp)
        if Gp[i] != 0
            push!(k, 1/length(find(Gp)))
        else
            push!(k,0)
        end
    end

    return powerSystemData(Sb,f,t,r,x,b,Y,bustype,Gp,Gq,Dp,Dq,Rp,Rq,Pmax,Pmin,Qmax,Qmin,Plim,Vg,Vceiling,Vfloor,busIdx,N,Nr,Ng,k)
end