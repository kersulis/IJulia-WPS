function lagrange_activeFlow_droop(psDL)
    # This block of code handles the case with "lagrange", "activeFlow", and "droop" settings.
    
    # unpack psDL (boilerplate):
    (Sb,f,t,r,x,b,Y,bustype,
    Gp,Gq,Dp,Dq,Rp,Rq,
    Pmax,Pmin,Qmax,Qmin,Plim,
    Vg,Vceiling,Vfloor,
    busIdx,N,Nr,Ng,k) = unpack_psDL(psDL)
    
    # Prepare all variables:
    sorted = [find(Rp); setdiff(1:N,find(Rp))]
    ρ0 = Rp[find(Rp)]
    Λ = speye(Nr)
    Nc = N - Nr

    Yr = Y[sorted[1:Nr],:]
    Yc = Y[sorted[Nr+1:end],:]

    G0 = Gp
    k = Float64[]
    for i = 1:length(psData.Gp)
        if psData.Gp[i] != 0
            push!(k, 1/length(find(psData.Gp)))
        else
            push!(k,0)
        end
    end

    kr = k[sorted[1:Nr]]
    kc = k[sorted[Nr+1:end]]

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
    Gpost = Array(Vector{Float64},0) # array of vectors with Float64 values

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

            # B = sparse(B)

            X = instInvert(A,B,i)

            rho = Float64[]
            j = 1
            for i=1:N
                if Rp[i] != 0
                    push!(rho, X[j])
                    j += 1
                else 
                    push!(rho,0)
                end
            end

            push!(ρ, rho)
            push!(θ, X[Nr+1:Nr+N]) # angles sorted by original bus numbers
            push!(α, X[Nr+N+1])

        if sum(abs(X[1:Nr])) == 0
            push!(score, NaN)
        else
            push!(score, [0.5(X[1:Nr] - ρ0)'*(X[1:Nr] - ρ0)][1])
        end
        # Compute conventional generation post-instanton:
        push!(Gpost, G0 + k.*(X[Nr+N+1]))
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
    rankedList = sortrows([zero2NaN(score) [ρ θ α constrIdx Gpost]]) # sort candidates; place NaNs at bottom
    
    # Pack results into instance of "instantonResults":
    return instantonResults(rankedList[:,1],rankedList[:,2],rankedList[:,3],rankedList[:,4],rankedList[:,5],rankedList[:,6])
end