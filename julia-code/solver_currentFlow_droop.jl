function solver_currentFlow_droop(psDL)
    # This function uses MOSEK to perform instanton analysis for the case with conventional generator droop resopnse.

    # unpack psDL (boilerplate):
    (Sb,f,t,r,x,b,Y,bustype,
    Gp,Gq,Dp,Dq,Rp,Rq,
    Pmax,Pmin,Qmax,Qmin,Plim,
    Vg,Vceiling,Vfloor,
        busIdx,N,Nr,Ng,k) = unpack_psDL(psDL)

    score = Float64[]
    rho = Array(Vector{Float64},0) # array of vectors with Float64 values
    Gpost = Array(Vector{Float64},0) # array of vectors with Float64 values
    theta = Array(Vector{Float64},0) # array of vectors with Float64 values
    alpha = Float64[]

    result = Symbol[]

    # Enforce each constraint in the f --> t direction
    for idx = 1:length(f)
        (m, ρ, θ, α, PowerBalance, AngleRef, Congestion, Participation) = droopCurrentModel(idx, 1, Rp, Gp, Dp, f, t, x, Y, bustype, Plim,k,Vg)
        status = solve(m);
        push!(score, getObjectiveValue(m))
        push!(result,status)
        push!(theta, getValue(θ)[:])
        push!(alpha, getValue(α))
        push!(rho,getValue(ρ)[:])
  
        # Compute conventional generation post-instanton:
        push!(Gpost, Gp + k.*getValue(α))
    end

    # Enforce each constraint in the t --> f direction
    for idx = 1:length(f)
        (m, ρ, θ, α, PowerBalance, AngleRef, Congestion, Participation) = droopCurrentModel(idx, -1, Rp, Gp, Dp, f, t, x, Y, bustype, Plim,k,Vg)
        status = solve(m);
        push!(score, getObjectiveValue(m))
        push!(result, status)
        push!(theta, getValue(θ)[:])
        push!(alpha, getValue(α))
        push!(rho, getValue(ρ)[:])

        # Compute conventional generation post-instanton:
        push!(Gpost, Gp + k.*getValue(α))
    end

    # Generate strings to tell user which constraint was violated:
    constrIdx = String[]
    for i = 1:2length(f)
        if i <= length(f)
            push!(constrIdx, "activeFlow: $(f[i]) -> $(t[i]) @ $(Plim[i])")
        else
            i2 = i - length(f)
            push!(constrIdx, "activeFlow: $(t[i2]) -> $(f[i2]) @ $(Plim[i2])")
        end
    end

    rankedList = sortrows([zero2NaN(score) [rho theta alpha constrIdx Gpost]]) # sort candidates; place NaNs at bottom

    # Pack results into instance of "instantonResults":
    return instantonResults(rankedList[:,1],rankedList[:,2],rankedList[:,3],rankedList[:,4],rankedList[:,5],rankedList[:,6])
end

function droopCurrentModel(idx, sense, Rp, Gp, Dp, f, t, x, Y, bustype, Plim,k, Vg)
    # 7-29-14
    # DROOP RESPONSE
    # Create model saturating line 'idx' in direction 'sense' (±1)
    # This function uses JuMP and Mosek
    # Line current constraints
    
    m = Model(solver = MosekSolver()) # Use MOSEK to solve model
    N = length(Rp)
    @defVar(m, ρ[1:N] >= 0) # Add decision variables to model (renewable gen)
    @defVar(m, θ[1:N]) # Add bus angles
    @defVar(m, α) # mismatch parameter
    setObjective(m, :Min, 0.5*sum([(ρ[i] - Rp[i])^2 for i in find(Rp)]))

    # add power balance constraints:
    @defConstrRef PowerBalance[1:N]
    for i in setdiff(1:N,find(bustype.==3))
        PowerBalance[i] = @addConstraint(m, sum([Y[i,j]*θ[j] for j in 1:N]) == Gp[i] + k[i]*α + ρ[i] - Dp[i])
    end

    # wind output must be 0 at non-wind buses:
    @addConstraint(m, NonWind, sum([ρ[i] for i in setdiff(collect(1:N),find(Rp))]) == 0)

    # angle reference bus must have angle 0:
    @addConstraint(m, AngleRef, θ[find(bustype.==3)[1]] == 0)
 
    # 'i' is f[idx]
    # 'k' is t[idx]
    #i = f[idx]
    #k = t[idx]

    Vi = Vg[f[idx]]
    Vk = Vg[t[idx]]

    # saturate a line (current constraint):
    @addConstraint(m, Congestion, θ[f[idx]] - θ[t[idx]] == sense*acos((Vi^2 + Vk^2 - (Plim[idx]*x[idx])^2)/(2Vi*Vk)))
    
    #@addConstraint(m, Congestion, θ[f[idx]] - θ[t[idx]] == sense*x[idx]*Plim[idx]) # saturate a line

    # @addConstraint(m, Congestion, (Vi^2 + Vk^2 - 2Vi*Vk*sense*cos(θ[i] - θ[k]))/(x[idx]^2) == Plim[idx]^2)

    @addConstraint(m, Participation, α == (sum(Dp) - sum([ρ[i] for i in find(Rp)])) - sum(Gp))
    return m, ρ, θ, α, PowerBalance, AngleRef, Congestion, Participation
end