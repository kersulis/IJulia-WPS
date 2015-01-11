function plotInstantonResults(	eventIdx,
								psData,
								results;
								plotType="graph",
								plotName="defaultPlotName",
								constraintType="activeFlow")
	# This function plots instanton results according to user-specified parameters.
	# name 					string; desired output file name
	# eventIdx				index of the extreme event the user wishes to plot
	# psData				instance of powerSystemData type
	# results 				instance of instantonResults type
	# plotType 				can be either "graph" or "bar"
	# constraintType		can be "activeFlow", "apparentFlow", or "voltageMagnitude"

	if constraintType == "activeFlow"
		if plotType == "graph"
			lineFlow = Float64[]
			for i = 1:length(psData.f)
			    push!(lineFlow, (results.θ[eventIdx][psData.f[i]]-results.θ[eventIdx][psData.t[i]])/psData.x[i])
			end

			writeDot(plotName,
			    psData.busIdx,
			    results.Gp[eventIdx] + results.ρ[eventIdx] - psData.Dp,
				results.ρ[eventIdx],
				psData.f,
				psData.t,
				lineFlow,
				psData.Plim,
				(8,12))
		end
	elseif constraintType == "currentFlow"
		if plotType == "graph"
			lineFlow = Float64[]
			for i = 1:length(psData.f)
				Vi = psData.Vg[psData.f[i]]
				Vk = psData.Vg[psData.t[i]]
				θi = results.θ[eventIdx][psData.f[i]]
				θk = results.θ[eventIdx][psData.t[i]]

			    push!(lineFlow, sign(θi-θk)*sqrt( ( Vi^2 + Vk^2 -2Vi*Vk*cos(θi-θk) )/(psData.x[i])^2) )
			end

			writeDot(plotName,
			    psData.busIdx,
			    results.Gp[eventIdx] + results.ρ[eventIdx] - psData.Dp,
				results.ρ[eventIdx],
				psData.f,
				psData.t,
				lineFlow,
				psData.Plim.^2,
				(8,12))
		end

	end

end