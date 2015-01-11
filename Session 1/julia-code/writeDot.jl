function writeDot(name, busIdx, busInj, renGen, f, t, lineFlow, lineLim, size=(11,14))
    # This function generates a graph that richly expresses the RTS96 system state.
    # name              a name for the graph and resulting dot file
    # busIdx            bus names (could be text) in order of busInj
    # busInj            injection at each bus
    # renGen            renewable generation at each bus (0 for non-wind buses)
    # f                 "from" node for each line
    # t                 "to" node for each line
    # lineFlow          flow on each line
    # lineLim           limit for each line
    # size              size of graphical output

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

    println("$(name).dot generated.\nUse \"run(`dot -Tsvg $(name).dot -o $(name).svg`)\" to create an SVG.")
end