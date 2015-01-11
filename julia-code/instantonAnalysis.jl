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
                return solver_activeFlow_droop(psDL)

                elseif genResponse=="proportional"
                # DONE: perform solver analysis with active power flows and proportional response
                return solver_activeFlow_proportional(psDL)

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
                return solver_apparentFlow_droop(psDL)

                elseif genResponse=="proportional"
                # perform solver analysis with apparent power flows and proportional response
                
                elseif genResponse=="none"
                # perform Lagrangian analysis with apparent power flows and no response
                
                else println("genResponse value invalid.")
            end
            
            elseif constraintType=="currentFlow"
                if genResponse=="droop"
                    # perform solver analysis with curret flow limits and droop response
                    return solver_currentFlow_droop(psDL)
                    
                else println("genResponse value invalid")
                end

            else println("constraintType value invalid.")
        end
        else println("method value invalid.")
    end
end
