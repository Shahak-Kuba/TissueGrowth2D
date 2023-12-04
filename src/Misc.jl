
function printInfo(simNum,simTotal,btype,N,kₛ,η,kf)
    println(@sprintf "----------------------------- Simulation %d/%d Complete -----------------------------" simNum simTotal)
    println(@sprintf "Boundary Type: %s, Cell count: %d, kₛ: %.5f, η: %.5f, kf: %.5f" btype N kₛ η kf)
    println(@sprintf "-----------------------------------------------------------------------------------")
end

function SaveData(data, SaveName, SaveFolder)
    # Check if the folder exists, create it if it doesn't
    if !isdir(SaveFolder)
        mkdir(SaveFolder)
    end
    
    # Construct the file path
    filepath = joinpath(SaveFolder, "$SaveName.jld2")
    
    # Check if the file already exists
    if isfile(filepath)
        error("File '$SaveName.jld2' already exists in folder '$SaveFolder'. Please choose a different name.")
    else
        # Save the data to the specified file
        save(filepath, "data", data)
    end
end

function LoadData(SaveName, SaveFolder)
    # Construct the file path
    filepath = joinpath(SaveFolder, "$SaveName.jld2")
    
    # Load the data from the specified file
    loaded_data = load(filepath)
    
    # Retrieve the data from the loaded file
    if haskey(loaded_data, "data")
        return loaded_data["data"]
    else
        error("No 'data' key found in the loaded file.")
    end
end