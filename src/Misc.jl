
function printInfo(simNum,simTotal,btype,N,kₛ,η,kf,M,D)
    println(@sprintf "----------------------------- Simulation %d/%d Complete -----------------------------" simNum simTotal)
    println(@sprintf "Boundary Type: %s, Cell count: %d, Springs per cell: %d" btype N Int(M/N))
    println(@sprintf "kₛ¹: %.5f, η¹: %.5f, kf¹: %.5f, Diffusivity: %.5f" kₛ η kf D)
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

function nonLinearRange(start, stop, length, dist_type)
    linear_range = LinRange(0, 1, length)

    # Applying different distribution types
    if dist_type == "exp"
        # Exponential scaling
        return start .+ (exp.(linear_range .* log(1 + stop - start)) .- 1)
    elseif dist_type == "sine"
        # Sinusoidal scaling
        return start .+ (sin.((π/2) .* linear_range) .* (stop - start))
    elseif dist_type == "cosine"
        # Cosine scaling
        return start .+ ((1 .- cos.((π/2) .* linear_range)) .* (stop - start))
    elseif dist_type == "quad"
        # Quadratic scaling
        return start .+ (linear_range .^ 2 .* (stop - start))
    elseif dist_type == "cubic"
        # Cubic scaling
        return start .+ (linear_range .^ 3 .* (stop - start))
    elseif dist_type == "sigmoid"
        # Sigmoid scaling
        linear_range = LinRange(-1, 1, length)
        k = 5; # slope steepness
        sigmoid_range = 1 ./ (1 .+ exp.(-k.*(linear_range)))
        return start .+ (sigmoid_range .* (stop - start))
    else
        error("Unsupported distribution type")
    end
end
