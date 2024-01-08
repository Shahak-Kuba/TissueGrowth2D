function remove_row!(matrix, row_to_remove)
    return matrix[setdiff(1:size(matrix, 1), row_to_remove), :]
end

function add_row!(matrix, new_row, pos)
    return vcat(matrix[1:pos-1, :], new_row', matrix[pos:end, :])
end

function calc_spring_densities(uᵢ)
    uᵢ₊₁ = circshift(uᵢ,1)
    return 1 ./ δ(uᵢ₊₁, uᵢ)
end

function calc_cell_densities(u,m)
    array = calc_spring_densities(u)
    return [sum(array[i:i+m-1]) for i in 1:m:length(array)-m+1]
end


P(event,ρ,α) = event ? α.*ρ : zeros(size(ρ))
A(event,ρ,β) = event ? β.*ρ : zeros(size(ρ))
E(event,ρ,γ) = event ? γ.*ρ : zeros(size(ρ))

function cell_probs(uᵢ,m,δt,prolif,death,embed,α,β,γ)
    ρ = calc_cell_densities(uᵢ,m)
    return (P(prolif,ρ,α).*δt, A(death, ρ,β).*δt, E(embed, ρ,γ).*δt)
end

function find_cell_index(arr::Vector{Float64}, threshold::Float64)
    cum_sum = cumsum(arr)
    index = findfirst(cum_sum .>= threshold)
    if index === nothing
        return nothing  # If the cumulative sum never reaches the threshold
    end
    return index
end


function event_callback!(integrator, t, u)
    (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ) = integrator.p
    (p,a,e) = cell_probs(u, m, δt, prolif, death, embed, α, β, γ)
    (r1,r2,r3) = rand(3)
    if r1 < (sum(p) + sum(a) + sum(e)) # check if event has occurred
        if r2 < sum(p) / (sum(p) + sum(a) + sum(e)) # prolif occurs
            idx = find_cell_index(p, r3 * sum(p))
            # Perform operations based on prolif occurrence if needed

        elseif sum(p) / (sum(p) + sum(a) + sum(e)) < r2 < (sum(p) + sum(a)) / (sum(p) + sum(a) + sum(e)) # death occurs
            idx = find_cell_index(a, r3 * sum(a))
            integrator.u = remove_row!(integrator.u, idx)
            # Remove the cell from u based on death occurrence

        else # embed occurs
            idx = find_cell_index(e, r3 * sum(e))
            # Perform operations based on embed occurrence if needed
        end
    end
    
    return integrator  # Return the updated u if modified
end

function affect!(integrator)
    println("test")
    (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ) = integrator.p
    u = integrator.u
    u_copy = copy(reshape(integrator.u, 2, Int(length(integrator.u)/2))')

    (p,a,e) = cell_probs(u_copy, m, δt, prolif, death, embed, α, β, γ)
    (r1,r2,r3) = rand(3)
    if r1 < (sum(p) + sum(a) + sum(e)) # check if event has occurred
        if r2 < sum(p) / (sum(p) + sum(a) + sum(e)) # prolif occurs
            idx = find_cell_index(p, r3 * sum(p))
            # Perform operations based on prolif occurrence if needed
 

        elseif sum(p) / (sum(p) + sum(a) + sum(e)) < r2 < (sum(p) + sum(a)) / (sum(p) + sum(a) + sum(e)) # death occurs
            resize!(integrator, length(u)-2)
            idx = find_cell_index(a, r3 * sum(a))
            u .= vec(remove_row!(u_copy, idx)')
            # Remove the cell from u based on death occurrence


        else # embed occurs
            idx = find_cell_index(e, r3 * sum(e))
            # Perform operations based on embed occurrence if needed
        end
    end
    nothing
end