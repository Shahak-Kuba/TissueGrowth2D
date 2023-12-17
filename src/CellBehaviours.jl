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

function cell_probs(uᵢ,m,δt)
    ρ = calc_cell_densities(uᵢ,m)
    return (P(true,ρ,0.01).*δt, A(false, ρ,0.02).*δt, E(true, ρ,0.03).*δt)
end