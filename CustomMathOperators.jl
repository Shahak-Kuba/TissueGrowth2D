function row_dot(A,B)
    return sum(abs2, (A.*B), dims=2)
end