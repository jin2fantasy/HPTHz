#=
    Find local extrema in Vector and return its index

    Date:2015-07-02
=#

function find_localextrema(A::Vector{Float64})
    difference = diff(A)
    maxima = Int[]
    minima = Int[]

    for i = 1:length(difference) - 1
        prev = difference[i]
        next = difference[i+1]
        if prev * next < 0
            if prev > next
                push!(maxima, i+1)
            else
                push!(minima, i+1)
            end
        end
    end

    return maxima, minima
end
