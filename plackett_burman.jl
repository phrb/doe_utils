function binary_search(array::Array{Int, 1}, target::Int)
    n::Int = length(array)

    left::Int = 1
    right::Int = n

    @inbounds while left <= right
        middle = (left + right) >>> 1

        if array[middle] < target
            left = middle + 1
        elseif array[middle] > target
            right = middle - 1
        else
            return true
        end
    end

    return false
end

function paley(matrix::Matrix{Int})
    prime::Int = size(matrix, 1)

    residues_length::Int = 1 + floor((prime - 1) / 2)    
    residues::Array{Int, 1} = Array{Int, 1}(residues_length)

    residues[1] = 0
    residues[2] = 1

    @inbounds for i = 2:(residues_length - 1)
        residues[i + 1] = (i * i) % prime
    end

    sort!(residues)

    @inbounds for i = 1:prime
       @inbounds for j = 1:prime
           offset_value::Int = (i + j - 1) % prime
           if binary_search(residues, offset_value)
               matrix[i, j] = -1
           else
               matrix[i, j] = 1
           end
       end
    end

    return matrix
end
