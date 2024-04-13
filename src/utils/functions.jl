function lu_descomposition(A::Tridiagonal{Float64, Vector{Float64}})
    
    n = length(A[1,:])
    L = zeros(n,n)
    U = zeros(n,n)

    
    L[1,1] = A[1,1]
    U[1,1] = 1
    U[1,2] = A[1,2]/L[1,1]

    for i in 2:n

        if i <n
            U[i,i] = 1
            L[i,i-1] = A[i-1,i]
            L[i,i] = A[i,i] - L[i,i-1]*U[i-1,i]
            U[i,i+1] =A[i,i+1]/L[i,i]
        end

        if i==n
            U[i,i] = 1
            L[i,i-1] = A[i-1,i]
            L[i,i] = A[i,i] - L[i,i-1]*U[i-1,i]
        end

    end
    
    return L, U
end

function crow_method(A::Tridiagonal{Float64, Vector{Float64}},b::Vector{Float64})
    
    n = length(b)
    L,U = lu_descomposition(A)
    
    z = zeros(n)
    z[1] = b[1]/L[1,1]

    for i in 2:n
        z[i] = b[i]/L[i,i] - L[i,i-1]*z[i-1]/L[i,i]
    end

    x = zeros(n)
    x[n] = z[n]

    for i in n:-1:2
        x[i-1] = z[i-1] - U[i-1,i]*x[i]
    end

    return x

end