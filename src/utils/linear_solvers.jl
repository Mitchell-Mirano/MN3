using LinearAlgebra

function Jacobi(A,b,x0=Nothing,iters = 10,tol= 0.01)

   
    D = Diagonal(A)
    L = -(LowerTriangular(A) - D )
    U = -(UpperTriangular(A) - D )
    D1 = inv(D)
    T = D1*(L+U)
    C = D1*b 
    xk = zeros(size(A)[1])

    n_iters = 0

    if x0 != Nothing
        xk = x0
    end

    for i in 1:iters
        new_x = T*xk +C
        if norm(new_x-xk) < tol
            xk = new_x
            n_iters +=1
            break
        end
        xk = new_x
        n_iters +=1
    end
    
    return iters,xk
end


function GaussSeidel(A,b,x0=Nothing,iters = 10,tol= 0.01)
    
    D = Diagonal(A)
    L = -(LowerTriangular(A) - D )
    U = -(UpperTriangular(A) - D )
    DL = inv(D-L)
    T = DL*U
    C = DL*b
    xk = zeros(size(A)[1])

    n_iters = 0

    if x0 != Nothing
        xk = x0
    end

    for i in 1:iters
        new_x = T*xk +C
        if norm(new_x-xk) < tol
            xk = new_x
            n_iters += 1
            break
        end
        xk = new_x
        n_iters +=1
    end
    
    return n_iters,xk
end