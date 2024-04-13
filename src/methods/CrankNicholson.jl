using  LinearAlgebra
using  Printf, PlotlyJS

include("../utils/functions.jl")

a,b = 0,1
t_min, t_max = 0,0.5

ci = x -> sin.(π*x)

xn = 10
tn  = 50

h = (b-a)/xn
k = (t_max-t_min)/tn 
α = 1

λ = k*α^2/h^2

Adi = fill(-λ/2,xn-2)
Ad  = fill(1+λ,xn-1)
Ads = fill(-λ/2,xn-2)
A  = Tridiagonal(Adi, Ad, Ads)  

Bdi = fill(λ/2,xn-2)
Bd  = fill(1-λ,xn-1)
Bds = fill(λ/2,xn-2)
B  = Tridiagonal(Bdi, Bd, Bds)  

S = zeros(tn+1,xn+1)
S[:,1] = zeros(tn+1)
S[:,xn+1] = zeros(tn+1)
S[1,2:xn] = ci(LinRange(a+h, b-h, xn-1))
S

for i in 2:tn+1
    #S[i,2:xn] = A\(B*S[i-1,2:xn])
    S[i,2:xn] = crow_method(A,B*S[i-1,2:xn])
end

println("x \t\t S(x,t) ")
for (index, value) in enumerate(S[tn+1,:])
    xi =  a + (index-1)*h
    @printf("x_%d = %.2f \t U(%.2f,%.2f) = %.2e \n", index, xi, xi, t_max, value)
end

x = LinRange(a, b, xn+1)
y = LinRange(t_min, t_max, tn+1)
Z = [vec(S[i, :]) for i in 1:length(S[:,1])]
layout = Layout(
    title="MN3 Implicit Finite Difference Method",
    autosize=false,
    width=600,
    height=400,
    margin=attr(l=65, r=50, b=65, t=90),
    scene=attr(
            xaxis_nticks=20,
            zaxis_nticks=4,
            camera_eye=attr(x=0, y=-1.5, z=0.5),
            aspectratio=attr(x=1, y=1, z=0.5)
        )
)
Plot(surface(x=x, y=y, z=Z), layout)