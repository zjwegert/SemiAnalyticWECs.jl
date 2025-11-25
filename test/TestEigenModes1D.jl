using SemiAnalyticWECs
using LinearAlgebra
using Test
using ForwardDiff

# free-free
N = 10
n = 100
L = 10
x = -L:L/n:L
μ,U,∂ₓ²U = eigenmodes_1d(:free,N,L,x)

# Test orthogonality
δx = L/n
diagws = 2*δx/3*ones(length(x));
diagws[1] = δx/3; diagws[2:2:end-1] .= 4*δx/3; diagws[end] = δx/3;
ws = Diagonal(diagws);
Orthog = U*ws*U'
@test maximum(abs,Orthog - I) < 1e-7
# Check BCs
@test maximum(abs,∂ₓ²U[:,1]) < 1e-10
@test maximum(abs,∂ₓ²U[:,end]) < 1e-10
∂ₓ³U = ForwardDiff.jacobian(x->eigenmodes_1d(:free,N,L,x)[3],x)
@test maximum(abs,∂ₓ³U[:,1]) < 1e-10
@test maximum(abs,∂ₓ³U[:,end]) < 1e-10
# Check PDE
for i ∈ 1:N
  ∂ₓ⁴U = diag(ForwardDiff.jacobian(x->diag(ForwardDiff.jacobian(x->eigenmodes_1d(:free,N,L,x)[3][i,:],x)),x))
  @test ∂ₓ⁴U ≈ U[i,:]*μ[i]^4
end

# clamped-clamped
N = 10
n = 100
L = 10
x = -L:L/n:L
μ,U,∂ₓ²U,α_hat = eigenmodes_1d(:clamped,N,L,x;debug=true)

# Test orthogonality
δx = L/n
diagws = 2*δx/3*ones(length(x));
diagws[1] = δx/3; diagws[2:2:end-1] .= 4*δx/3; diagws[end] = δx/3;
ws = Diagonal(diagws);
Orthog = U*ws*U'
@test maximum(abs,Orthog - I) < 1e-7
# Check BCs
u(x,μ,α) = @. α[1]*exp(μ*x) + α[2]*exp(-μ*x) + α[3]*cos(μ*x) + α[4]*sin(μ*x);
@test maximum(abs,U[:,1]) < 1e-14
@test maximum(abs,U[:,end]) < 1e-14
∂ₓU = map(i->map(y->ForwardDiff.derivative(x->u(x,μ[i],α_hat[i]),y),x),eachindex(μ))
∂ₓU = reduce(hcat,∂ₓU)'
@test maximum(abs,∂ₓU[:,1]) < 1e-14
@test maximum(abs,∂ₓU[:,end]) < 1e-14
# Check PDE
∂ₓ⁴U = map(i->map(y->ForwardDiff.derivative(x->ForwardDiff.derivative(x->ForwardDiff.derivative(x->ForwardDiff.derivative(x->u(x,μ[i],α_hat[i]),x),x),x),y),x),eachindex(μ))
for i ∈ 1:N
  @test ∂ₓ⁴U[i] ≈ U[i,:]*μ[i]^4
end