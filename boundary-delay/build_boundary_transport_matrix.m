function [Azy,Azz] = build_boundary_transport_matrix(Nrho,Ny,dt,dx,drho, ...
    tau_t,dtau_t,rho)
%BUILD_BOUNDARY_TRANSPORT_MATRIX Implicit upwind discretization in rho.
%
% (1+q_k) z_k^{n+1}-q_k z_{k-1}^{n+1}=z_k^n,
% q_k=dt*(1-dtau(t_{n+1})*rho_k)/(tau(t_{n+1})*drho).
% At k=1, z_0^{n+1}=y_1^{n+1}/dx.

rho_nodes = rho(2:end).';
q = dt*(1-dtau_t*rho_nodes)/(tau_t*drho);

Azz = diag(1+q);
if Nrho > 1
    Azz = Azz + diag(-q(2:end),-1);
end
Azz = sparse(Azz);

Azy = sparse(Nrho,Ny);
Azy(1,1) = -q(1)/dx;
end
