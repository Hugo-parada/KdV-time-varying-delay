function [Azy,Azz] = build_internal_transport_matrix(Nrho,Ny,dt,drho, ...
    tau_t,dtau_t,rho)
%BUILD_INTERNAL_TRANSPORT_MATRIX Implicit upwind discretization in rho.
% z(x,0,t_{n+1})=y(x,t_{n+1}) is imposed in the coupled solve.
rho_nodes = rho(2:end).';
q = dt*(1-dtau_t*rho_nodes)/(tau_t*drho);
Trho = diag(1+q);
if Nrho > 1
    Trho = Trho+diag(-q(2:end),-1);
end
Azz = kron(sparse(Trho),speye(Ny));
Azy = sparse(Ny*Nrho,Ny);
Azy(1:Ny,:) = -q(1)*speye(Ny);
end
