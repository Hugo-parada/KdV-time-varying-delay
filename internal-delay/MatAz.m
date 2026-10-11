function[Azy,Azz]=MatAz(Nrho,Ny,dt,drho,tau,dtau,tnew,rho)
% Implicit upwind matrix for the internal memory variable.
% Z stores rho_1,...,rho_Nrho at all interior x nodes.

rhon=rho(2:end);
q=dt*(1-dtau(tnew)*rhon)/(tau(tnew)*drho);

Trho=diag(1+q);
if Nrho>1;
    Trho=Trho+diag(-q(2:end),-1);
end

Azz=kron(Trho,eye(Ny));
Azy=zeros(Ny*Nrho,Ny);
Azy(1:Ny,:)=-q(1)*eye(Ny);
end
