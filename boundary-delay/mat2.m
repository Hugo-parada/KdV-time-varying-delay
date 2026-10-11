function[A2_1,A2_2]=mat2(Nx,Nrho,dx,dt,drho,tau,dtau,tnew)
% Implicit upwind matrix for
% tau(t)z_t+(1-tau'(t)rho)z_rho=0. 

Ny=Nx-2;
rho=(1:Nrho)*drho;
q=dt*(1-dtau(tnew)*rho)/(drho*tau(tnew));

A2_2=diag(1+q);
if Nrho>1;
    A2_2=A2_2+diag(-q(2:end),-1);
end

A2_1=zeros(Nrho,Ny);
% z(0,tnew)=y_x(0,tnew)=y_1(tnew)/dx.
% In the recovered code this coefficient accidentally used k=Nrho.
A2_1(1,1)=-q(1)/dx;
end
