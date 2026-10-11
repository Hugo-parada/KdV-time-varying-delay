function[A2_1,A2_2]=mat3(Nx,Nrho,dx,dt,drho,h)
% Constant-delay transport matrix.

Ny=Nx-2;
q=dt/(h*drho);
A2_2=diag((1+q)*ones(Nrho,1));
if Nrho>1;
    A2_2=A2_2+diag(-q*ones(Nrho-1,1),-1);
end
A2_1=zeros(Nrho,Ny);
A2_1(1,1)=-q/dx;
end
