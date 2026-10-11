function[A1_1,A1_2]=mat1(Nx,Nrho,dx,dt,alpha,beta)
% Matrix for the implicit linear KdV part.
 
Ny=Nx-2;
a1=1+3*dt/(dx^3);
a2=-3*dt/(dx^3)+dt/(2*dx);
a3=dt/(dx^3);
a4=-a3-dt/(2*dx);

A1_1=diag(a1*ones(Ny,1))+diag(a2*ones(Ny-1,1),1)...
    +diag(a4*ones(Ny-1,1),-1);
if Ny>=3;
    A1_1=A1_1+diag(a3*ones(Ny-2,1),2);
end

% y_{Nx-1}=-alpha*y_1-beta*dx*z(1)
A1_1(Ny-1,1)=A1_1(Ny-1,1)-alpha*a3;
A1_1(Ny,1)=A1_1(Ny,1)-alpha*a2;

A1_2=zeros(Ny,Nrho);
A1_2(Ny-1,Nrho)=-beta*dx*a3;
A1_2(Ny,Nrho)=-beta*dx*a2;
end
