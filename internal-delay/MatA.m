function[A]=MatA(Nx,dx,dt,a)
% Implicit KdV matrix with instantaneous internal damping.

Ny=Nx-2;
a1=1+3*dt/(dx^3);
a2=-3*dt/(dx^3)+dt/(2*dx);
a3=dt/(dx^3);
a4=-a3-dt/(2*dx);

A=diag(a1*ones(Ny,1))+diag(a2*ones(Ny-1,1),1)...
  +diag(a4*ones(Ny-1,1),-1);
if Ny>=3;
    A=A+diag(a3*ones(Ny-2,1),2);
end
A=A+dt*diag(a);
end
