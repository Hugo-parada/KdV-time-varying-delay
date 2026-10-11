function[t,E]=kdv_cte_delay(alpha,beta,h)
% Nonlinear KdV with constant delayed boundary feedback.

L=1;
T=10;
y0=@(x) 0.5*(1-cos(2*pi*x));
z0=@(rho) -0.5*sin(2*pi*rho);

Nx=100;
dx=L/Nx;
x=0:dx:L;
Y0=y0(x(2:Nx-1))';

Nrho=100;
drho=1/Nrho;
rho=0:drho:1;
Z0=z0(rho(2:end))';

Nt=100;
dt=T/Nt;
t=0:dt:T;

[A1_1,A1_2]=mat1(Nx,Nrho,dx,dt,alpha,beta);
[A2_1,A2_2]=mat3(Nx,Nrho,dx,dt,drho,h);
A=[A1_1 A1_2;A2_1 A2_2];

E=zeros(1,Nt+1);
Yfull=zeros(Nx+1,1);
Yfull(2:Nx-1)=Y0;
Yfull(Nx)=-alpha*Y0(1)-beta*dx*Z0(end);
E(1)=0.5*trapz(x,Yfull.^2)+0.5*abs(beta)*h*trapz(rho,[0;Z0].^2);

for n=1:Nt;
    Yfull=zeros(Nx+1,1);
    Yfull(2:Nx-1)=Y0;
    Yfull(Nx)=-alpha*Y0(1)-beta*dx*Z0(end);
    Yplus=Yfull(3:Nx);
    NL=Y0.*(Yplus-Y0)/dx;

    V1=A\[Y0-dt*NL;Z0];
    Y0=V1(1:Nx-2);
    Z0=V1(Nx-1:end);

    Yfull=zeros(Nx+1,1);
    Yfull(2:Nx-1)=Y0;
    Yfull(Nx)=-alpha*Y0(1)-beta*dx*Z0(end);
    E(n+1)=0.5*trapz(x,Yfull.^2)+0.5*abs(beta)*h*trapz(rho,[Y0(1)/dx;Z0].^2);
end
end
