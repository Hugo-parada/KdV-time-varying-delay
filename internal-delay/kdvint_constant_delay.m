function[t,E]=kdvint_constant_delay(h)
% KdV with localized constant internal delay.

L=1;
T=20;
Nx=10;
dx=L/Nx;
x=0:dx:L;
xint=x(2:Nx-1)';
Ny=Nx-2;

Nrho=20;
drho=1/Nrho;
rho=0:drho:1;
Nt=100;
dt=T/Nt;
t=0:dt:T;

Y0=1-cos(2*pi*xint);
Z0=Y0*cos(2*pi*rho(2:end));

a=zeros(Ny,1);b=zeros(Ny,1);xi=zeros(Ny,1);
ind=find(xint<L/2+10*eps);
a(ind)=2;b(ind)=1;xi(ind)=2.1;

A1=MatA(Nx,dx,dt,a);
[A2_1,A2_2]=MatAz_constant(Nrho,Ny,dt,drho,h);
B=zeros(Ny,Ny*Nrho);
B(:,(Nrho-1)*Ny+1:Nrho*Ny)=dt*diag(b);
A=[A1 B;A2_1 A2_2];

E=zeros(1,Nt+1);
Zfull=[Y0 Z0];
E(1)=0.5*dx*sum(Y0.^2)+0.5*h*dx*sum(xi.*trapz(rho,Zfull.^2,2));

for n=1:Nt;
    Yplus=[Y0(2:end);0];
    NL=Y0.*(Yplus-Y0)/dx;
    V1=A\[Y0-dt*NL;Z0(:)];
    Y0=V1(1:Ny);
    Z0=reshape(V1(Ny+1:end),Ny,Nrho);
    Zfull=[Y0 Z0];
    E(n+1)=0.5*dx*sum(Y0.^2)+0.5*h*dx*sum(xi.*trapz(rho,Zfull.^2,2));
end
end
