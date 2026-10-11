function[t,E]=kdv_time_delay(alpha,beta,d)
% Nonlinear KdV with time-varying delayed boundary feedback.

L=1;
T=10;
tau=@(t) d*(1.5+sin(t));
dtau=@(t) d*cos(t);

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

E=zeros(1,Nt+1);

% initial energy
Yfull=zeros(Nx+1,1);
Yfull(2:Nx-1)=Y0;
Yfull(Nx)=-alpha*Y0(1)-beta*dx*Z0(end);
Zfull=[0;Z0];
E(1)=0.5*trapz(x,Yfull.^2)+0.5*abs(beta)*tau(0)*trapz(rho,Zfull.^2);

warned=0;
for n=1:Nt;
    tnew=t(n+1);

    if tau(tnew)<=0;
        error('tau(t) must stay positive')
    end
    if min(1-dtau(tnew)*rho)<=0 && warned==0;
        warning('1-tau''(t)rho changes sign')
        warned=1;
    end

    [A2_1,A2_2]=mat2(Nx,Nrho,dx,dt,drho,tau,dtau,tnew);
    A=[A1_1 A1_2;A2_1 A2_2];

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
    Zfull=[Y0(1)/dx;Z0];
    E(n+1)=0.5*trapz(x,Yfull.^2)+0.5*abs(beta)*tau(tnew)*trapz(rho,Zfull.^2);
end
end
