function[t,E]=kdvint_timedelay(M)
% KdV with localized time-varying internal delay.
% M is the maximum delay, so tau(t)=(M-0.5)+0.5*sin(t).

L=1;
T=20;
tau=@(t) (M-0.5)+0.5*sin(t);
dtau=@(t) 0.5*cos(t);

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

% a=2, b=1 and xi=2.1 on (0,L/2)
a=zeros(Ny,1);b=zeros(Ny,1);xi=zeros(Ny,1);
ind=find(xint<L/2+10*eps);
a(ind)=2;b(ind)=1;xi(ind)=2.1;

A1=MatA(Nx,dx,dt,a);
B=zeros(Ny,Ny*Nrho);
B(:,(Nrho-1)*Ny+1:Nrho*Ny)=dt*diag(b);

E=zeros(1,Nt+1);
Zfull=[Y0 Z0];
Erho=trapz(rho,Zfull.^2,2);
E(1)=0.5*dx*sum(Y0.^2)+0.5*tau(0)*dx*sum(xi.*Erho);

warned=0;
for n=1:Nt;
    tnew=t(n+1);
    if tau(tnew)<=0;
        error('tau(t) must stay positive')
    end
    if min(1-dtau(tnew)*rho)<=0 && warned==0;
        warning('1-tau''(t)rho changes sign: outside the upwind regime')
        warned=1;
    end

    [A2_1,A2_2]=MatAz(Nrho,Ny,dt,drho,tau,dtau,tnew,rho);
    A=[A1 B;A2_1 A2_2];

    Yplus=[Y0(2:end);0];
    NL=Y0.*(Yplus-Y0)/dx;

    V1=A\[Y0-dt*NL;Z0(:)];
    Y0=V1(1:Ny);
    Z0=reshape(V1(Ny+1:end),Ny,Nrho);

    Zfull=[Y0 Z0];
    Erho=trapz(rho,Zfull.^2,2);
    E(n+1)=0.5*dx*sum(Y0.^2)+0.5*tau(tnew)*dx*sum(xi.*Erho);
end
end
