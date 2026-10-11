function [t,E] = solve_internal_delay(tau,dtau,opts)
%SOLVE_INTERNAL_DELAY Finite-difference solver for delayed internal feedback.
if nargin < 3 || isempty(opts)
    opts = struct();
end
L = get_option(opts,'L',1);
T = get_option(opts,'T',20);
Nx = get_option(opts,'Nx',10);
Nrho = get_option(opts,'Nrho',20);
Nt = get_option(opts,'Nt',100);
a_value = get_option(opts,'a_value',2);
b_value = get_option(opts,'b_value',1);
xi_value = get_option(opts,'xi_value',2.1);
omega_end = get_option(opts,'omega_end',L/2);
if Nx < 4 || Nrho < 1 || Nt < 1
    error('Nx must be >= 4 and Nrho,Nt must be positive integers.');
end

x = linspace(0,L,Nx+1);
xint = x(2:Nx-1).';
rho = linspace(0,1,Nrho+1);
t = linspace(0,T,Nt+1);
dx = L/Nx; drho = 1/Nrho; dt = T/Nt; Ny = Nx-2;

mask = xint <= omega_end+10*eps(L);
a = zeros(Ny,1); b = zeros(Ny,1); xi = zeros(Ny,1);
a(mask) = a_value; b(mask) = b_value; xi(mask) = xi_value;

Y = 1-cos(2*pi*xint);
Z = Y*cos(2*pi*rho(2:end));
Akdv = build_internal_kdv_matrix(Nx,dx,dt,a);
Cdelay = sparse(Ny,Ny*Nrho);
Cdelay(:,(Nrho-1)*Ny+(1:Ny)) = dt*spdiags(b,0,Ny,Ny);

E = zeros(1,Nt+1);
E(1) = internal_energy(Y,Z,tau(0),xi,dx,drho);
warned = false;
for n = 1:Nt
    tnew = t(n+1);
    tau_new = tau(tnew);
    dtau_new = dtau(tnew);
    if tau_new <= 0
        error('The delay tau(t) must stay strictly positive.');
    end
    if any(1-dtau_new*rho(2:end) <= 0) && ~warned
        warning(['The condition 1-dtau(t)*rho>0 is violated. ', ...
            'This run is outside the theoretical/upwind regime used in the paper.']);
        warned = true;
    end

    [Azy,Azz] = build_internal_transport_matrix(Nrho,Ny,dt,drho, ...
        tau_new,dtau_new,rho);
    A = [Akdv,Cdelay;Azy,Azz];
    yplus = [Y(2:end);0];
    nonlinear = Y.*((yplus-Y)/dx);
    V = A\[Y-dt*nonlinear;Z(:)];
    Y = V(1:Ny);
    Z = reshape(V(Ny+1:end),Ny,Nrho);
    E(n+1) = internal_energy(Y,Z,tau_new,xi,dx,drho);
end
end

function value = get_option(opts,name,default_value)
if isfield(opts,name)
    value = opts.(name);
else
    value = default_value;
end
end
