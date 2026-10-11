function [t,E] = solve_boundary_delay(alpha,beta,tau,dtau,opts)
%SOLVE_BOUNDARY_DELAY Finite-difference solver for the KdV equation with
%time-dependent delayed boundary feedback.
%
% z(rho,t)=y_x(0,t-tau(t)*rho) satisfies
% tau(t) z_t+(1-dtau(t)*rho) z_rho=0, z(0,t)=y_x(0,t).

if nargin < 5 || isempty(opts)
    opts = struct();
end
L = get_option(opts,'L',1);
T = get_option(opts,'T',10);
Nx = get_option(opts,'Nx',100);
Nrho = get_option(opts,'Nrho',100);
Nt = get_option(opts,'Nt',100);
if Nx < 4 || Nrho < 1 || Nt < 1
    error('Nx must be >= 4 and Nrho,Nt must be positive integers.');
end

x = linspace(0,L,Nx+1);
rho = linspace(0,1,Nrho+1);
t = linspace(0,T,Nt+1);
dx = L/Nx; drho = 1/Nrho; dt = T/Nt; Ny = Nx-2;

y0 = @(xx) 0.5*(1-cos(2*pi*xx));
z0 = @(rr) -0.5*sin(2*pi*rr);
Y = y0(x(2:Nx-1)).';
Z = z0(rho(2:end)).';

[Akdv,Cdelay] = build_boundary_kdv_matrix(Nx,Nrho,dx,dt,alpha,beta);
E = zeros(1,Nt+1);
E(1) = boundary_energy(Y,Z,0,alpha,beta,tau(0),dx,x,rho);

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

    [Azy,Azz] = build_boundary_transport_matrix(Nrho,Ny,dt,dx,drho, ...
        tau_new,dtau_new,rho);
    A = [Akdv,Cdelay;Azy,Azz];

    yfull = recover_boundary_state(Y,Z,alpha,beta,dx,Nx);
    yplus = yfull(3:Nx);
    nonlinear = Y.*((yplus-Y)/dx);

    V = A\[Y-dt*nonlinear;Z];
    Y = V(1:Ny);
    Z = V(Ny+1:end);
    E(n+1) = boundary_energy(Y,Z,Y(1)/dx,alpha,beta,tau_new,dx,x,rho);
end
end

function value = get_option(opts,name,default_value)
if isfield(opts,name)
    value = opts.(name);
else
    value = default_value;
end
end
