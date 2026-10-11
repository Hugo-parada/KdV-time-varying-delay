function E = boundary_energy(Y,Z,zrho0,alpha,beta,tau_t,dx,x,rho)
%BOUNDARY_ENERGY Discrete form of
%  1/2 int_0^L y^2 dx + |beta| tau(t)/2 int_0^1 z^2 d rho.

Nx = round(x(end)/dx);
yfull = recover_boundary_state(Y,Z,alpha,beta,dx,Nx);
zfull = [zrho0; Z];

E = 0.5*trapz(x(:),yfull.^2) ...
    + 0.5*abs(beta)*tau_t*trapz(rho(:),zfull.^2);
end
