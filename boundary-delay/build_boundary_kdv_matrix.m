function [Akdv,Cdelay] = build_boundary_kdv_matrix(Nx,Nrho,dx,dt,alpha,beta)
%BUILD_BOUNDARY_KDV_MATRIX Implicit matrix for D_x+D_x^+D_x^+D_x^-.
%
% Unknowns are y_1,...,y_{Nx-2}. The right boundary condition gives
%   y_{Nx-1} = -alpha*y_1-beta*dx*z_{Nrho}.

Ny = Nx-2;

a1 = 1+3*dt/dx^3;
a2 = -3*dt/dx^3+dt/(2*dx);
a3 = dt/dx^3;
a4 = -dt/dx^3-dt/(2*dx);

Akdv = diag(a1*ones(Ny,1)) ...
    + diag(a2*ones(Ny-1,1),1) ...
    + diag(a4*ones(Ny-1,1),-1);
if Ny >= 3
    Akdv = Akdv + diag(a3*ones(Ny-2,1),2);
end

Akdv(Ny-1,1) = Akdv(Ny-1,1)-alpha*a3;
Akdv(Ny,1)   = Akdv(Ny,1)-alpha*a2;

Cdelay = sparse(Ny,Nrho);
Cdelay(Ny-1,Nrho) = -beta*dx*a3;
Cdelay(Ny,Nrho)   = -beta*dx*a2;

Akdv = sparse(Akdv);
end
