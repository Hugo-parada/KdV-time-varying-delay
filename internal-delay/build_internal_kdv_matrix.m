function Akdv = build_internal_kdv_matrix(Nx,dx,dt,a)
%BUILD_INTERNAL_KDV_MATRIX I+dt*(D_x+D_x^+D_x^+D_x^-+a(x)).
Ny = Nx-2;
a1 = 1+3*dt/dx^3;
a2 = -3*dt/dx^3+dt/(2*dx);
a3 = dt/dx^3;
a4 = -dt/dx^3-dt/(2*dx);
Akdv = diag(a1*ones(Ny,1)) ...
    + diag(a2*ones(Ny-1,1),1) ...
    + diag(a4*ones(Ny-1,1),-1);
if Ny >= 3
    Akdv = Akdv+diag(a3*ones(Ny-2,1),2);
end
Akdv = sparse(Akdv+dt*diag(a));
end
