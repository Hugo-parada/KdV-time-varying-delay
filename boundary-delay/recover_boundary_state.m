function yfull = recover_boundary_state(Y,Z,alpha,beta,dx,Nx)
%RECOVER_BOUNDARY_STATE Recover all spatial grid values from the unknowns.

yfull = zeros(Nx+1,1);
yfull(2:Nx-1) = Y;
yfull(Nx) = -alpha*Y(1)-beta*dx*Z(end); % y_{Nx-1}
end
