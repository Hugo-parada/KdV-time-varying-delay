function [t,E] = kdv_cte_delay(alpha,beta,h,opts)
%KDV_CTE_DELAY Boundary constant-delay comparison tau(t)=h.
if nargin < 4
    opts = struct();
end
if h <= 0
    error('h must be positive.');
end
tau = @(tt) h+0*tt;
dtau = @(tt) 0*tt;
[t,E] = solve_boundary_delay(alpha,beta,tau,dtau,opts);
end
