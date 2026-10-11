function [t,E] = kdv_time_delay(alpha,beta,d,opts)
%KDV_TIME_DELAY Boundary time-varying delay tau(t)=d*(1.5+sin(t)).
if nargin < 4
    opts = struct();
end
if d <= 0
    error('d must be positive.');
end
tau = @(tt) d*(1.5+sin(tt));
dtau = @(tt) d*cos(tt);
[t,E] = solve_boundary_delay(alpha,beta,tau,dtau,opts);
end
