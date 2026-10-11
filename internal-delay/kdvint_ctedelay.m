function [t,E] = kdvint_ctedelay(h,opts)
%KDVINT_CTEDELAY Internal constant-delay comparison tau(t)=h.
if nargin < 2
    opts = struct();
end
if h <= 0
    error('h must be positive.');
end
tau = @(tt) h+0*tt;
dtau = @(tt) 0*tt;
[t,E] = solve_internal_delay(tau,dtau,opts);
end
