function [t,E] = kdvint_timedelay(M,opts)
%KDVINT_TIMEDELAY Internal time-varying delay with maximum M.
% tau(t)=(M-0.5)+0.5*sin(t), matching the recovered runs and figure labels.
if nargin < 2
    opts = struct();
end
if M <= 1
    error('M must be > 1 so that the delay stays positive.');
end
tau = @(tt) (M-0.5)+0.5*sin(tt);
dtau = @(tt) 0.5*cos(tt);
[t,E] = solve_internal_delay(tau,dtau,opts);
end
