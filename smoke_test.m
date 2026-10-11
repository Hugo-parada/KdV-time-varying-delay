% Short consistency test for the corrected recovered codes.
clear;clc;

root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'boundary-delay'));
addpath(fullfile(root,'internal-delay'));

[t1,E1]=kdv_time_delay(0.1,0.1,0.2);
assert(length(t1)==length(E1));
assert(all(isfinite(E1)));
assert(all(E1>0));

[t2,E2]=kdvint_timedelay(1.5);
assert(length(t2)==length(E2));
assert(all(isfinite(E2)));
assert(all(E2>0));

disp('Smoke test passed.')
