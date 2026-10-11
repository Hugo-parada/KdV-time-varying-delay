%SMOKE_TEST Basic consistency checks for both corrected solvers.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'boundary-delay'));
addpath(fullfile(root,'internal-delay'));
opts_b = struct('T',0.5,'Nt',10,'Nx',30,'Nrho',20);
[tb,Eb] = kdv_time_delay(0.1,0.1,0.2,opts_b);
assert(numel(tb)==numel(Eb));
assert(all(isfinite(Eb)) && all(Eb>0));
opts_i = struct('T',0.5,'Nt',10,'Nx',10,'Nrho',10);
[ti,Ei] = kdvint_timedelay(1.5,opts_i);
assert(numel(ti)==numel(Ei));
assert(all(isfinite(Ei)) && all(Ei>0));
disp('Smoke test passed.');
