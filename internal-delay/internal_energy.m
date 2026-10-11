function E = internal_energy(Y,Z,tau_t,xi,dx,drho)
%INTERNAL_ENERGY Discrete physical energy used for the internal delay.
Nrho = size(Z,2);
Zfull = [Y,Z];
wrho = ones(Nrho+1,1);
wrho([1,end]) = 0.5;
int_rho = drho*(Zfull.^2)*wrho;
E = 0.5*dx*sum(Y.^2) + 0.5*tau_t*dx*sum(xi.*int_rho);
end
