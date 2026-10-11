function[Azy,Azz]=MatAz_cte(Nrho,Ny,dt,drho,h)
% Constant-delay transport matrix for the internal memory variable.

q=dt/(h*drho);
Trho=diag((1+q)*ones(Nrho,1));
if Nrho>1;
    Trho=Trho+diag(-q*ones(Nrho-1,1),-1);
end
Azz=kron(Trho,eye(Ny));
Azy=zeros(Ny*Nrho,Ny);
Azy(1:Ny,:)=-q*eye(Ny);
end
