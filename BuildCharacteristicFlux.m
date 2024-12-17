function [fj] = BuildCharacteristicFlux(uj,uj1,fj)
%Characteristic function is used to transfer the original flux to the characteristic flux
fj = fj';
gamma = 1.4;
Umid = 0.5*(uj+uj1);             %compute U i+0.5
rho = Umid(1,1);                 %decouple the solution vector U
rhou = Umid(1,2);
rhoe = Umid(1,3);
u = rhou / rho;                  % compute lambda matrix
e = rhoe / rho - 0.5 * u ^ 2;
p = (gamma - 1) * rho * e;
a =  (gamma*p/rho)^0.5;

Rinv=[u*u/2-a*a/(gamma-1),   -u,               1;
        -u-(gamma-1)/a*u*u/2,  1+(gamma-1)/a*u,  -(gamma-1)/a;
        -u+(gamma-1)/a*u*u/2,  1-(gamma-1)/a*u,  (gamma-1)/a];     % compute diag matrix R and R^-1


fj = Rinv * fj;
fj = fj'; 
end

