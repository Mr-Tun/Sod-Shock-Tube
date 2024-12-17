function [f] = ReconstructFlux(uj,uj1,cf)
%RECONSTRUCTFLUX 此处显示有关此函数的摘要
%   此处显示详细说明
cf = cf';
gamma = 1.4;
Umid = 0.5*(uj+uj1);             %compute U i+0.5
rho = Umid(1,1);                 %decouple the solution vector U
rhou = Umid(1,2);
rhoe = Umid(1,3);
u = rhou / rho;                  % compute lambda matrix
e = rhoe / rho - 0.5 * u ^ 2;
p = (gamma - 1) * rho * e;
a =  (gamma*p/rho)^0.5;
H = a*a/(gamma-1) + 0.5 * u ^ 2;
R = [-(gamma-1)/a/a,        -1/2/a,                          1/2/a;
       -(gamma-1)/a/a*u,      -(u-a)/2/a,                      (u+a)/2/a;
       -(gamma-1)/a/a*u*u/2,  -(u*u/2+a*a/(gamma-1)-u*a)/2/a,  (u*u/2+a*a/(gamma-1)+u*a)/2/a];    % compute diag matrix R and R^-1
f = R *cf;
f = f'; 
end

