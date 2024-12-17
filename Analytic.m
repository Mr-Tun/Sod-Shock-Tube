function [space,U] = Analytic(t)
gamma=1.4;
L = 2;                              % the length of the shock tube -m
p0 = 1*10^5;                        % the initial pressure --Pa
rho0 = 1;                           % the initial density   --kg/m³
p4 = 2*p0;                          % the pressure in the region 4
rho4 = 2*rho0;                      % the density in the region 4
u4 = 0;
p1 = p0;                            % the pressure in the region 1
rho1 = rho0;                        % the density in the region 1
u1 = 0;
a1 = sqrt(gamma*p1/rho1);           % sonic velocity in the region 1
a4 = sqrt(gamma*p4/rho4);           % sonic velocity in the region 4
pstar = Newton (p1,p4,rho1,rho4);   % the pressure in the region 3 and 2
g = gamma*pstar/(pstar-p1)-(gamma-1)/2;
shock_wave = sqrt(g*(pstar-p1)/rho1);
wave_head = a4;
p2 = pstar; u2 = 1/g*shock_wave; rho2 = g/(g-1)*rho1; % the state in the region 2
p3= pstar; u3 = u2; rho3 = rho4*(p3/p4)^(1/gamma);    % the state in the region 3 
a3 = sqrt(gamma*p3/rho3);
wave_tail = u3-a3;
N = 2000;
dx = L/N;
space = 0:dx:2;
N1 = round((1-t*abs(wave_head))/dx);
N2 = round((1-t*abs(wave_tail))/dx) -N1;
N3 = round((1+t*abs(u2))/dx)-N2-N1;
N4 = round((1+t*abs(shock_wave))/dx)-N3-N2-N1;
N5 = N+1-N1-N2-N3-N4;
u5 = linspace(0, u3, N2); u5 = u5';
p5 = p4.*(1-(gamma-1)/2*u5/a4).^(2*gamma/(gamma-1));
rho5 = rho4.*(1-(gamma-1)/2*u5/a4).^(2/(gamma-1));
u = [u4*ones(N1,1);u5.*ones(N2,1);u3*ones(N3,1);u2*ones(N4,1);u1*ones(N5,1)];
p = [p4*ones(N1,1);p5.*ones(N2,1);p3*ones(N3,1);p2*ones(N4,1);p1*ones(N5,1)];
rho = [rho4*ones(N1,1);rho5.*ones(N2,1);rho3*ones(N3,1);rho2*ones(N4,1);rho1*ones(N5,1)];
U = [rho, u , p];
end

function [f] = F(pstar,pi,rhoi)
%F function is used to compute f value in the sod- shock tube
gamma=1.4;
if pstar>pi
f=(pstar-pi)/(rhoi*(gamma*pi/rhoi)^0.5*((gamma+1)/2/gamma*pstar/pi+(gamma-1)/2/gamma)^0.5);
else
f=2/(gamma-1)*(gamma*pi/rhoi)^0.5*((pstar/pi)^((gamma-1)/2/gamma)-1);
end
end

function pstar = Newton(p1,p2,rho1,rho2)
%Newton function is used to solve F based on  Newton’s iteration 
pstar = 1.5;                       
pstar_new = 1.5+ 0.1;  % initial solution 
dp= 0.001;
while abs(pstar-pstar_new)>0.001
pstar = pstar_new;
f1=F(pstar*10^5,p1,rho1)+F(pstar*10^5,p2,rho2);    % compute The derivative of F at pstar 
f2=F((pstar+dp)*10^5,p1,rho1)+F((pstar+dp)*10^5,p2,rho2);
dfdp=(f2-f1)/dp;
pstar_new=pstar-f1/dfdp;                 % get the new pstar for next iteration
end
pstar = pstar_new*10^5;
end