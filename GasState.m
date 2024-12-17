function [T] = GasState(p,rho)
%GasState function uses the ideal gas equation of state to solve the air's temperature by
%pressure and density 
% global Rg
Rg = 287.05;
T = p/rho/Rg;  % temperature of the gas --k
end

