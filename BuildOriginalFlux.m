function [Fluxp, Fluxn] = BuildOriginalFlux(U)
% BuildFulx function is used to build the spilting Flux vector F+, F- based
% on the Solution vector U
% glabal gamma if glabal variable is allowed
    
    %decouple the solution vector U
    gamma = 1.4;
    rho = U(:,1);
    rhou = U(:,2);
    rhoe = U(:,3);
    % compute lambda matrix
    u = rhou ./ rho;
    e = rhoe ./ rho - 0.5 * u .^ 2;
    p = (gamma - 1) * rho .* e;
    a =  (gamma.*p./rho).^0.5;
    H = a.*a./(gamma-1) + 0.5 * u .^ 2;
    lambda1 = u;
    lambda2 = u - a;
    lambda3 = u + a;
    % split the lambda to lambda+ and lambda-
    [lambda1p,lambda1n] = SplitViscosity(lambda1);
    [lambda2p,lambda2n] = SplitViscosity(lambda2);
    [lambda3p,lambda3n] = SplitViscosity(lambda3);
    % compute flux+ (Steger-Warming)
    Flux1 = 2*(gamma-1).*lambda1p + lambda2p + lambda3p;
    Flux2 = 2*(gamma-1).*lambda1p.*lambda1 + lambda2p.*lambda2 + lambda3p.*lambda3; 
    Flux3 = (gamma-1)*u.^2.*lambda1p +(H-u.*a).*lambda2p+(H+u.*a).*lambda3p;
    %Flux3 = (gamma-1)*lambda1p.*lambda1.^2 + 0.5*lambda2p.*lambda2.^2 + 0.5*lambda3p.*lambda3.^2 +(3-gamma)/2/(gamma-1).*(lambda2p+lambda3p).*a.^2;
    Fluxp = 1/2/gamma.*[rho.*Flux1,rho.*Flux2,rho.*Flux3];
    % compute flux- (Steger-Warming)
    Flux4 = 2*(gamma-1).*lambda1n + lambda2n + lambda3n;
    Flux5 = 2*(gamma-1).*lambda1n.*lambda1 + lambda2n.*lambda2 + lambda3n.*lambda3; 
    %Flux6 = (gamma-1)*lambda1n.*lambda1.^2 + 0.5*lambda2n.*lambda2.^2 + 0.5*lambda3n.*lambda3.^2 +(3-gamma)/2/(gamma-1).*(lambda2n+lambda3n).*gamma.*p./rho;
    Flux6 = (gamma-1)*u.^2.*lambda1n +(H-u.*a).*lambda2n+(H+u.*a).*lambda3n;
    Fluxn = 1/2/gamma.*[rho.*Flux4,rho.*Flux5,rho.*Flux6];
    
end

function [lambda_postive,lambda_negative] = SplitViscosity(lambda)
%SpiltViscosity fuction is used to spilt the eigenvalues vector by Steger-Warming method
% for In order to enable the reduction of interruptions and oscillations, 
% artificial viscosity is introduced
epsilon = 1*10^-30;            % artificial viscosity can be select from 0.000001-0.001
lambda_postive = 0.5.*((lambda.^2+epsilon^2).^0.5+lambda);
lambda_negative = 0.5.*(-(lambda.^2+epsilon^2).^0.5+lambda);
end