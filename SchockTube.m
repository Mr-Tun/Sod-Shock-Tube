clear; close all

%% Define the parameter of the test gas (air)
% gobal gamma Rg Cv
gamma=1.4;                      % the adiabatic index γ for air is approximately 1.4
Rg=287.05;                      % the specific gas constant for air --J/(kg·K)
Cv = 1/(gamma-1)*Rg;            % compute the specific heat capacity at constant volume  --J/(kg·K)

%% Define initial conditions of shock tube
L = 2;                              % the length of the shock tube --m
p0 = 1*10^5;                        % the initial pressure --Pa
rho0 = 1;                           % the initial density   --kg/m³
p1 = 2*p0;                          % the pressure in the x<=1
rho1 = 2*rho0;                      % the density in the x<=1
T1 = GasState(p1,rho1) ;            % the temperature in the x<=1
p2 = p0;                            % the pressure in the 1<x<=2
rho2 = rho0;                        % the density in the 1<x<=2
T2 = GasState(p2,rho2);             % the temperature in 1<x<=2

%% Grid generation and initializing
N = 100;                                        % the number of the grid points
dx = L/N;                                       % the length of the cell (deltax)
p = [p1 * ones(N/2,1); p2*ones(N/2,1)];         % pressure distributions when t=0
rho = [rho1 * ones(N/2,1); rho2*ones(N/2,1)];   % density distributions when t=0
T = [T1* ones(N/2,1); T2*ones(N/2,1)];          % temperature distributions when t=0
u = zeros(N,1);                                 % velocity distributions when t=0
U = [rho, rho.*u, rho.*(Cv.*T+0.5.*u.^2)];      % initialize the solution vector
Uhis(:,:,1) = U;                                % Initialize Uhis x-axis is the space, z-axis is time, and y-axis is U1,U2,U3
t = 0;                                          % initial time --s
timestep = 0;                                   % initial timestep
time = [];                                      % initial history vector of time
time(1,1)=t;
space = 0:L/(N-1):L;                            % initial schock tube space
space = space';

%% Start integral process
while timestep<1000&&t<0.002                    % set end conditions
    u = U(:,2)./U(:,1);
    e = U(:,3)./U(:,1) - 0.5*u.^2;
    T = e./Cv;
    dt = 0.1*dx/max(abs(u) + (gamma*Rg*T).^0.5); % CFL condition
    t = t+dt;                                       % current time in shock tube
    timestep = timestep+1;                          % current timestep
    U = Uhis(:,:,timestep);                         % previous solution vector
    [Fluxp, Fluxn] = BuildOriginalFlux(U);                  % compute Flux+ and Flux-
    CFp = ones(length(Fluxp),3);CFn = ones(length(Fluxn),3);
    for i=1:length(Fluxp)-1
    CFp(i,:) = BuildCharacteristicFlux(U(i,:),U(i+1,:),Fluxp(i,:));
    CFn(i,:) = BuildCharacteristicFlux(U(i,:),U(i+1,:),Fluxn(i,:));
    end
    for j = 4:1:N-3
        %% A. Original Flux (Conservative scheme)
        
        
        % fpp=Fluxp(j,:); fpn=Fluxp(j-1,:);fnn=Fluxn(j,:);fnp=Fluxn(j+1,:);   
        
        fpp = VanLeer(Fluxp(j-1,:),Fluxp(j,:),Fluxp(j+1,:));
        fpn = VanLeer(Fluxp(j-2,:),Fluxp(j-1,:),Fluxp(j,:));
        fnp = VanLeer(Fluxn(j+2,:),Fluxn(j+1,:),Fluxn(j,:));
        fnn = VanLeer(Fluxn(j+1,:),Fluxn(j,:),Fluxn(j-1,:));
        
        
        %% B. Characteristic Flux
        % cfpp = MinMod(CFp(j-1,:),CFp(j,:),CFp(j+1,:));
        % cfpn = MinMod(CFp(j-2,:),CFp(j-1,:),CFp(j,:));
        % cfnp = MinMod(CFn(j+2,:),CFn(j+1,:),CFn(j,:));
        % cfnn = MinMod(CFn(j+1,:),CFn(j,:),CFn(j-1,:));
        % fpp = ReconstructFlux(U(j,:),U(j+1,:),cfpp);
        % fnp = ReconstructFlux(U(j,:),U(j+1,:),cfnp);       
        % fpn = ReconstructFlux(U(j,:),U(j-1,:),cfpn);
        % fnn = ReconstructFlux(U(j,:),U(j-1,:),cfnn);
        
        %% Commented out code A when using B, Commented out code B when using A
        phij = ((fpp+fnp) - (fpn+fnn))/dx;
        U(j,:) = U(j,:) - phij*dt;                       % 1st order in time
    end

    Uhis(:,:,timestep+1) = U;
    time(timestep,1)= t;
end

%% Decouple the solution vector
rhohis(:,:) = Uhis(:,1,:);
rhouhis(:,:) = Uhis(:,2,:);
rhoehis(:,:) = Uhis(:,3,:);
uhis = rhouhis./rhohis;
ehis = rhoehis ./ rhohis - 0.5 * uhis .^ 2;
This = ehis./Cv;
phis = (gamma - 1) * rhohis .* ehis;

%% Visualization of results
figure(1)                          %Pressure Distribution on the Space at Different Time
subplot(2,2,1)
epoch = round(timestep*0.5);                   %select the epoch time
[space1,Uexact] = Analytic(time(epoch));
legend_str =[num2str(time(epoch)*1000),'\fontname{Times New Roman}ms'];
plot(space,phis(:,epoch)/1000000,'LineWidth', 2)
xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
ylabel('\fontname{Times New Roman}\itp/\rmMPa','FontSize', 16);
hold on
plot (space1, Uexact(:,3)/1000000,'LineStyle','--','LineWidth',2);
legend('Numerical Solution','Analytic(Exact) Solution','location','best');
title('\fontname{Times New Roman}Pressure Distribution in Shock Tube at',legend_str);

subplot(2,2,2)                           % Density Distribution on the Space at Different Time
    plot(space,rhohis(:,epoch),'LineWidth', 2)
    xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
    ylabel('\fontname{Times New Roman}\itρ/\rmkg/m³','FontSize', 16);
    hold on
plot (space1, Uexact(:,1),'LineStyle','--','LineWidth',2);
legend('Numerical Solution','Analytic(Exact) Solution','location','best');
title('\fontname{Times New Roman} Density Distribution in Shock Tube at',legend_str);

subplot(2,2,3)                           % Velocity Distribution on the Space at Different Time
    plot(space,uhis(:,epoch),'LineWidth', 2)
    xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
    ylabel('\fontname{Times New Roman}\itu/\rmm/s','FontSize', 16);
    hold on
    plot (space1, Uexact(:,2),'LineStyle','--','LineWidth',2);
legend('Numerical Solution','Analytic(Exact) Solution','location','best');
title('\fontname{Times New Roman} Velocity Distribution in Shock Tube at',legend_str);