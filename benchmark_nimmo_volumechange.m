% Script to solve coupled ice shell thermal and stress evolution
% This version of the code is meant to reproduce Nimmo (2004) - Volume
% change effect. The basic idea is that the ice shell begins as 2.4 km
% thick and with a temperature evolution given by the stefan solution.
% I was basically unable to reproduce Nimmo's results
%
%
% Max Rudolph, March 19, 2020
clear;
close all;

% Numerical parameters
nr = 101; % number of grid points
relaxation_parameter=.05; % used in nonlinear loop.
maxiter=300;
% Define physical constants and parameters
% Physical constants
seconds_in_year = 3.1558e7;
R=8.314e-3;     % in kJ/mol/K
% Boundary conditions and internal heating
H=0; % internal heating rate.
Tb=270;
Ts=100;
Ro = 1500*1000;             % outer radius of ice shell (m)
initial_thickness = 2400;
Ri = Ro-initial_thickness;         % inner radius of ice shell (m)
Rc = 1.60e5;             % core radius (m)
% Elastic and Viscous properties
E = 5e9;        % shear modulus of ice (Pa)
nu = 0.3;      % Poisson ratio of ice (-)
beta_w = 4e-10; % Compressibility of water (1/Pa)
alpha_l = 1e-4*0; % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
rho_i=900;      % density of ice (kg/m^3)
rho_w=1000;     % density of water (kg/m^3)
Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=1e15;       % basal viscosity (Pa-s)
mu = @(T) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
% Thermal properties
Cp = 2100; %heat capacity of ice J/kg/K
Lf = 334*1000; % latent heat of fusion (J/kg)
kappa = 1e-6;% m/s/s
k=kappa*rho_i*Cp;

% Benchmark-specific things
% temperature solution
lam1 = 0.65;
cooling_age = @(zm) (zm/2/lam1)^2/kappa;
initial_cooling_age = cooling_age(initial_thickness);
zm = @(t) 2*lam1*sqrt(kappa*t);
dTdt = @(z,t) -(Tb-Ts)/erf(lam1)*exp(-z.^2/4/kappa/t).*z/(2*sqrt(pi*kappa*t^3));
stefan_T = @(z,t)  erf( (z)/2./sqrt(kappa*t) )/erf(lam1)*(Tb-Ts)+Ts

% calculate maxwell time at 100, 270
fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(100)/E,mu(Tb)/E);
fprintf('Thermal diffusion timescale %.2e\n',(4e4)^2/kappa);
% set end time and grid resolution
t_end = 90e6*seconds_in_year;
% dt = 1e4*seconds_in_year; % time step in seconds
dtmax = 1e4*seconds_in_year;
dtmin = seconds_in_year;
% dt1 = 3600; % size of first timestep
% times = logspace(log10(dt1),log10(t_end+dt1),1e4)-dt1;
plot_interval = t_end;
save_interval = 1e5*seconds_in_year;
save_depths = [0 1 2 3 5]*1000;
nsave = ceil(t_end/save_interval) + 1; nsave_depths = length(save_depths);
sigma_t_store = zeros(nsave_depths,nsave);

% set up the grid
grid_r = linspace(Ri,Ro,nr); % set up the grid

% initialize solution vectors (IC)
sigma_r_last = zeros(nr,1); % initial stresses
sigma_t_last = zeros(nr,1); % initial stresses
T_last = linspace(Tb,Ts,nr)';
% stefan initial condition
eta = (Ro-grid_r)/2/sqrt(kappa*initial_cooling_age);
theta = erf(eta)/erf(lam1);
T_last(:) = theta*(Tb-Ts)+Ts;
% initialize with a stefan solution corresponding to initial thickness
er_last = zeros(nr,1); % strains
et_last = zeros(nr,1);
ur_last = zeros(nr,1);      % displacement
z_last = 0; % total amount of thickening
Pex_last = 0; %initial overpressure

% Set up plot
hf2=figure();

plot_times = [0.0 0.8 1.6 3.2 6.4 12 24 46 90]*1e6*seconds_in_year; iplot=2;
hf=figure(2);
subplot(1,4,1); % sigma_r and sigma_t
h=plot(sigma_r_last,Ro-grid_r); hold on;
plot(sigma_r_last,Ro-grid_r,'--','Color',h.Color);
% h=legend('\sigma_r','\sigma_t','Interpreter','tex'); h.AutoUpdate=false;
title('Stress (Pa)','Interpreter','tex');
ylabel('r (m)');
set(gca,'YDir','reverse');
subplot(1,4,2); % e_r and e_t
h=plot( sigma_r_last,Ro-grid_r); hold on;
plot( sigma_r_last,Ro-grid_r,'--','Color',h.Color); hold on;
% h=legend('r','t'); h.AutoUpdate=false;
title('Strain (-)','Interpreter','tex');
set(gca,'YDir','reverse');
subplot(1,4,3); % temperature
plot(T_last,Ro-grid_r); hold on; title('T (K)','Interpreter','tex'); set(gca,'YDir','reverse');
subplot(1,4,4); % radial displacement (u)
plot(ur_last,Ro-grid_r); hold on; title('u_r');
set(gca,'YDir','reverse');
last_plot_time = 0;

fig1a.h = figure(); % Nimmo's Figure 1a
subplot(2,1,1);
[ax,h1,h2]=plotyy((Ro-grid_r)/1e3,sigma_t_last/1e6,(Ro-grid_r)/1e3,T_last);
h1.DisplayName = '0.0';
fig1a.ax = ax;
axis(ax(1));
h2.Color = h1.Color;
h2.LineStyle = '--';
hold(ax(1)); hold(ax(2));
set(ax,'Xlim',[0 20]);
set(ax(1),'YLim',[-5 40]);
set(ax(1),'YTick',[-5:5:40]);
set(ax(2),'YTick',[100:20:180]);
set(ax(1),'YTickLabelMode','auto');
ylabel(ax(1),'Tangential Stress (MPa)');
xlabel(ax(1),'Depth (km)');
ylabel(ax(2),'Temperature (K)');
set(ax(2),'YLim',[100 180]);

time=0; itime=1;
% save initial state
isave = 1;
sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
time_store(isave) = time;
last_store = time; isave = isave+1;


while time < t_end
    % In each timestep, we do the following
    % 1. Calculate the amount of basal freeze-on and advance the mesh
    % 2. Solve the heat equation using an implicit method
    % 3. Solve for sigma_r
    % 4. Calculate sigma_t
    % 5. Calculate the radial displacements u(r)
    
    % 1. form discrete operators and solve the heat equation
    %     L = zeros(nr,nr);
    %     R = zeros(nr,1);
    %     for i=1:nr
    %         r = grid_r(i);
    %         if i==1
    %             drm = grid_r(i+1)-grid_r(i);
    %         else
    %             drm = grid_r(i)-grid_r(i-1);
    %         end
    %         if i==nr
    %             drp = drm;
    %         else
    %             drp = grid_r(i+1)-grid_r(i);
    %         end
    %         rA = r + drp/2;
    %         rB = r - drm/2;
    %         kA = k;% thermal conductivities
    %         kB = k;
    %         dr = rA-rB;
    %         coef_plus = -kA*rA^2/r^2/drp/dr;
    %         coef_center = rho_i*Cp/dt + kA*rA^2/r^2/drp/dr + kB*rB^2/r^2/drm/dr;
    %         coef_minus = -kB*rB^2/r^2/drm/dr;
    %         R(i) = rho_i*Cp/dt*T_last(i) + H;
    %         L(i,i) =  coef_center;
    %         if( i==1 )
    %             L(i,i+1) = coef_plus-coef_minus;
    %             R(i) = R(i) - 2*Tb*coef_minus;
    %         elseif i==nr
    %             L(i,i-1) = coef_minus-coef_plus;
    %             R(i) = R(i) - 2*Ts*coef_plus;
    %         else
    %             L(i,i-1) = coef_minus;
    %             L(i,i+1) = coef_plus;
    %         end
    %     end
    %     T = L\R;
    
    % 2. Calculate thickening
    % basal heat flux
    %     qb = -k*(T(2)-T(1))/(grid_r(2)-grid_r(1)); % W/m^2
    % thickening would be dx/dt = qb/(L*rho_i)
    %     delta_rb = dt*qb/Lf/rho_i;
    %     z = z_last + delta_rb;
    
    % calculate the timestep
    delta = Inf;
    dt=dtmax*2;
    while abs(delta) > (grid_r(2)-grid_r(1))/2 && dt >= dtmin
        dt = dt/2;
        new_thickness = zm(initial_cooling_age+time+dt);
        delta = new_thickness - (Ro-grid_r(1));
    end
    fprintf('Timestep %e elapsed time %e\n',dt,time+dt);
    z=new_thickness-initial_thickness;
    delta_rb = new_thickness-initial_thickness;
    f = (1-2*initial_thickness/Ro);
    delta_rho = (rho_w-rho_i);
    surface_uplift = (delta_rb)*(delta_rho)/rho_w*f/(1+(f*(delta_rho)/rho_w));
    
    % calculate new ocean pressure (Manga and Wang 2007, equation 5)
    %     Pex_pred = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur_last(1));% ur_last because we don't yet know the uplift
    Pex_pred = surface_uplift*E/(1-nu)*2/Ro^2*(new_thickness/4); % last term assumes elastic thickness ~1/4 total thickness
    % re-mesh onto new grid
    new_grid_r = linspace(Ri-delta_rb,Ro,nr);
    [T_last,sigma_r_last,sigma_t_last,er_last,et_last] = interpolate_solution(new_grid_r,grid_r,T_last,sigma_r_last,sigma_t_last,er_last,et_last,Tb);
    grid_r = new_grid_r; % end interpolation step
    
    % calculate d/dr(Tdot)
    %     Tdot = (T-T_last)/dt;
%     Tdot = zeros(size(T_last));
%     Tdot(:) = dTdt(Ro-new_grid_r,initial_cooling_age+time+dt);
%     T = T_last + Tdot*dt;
    T = stefan_T( Ro-new_grid_r',initial_cooling_age+time+dt);
    Tdot = (T-T_last)/dt;
    
    dTdotdr = zeros(nr,1);
    for i=2:nr-1
        dTdotdr(i) = (Tdot(i+1)-Tdot(i))/(grid_r(i+1)-grid_r(i));
    end
    dTdotdr(1) = (Tdot(2)-Tdot(1))/(grid_r(i+1)-grid_r(i));
    dTdotdr(nr) = (Tdot(nr)-Tdot(nr-1))/(grid_r(nr)-grid_r(nr-1));
    
    % 3. Nonlinear loop over pressure.
    % because the ocean pressure depends on the uplift, we make a guess
    % (above). Using this guess, we calculate stresses, strains, and
    % displacements. Then we re-calculate the pressure using the new value
    % of radial displacement. We continue until the pressure used in the
    % calculations has converged to the pressure consistent with the
    % calculated displacement;
    for iter=1:maxiter
        if iter>1
            %             Pex = Pex + 0.01*(Pex_post-Pex);
            durdP = (1-nu)/E*Ro^2/2*1/new_thickness;
            dP = (surface_uplift-ur(end))/durdP;
            Pex = Pex + 0.1*dP;
            
        else
            Pex = Pex_last;
        end
        % 3a. Assemble and solve equations for radial stress
        % calculate viscosity at each node (needed by solver)
        mu_node = zeros(nr,1);
        for i=1:nr
            mu_node(i) = mu(T(i));
        end
        [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-Pex,E,nu,dt);
        
        dT = T-T_last;
        
        dsigma_t = sigma_t - sigma_t_last;
        dsigma_r = sigma_r - sigma_r_last;
        
        de_t = 1/E*(dsigma_t-nu*(dsigma_t+dsigma_r))+alpha_l*dT + dt/2*(sigma_tD./mu_node); % change in tangential strain
        de_r = 1/E*(dsigma_r-2*nu*dsigma_t)         +alpha_l*dT + dt/2*(sigma_rD./mu_node); % HS91 equations A5-6
        er = er_last + de_r;
        et = et_last + de_t;
        ur = grid_r'.*et; %radial displacement
        % re-calculate excess pressure using new uplift
        % Pex_post = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur(1));% ur_last because we don't yet know the uplift
%         fprintf('end of nonlinear loop. actual uplift %.2e desired uplift %.2e\n',ur(end),surface_uplift);
        
        
        % check for convergence
        if abs( ur(end)-surface_uplift )/abs(surface_uplift) < 1e-6
            fprintf('time %.2e: converged in %d iterations. actual uplift %.2e desired uplift %.2e\n',time/seconds_in_year,iter,ur(end),surface_uplift);
            break;
        elseif iter==maxiter
            error('Nonlinear loop failed to converge');
        end
    end%end nonlinear loop
    
    % 6. advance to next time step and plot (if needed)
    sigma_r_last = sigma_r;
    sigma_t_last = sigma_t;
    T_last = T;
    er_last = er;
    et_last = et;
    z_last = z;
    Pex_last = Pex;
    
    time = time + dt;
    
    if (time >= plot_times(iplot) || time >= t_end)
        iplot = iplot+1;
        
        figure(hf);
        subplot(1,4,1);
        h=plot(sigma_r,Ro-grid_r);
        plot(sigma_t,Ro-grid_r,'--','Color',h.Color);
        subplot(1,4,2);
        h=plot(er,Ro-grid_r);
        plot(et,Ro-grid_r,'--','Color',h.Color);
        subplot(1,4,3);
        plot(T,Ro-grid_r);
        subplot(1,4,4);
        plot(ur,Ro-grid_r);
        
        
        figure(hf2);
        plot(ur(end),sigma_t(end),'.'); hold on;
        
        figure(fig1a.h); % Nimmo's Figure 1a
%         str1 = sprintf('%.1f',plot_times(iplot)/seconds_in_year/1e6);
        h=plot(fig1a.ax(1),(Ro-grid_r)/1e3,sigma_t_last/1e6);
%         fig1a.lh(end+1) = h;
        plot(fig1a.ax(2),(Ro-grid_r)/1e3,T_last,'--','Color',h.Color);
        
        
        last_plot_time = time;
        drawnow();
    end
    if (time-last_store >= save_interval || time >= t_end)
        sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
        time_store(isave) = time;
        last_store = time; isave = isave+1;
    end
    
    
end
%% add legends
for i=1:length(plot_times)
    labels{i} = sprintf('%.1f',plot_times(i)/seconds_in_year/1e6);
end
figure( fig1a.h );
axis(fig1a.ax(1));
legend(labels,'Location','southeast');


%% Nimmo's figure 1b
figure(fig1a.h);
labels = {};
for i=1:length(save_depths)
    labels{i} = sprintf('%.1f km',save_depths(i)/1000);
end
subplot(2,1,2);
plot(time_store(1:isave-1)/seconds_in_year/1e6,sigma_t_store(:,1:isave-1)/1e6);
legend(labels,'Location','southeast');
xlabel('Time (Myr)');
ylabel('Tangential Stress (MPa)');
set(gcf,'Color','w');
h=gcf;
h.Position(3:4) = [390 580];
saveas(gcf,'Nimmo_Figure2.eps','psc2')


