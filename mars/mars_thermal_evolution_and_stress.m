function results = mars_thermal_evolution_and_stress(parameters)
% function to solve coupled ice shell thermal and stress evolution
% Max Rudolph, March 19, 2020
% adapted for other moons by Alyssa Rhoden, 6/2021
%
% Cite:
% Rudolph, M.L., Manga, M., Walker, M., and Rhoden, A. Cooling Crusts
% Create Concommitant Cryovolcanic Cracks. Geophysical Research Letters
% 49(5) e2021GL094421

% clear;
% close all;
% addpath core; % this is where the helper functions live.
% addpath mars;
% addpath ~/sw/matlab/crameri

% parameters to vary across models
% crustal thickness
% crustal heat production 0.3-0.7
% length of initial stress-free time
% rheology (arrhenius vs. power law)
% rheology (reference viscosity)

% compute
% depth of maximum differential stress
% depth of transition from compression to tension
% for crust_heat_fraction = [0.3 0.5 0.7]


nrs = [512]; % number of points used in the radial direction

% Settings related to numerics
label='Mars';
seconds_in_year = 3.1558e7;
max_depth = 4e5; % maximum depth for saving solution values (m)
relaxation_parameter = 1e-3;%1e-3; % used for fixed point iteration in pressure convergence loop.
t_end = 4500e6*seconds_in_year;%  3*perturbation_period; 5e8*seconds_in_year;
dtmax = 5e6*seconds_in_year;
dtmin = 100*seconds_in_year;%*seconds_in_year;
% no_stress_time = 1.5e9*seconds_in_year; % time before which stresses are not allowed to increase
no_stress_time = parameters.no_stress_time*seconds_in_year;

% Stuff related to the Mars thermal evolution model
arh =2.0;   % constant from Michaut equation 12
C   =0.5;   % Davaille and Jaupart 1993 constant for heat flux

% Rheology
viscosity_model = 2;    % 0 = Nimmo (2004), 1 = Goldsby and Kohlstedt (2001), 2=Arrhenius
% viscosity.d = NaN;      %1e-3; % grain size in m used to calculate the viscosity (for G-K)
% viscosity.P = NaN;      %1e5; % Pressure in MPa used to calculate the viscosity (for G-K)
% mub=3e20;               % Reference viscosity (at reference temperature)
mub = parameters.viscosity;
Tref = 1600;            % Reference temperature, Kelvin.
Q=300;                  % value from Michaut et al. 2025, kJ/mol
R=8.314e-3;             % in kJ/mol/K
if viscosity_model == 2
    mub=3e20;               % Reference viscosity (at reference temperature)
    Q=300;                  % value from Michaut et al. 2025, kJ/mol
    mu = @(T,P,stress) mub*exp(Q/R*(1./T - 1./Tref)); % Michaut et al. 2025 - Arrhenius form
    dTnu = @(T) R/Q*T^2; % rheological temperature scale (positive sign??)
elseif viscosity_model==3
    viscosity.d = 7.08e-3;  %grain size in m. 7.08e-3 gives 3e20 Pa-s at 1 GPa pressure and 0 stress
    % viscosity.P = NaN;      %1e5; % Pressure in MPa used to calculate the viscosity (for G-K)
    Q=375; % Q used in hirth and kohlstedt model - use it for the mantle too?
    mu = @(T,P,stress) hirth_kohlstedt(stress,T,viscosity.d,P);
    dTnu = @(T) R/Q*T^2;
end

% crust properties
rhoc=2900;
kcrust=3.0;

% Mechanical properties
nu = 0.25;              % Poisson ratio of lithosphere (-)
E = 0.8e11;             % shear modulus of lithosphere (Pa) (T&S Appendix B5, for basalt/gabbro)
K_eff = 4e11;           % effective bulk modulus of mantle+core (Pa)
alpha_v = 2.5e-5;       % volumetric thermal expansivity (1/K)
alpha_l = alpha_v/3;    % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
alpha_v_bl = alpha_v;  % define the alpha_v for the lid evolution separatey so that we can isolate it later.
% alpha_v = 0; % eliminate mantle shrinkage...

% Heat transport properties:
Cp = 1150;            % specific heat capacity, J/kg/K
k = @(T) 4;           % Thermal conductivity, W/m/K (Thieriet et al., 2019)

% Initial and boundary conditions
thickness = 70e3;       % Initial lithosphere thickness (m)
Ts = 250;               % Surface temperature, K.
tstart = 0*seconds_in_year; % starting time of the model, used for radioactive heating.
% note - initial temperature profile is steady state, calculated later
% Initial basal temperature
Tm0 = 1700;             % initial mantle temperature
DTbl = arh*dTnu(Tm0);
Tb = Tm0-DTbl;          % Temperature at base of lid

% Planet properties
Ro = 3.3895e6;      % outer radius of lithosphere
Ri = Ro-thickness;  % initial inner radius of lithosphere
Rc = 1.830e6;       % core radius, m (Samuel et al., 2023)
h_crust = 60e3;     % crust thickness (assumed constant)
% crust_heat_fraction = 0.5; % fraction of primitive mantle heat production concentrated within crust
crust_heat_fraction =parameters.crust_heat_fraction;
mars_mass = 6.4169e23;
silicate_mass = 0.75*mars_mass; % assuming core 25% as in Khan et al. 2022 EPSL
silicate_density = silicate_mass / (4/3*pi*(Ro^3-Rc^3));% density of bulk silicate mars
% rho = silicate_density; % uniform density approximation
crust_mass = rhoc*4/3*pi*(Ro^3-(Ro-h_crust)^3);     % mass of the crust
mantle_mass = silicate_mass - crust_mass; % mass of the mantle
mantle_density = mantle_mass/( 4/3*pi*((Ro-h_crust)^3-Rc^3));
rho=mantle_density;
crust_mass_fraction = crust_mass/mantle_mass;
% compute the heating per unit mass in the crust
% mantle heating = [h]*rho*V
crustal_heating_factor = crust_heat_fraction/crust_mass_fraction; % this is the enrichment in volumetric heating relative to primitive mantle material
mantle_heating_factor = (1-crust_heat_fraction)/(1-crust_mass_fraction);
g = 3.73;           % surface gravity (m/s^2)

kappa = k(Tb)/rho/Cp;           % thermal diffusivity m^2/s

% make a structure to store all of the parameters defined up to this point
vars = whos();
for i=1:length(vars)
    state.(vars(i).name) = eval(vars(i).name);
end
results.state=state;

% Failure criterion:
tensile_strength = 1e99; % tensile strength, Pa
cohesion = 1e99;  % plastic yield strength, MPa
friction = 0.0; % friction angle for plastic yielding


if viscosity_model == 0
    label = [label '-nimmovisc'];
elseif viscosity_model == 1
    label = [label '-goldsbykohlstedt'];
elseif viscosity_model == 2
    label = [label '-arrhenius']
else
    error('not implemented');
end

for inr=1:length(nrs) % loop over nr values for resolution tests
    ifail = 1; % index into list of times at which failure occurred.
    nr = nrs(inr); % number of grid points
    maxiter=1000;
    time=0;

    % calculate maxwell time at Ts, Tb.
    fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(Ts,0)/E,mu(Tb,0)/E);
    fprintf('Thermal diffusion timescale %.2e\n',(Ro-Ri)^2/kappa);

    plot_interval = 5e6*seconds_in_year;
    save_interval = 1e6*seconds_in_year;
    save_depths = linspace(0,max_depth,500);

    nsave = ceil(t_end/save_interval) + 1;
    nsave_depths = length(save_depths);
    sigma_t_store = zeros(nsave_depths,nsave);

    results.time = NaN*zeros(nsave,1); results.time(1) = 0;
    results.thickness = zeros(nsave,1); results.thickness(1) = Ro-Ri;
    results.z = zeros(nsave,1);
    results.Ri = zeros(nsave,1); results.Ri(1) = Ri;
    results.qb = zeros(nsave,1);
    results.sigma_t = NaN*zeros(nsave_depths,nsave);
    results.sigma_r = NaN*zeros(nsave_depths,nsave);
    results.e_t = NaN*zeros(nsave_depths,nsave);
    results.e_r = NaN*zeros(nsave_depths,nsave);
    results.Pex = zeros(nsave,1);
    results.Tm = zeros(nsave,1); results.Tm(1) = Tm0;
    results.Pex_crit = zeros(nsave,1);
    results.dTdr = zeros(nsave_depths,nsave);
    results.T = zeros(nsave_depths,nsave);
    results.Tb = zeros(nsave,1);
    results.ur = zeros(nsave_depths,nsave);
    results.ur_base = NaN*zeros(1,nsave);
    results.failure_time = zeros(1,nsave);
    results.failure_P = zeros(1,nsave);
    results.failure_Pex_crit = zeros(1,nsave);
    results.failure_dP = zeros(1,nsave);
    results.failure_thickness = zeros(1,nsave);
    results.failure_top = zeros(1,nsave);
    results.failure_bottom = zeros(1,nsave);
    results.failure_erupted_volume = NaN*zeros(1,nsave);
    results.failure_erupted_volume_pressurechange = NaN*zeros(1,nsave);
    results.failure_erupted_volume_volumechange = NaN*zeros(1,nsave);
    results.stresss_crossover_depth = NaN*zeros(1,nsave);
    results.maximum_differential_stress = NaN*zeros(1,nsave);
    results.maximum_stress_depth = NaN*zeros(1,nsave);
    results.maximum_stress_depth = NaN*zeros(1,nsave);

    erupted_volume = 0;
    erupted_volume_pressurechange = 0;
    erupted_volume_volumechange = 0;

    % set up the grid
    grid_r = linspace(Ri,Ro,nr); % set up the grid

    % initialize solution vectors (IC)
    % eccentricity_last = e0;
    sigma_r_last = zeros(nr,1); % initial stresses
    sigma_t_last = zeros(nr,1); % initial stresses
    siiD_last = zeros(nr,1); % deviatoric stress invariant - used for viscosity
    % T_last = zeros(nr,1);
    % Temperature initial condition
    % option 1: Initialize T with stefan solution
    % T_last(:) = solve_stefan_analytic(grid_r(end)-grid_r,k(Tb),rho_lith,Cp,Lf,Tb,Ts);
    % option 2: linear temperature profile
    % T_last = 1600*ones(nr,1);% add temperature solution here!
    % T_last(:) = linspace(Tb,Ts,nr);
    % option 3 - solve the temperature equation to obtain steady solution for current heating:
    T_last = linspace(Tb,Ts,nr);
    H = zeros(nr,1);
    iscrust = grid_r>=(Ro-h_crust);
    H( iscrust ) = mars_heating((time+tstart)/seconds_in_year)*crustal_heating_factor*rhoc;
    H(~iscrust ) = mars_heating((time+tstart)/seconds_in_year)*mantle_heating_factor*rho;
    kvec = k(Tm0)*ones(nr,1);
    kvec(iscrust) = kcrust;
    Cpvec = Cp*ones(nr,1);
    rhovec = rho*ones(nr,1); rhovec(iscrust) = rhoc;

    [T_last,dTdotdr] = solve_temperature_shell_mars(grid_r,T_last,Tb,Ts,kvec,rhovec,Cpvec,H,Inf,0.0); % call solver with infinite timestep
    Tm = Tm0; % mantle temperature

    er_last = zeros(nr,1); % strains
    et_last = zeros(nr,1);
    ur_last = zeros(nr,1); % displacement
    z_last = 0;    % total amount of thickening
    dzdt_last = 0; % thickening rate
    Pex_last = 0; %initial overpressure

    % Set up plot
    if parameters.do_plots
        hf2=figure();

        plot_times = linspace(0,t_end,5); iplot=2;
        hf=figure();
        subplot(1,4,1); % sigma_r and sigma_t
        h=plot(sigma_r_last,Ro-grid_r); hold on;
        plot(sigma_t_last,Ro-grid_r,'--','Color',h.Color);
        % h=legend('\sigma_r','\sigma_t','Interpreter','tex'); h.AutoUpdate=false;
        title('Stress (Pa)','Interpreter','tex');
        ylabel('Depth (m)');
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
        fig1a.ax = ax;
        h2.Color = h1.Color;
        h2.LineStyle = '--';
        hold(ax(1)); hold(ax(2));
        % set(ax,'Xlim',[0 10]);
        % set(ax(1),'YLim',[-10 40]);
        % set(ax(1),'YTick',[-10:5:40]);
        % set(ax(2),'YTick',[100:20:180]);
        set(ax(1),'YTickLabelMode','auto');
        ylabel(ax(1),'Tangential Stress (MPa)');
        xlabel(ax(1),'Depth (km)');
        ylabel(ax(2),'Temperature (K)');
        % set(ax(2),'YLim',[100 180]);
    end
    itime=1;
    % save initial state
    isave = 1;
    sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
    time_store(isave) = time;
    last_store = time; isave = isave+1;

    failure_mask = false(size(grid_r)); % stores whether failure occurred
    failure_time = zeros(size(grid_r)); % stores the time at which failure occurred

    while time < t_end && (Ri-z_last > Rc)
        % In each timestep, we do the following
        % 1. Calculate the amount of LID THICKENING and advance the mesh
        % 2. Solve the heat equation using an implicit method
        % 3. Solve for sigma_r
        % 4. Calculate sigma_t
        % 5. Calculate the radial displacements u(r)

        % 1. Calculate LID THICKENING and interpolate old solution onto new mesh
        % calculate heat flux
        dt = dtmax;
        Tg = Tb-(T_last(2)-Tb);
        dTdr_b_last = (T_last(2)-Tg)/2/(grid_r(2)-grid_r(1));
        qlid = -k(Tb)*dTdr_b_last; % this is the upward conducted heat flow from the last timestep

        % [tidal_heating,total_heating] = Qbelow(grid_r(end)-grid_r(1),eccentricity_last);
        % total_heating = 0; % for now, to obtain a solution.
        % qb_net = qb - total_heating; % first term is conducted heat. second term is heat supplied from below.

        % Implement the thermal evolution model...
        D = Ro-Ri-z_last; % z is the amount by which the lid has thickened
        mantle_volume = 4/3*pi*((Ro-D)^3-Rc^3);
        Cm = rho*Cp*mantle_volume; % mantle heat capacity
        Slid = 4*pi*(Ro-D)^2;
        h_conv = mantle_heating_factor*rho*mars_heating(time/seconds_in_year); % mantle volumetric heat production
        % boundary layer heat transport into the lid:
        qbl = C*k(Tm)*(alpha_v_bl*rho*g/kappa/mu(Tm,0))^(1/3)*dTnu(Tm)^(4/3);
        % temperature difference across the boundary layer:
        DTbl = arh*dTnu(Tm);
        Tl = Tm-DTbl;% temp at base of conductive layer
        dH = rho*Cp*DTbl; % enthalpy change across the lid
        % rate of change of mantle temperature
        dTmdt = 1/Cm * (-Slid*qlid + h_conv*mantle_volume); % Michaut et al. Equation 18
        % rate of change of lid thickness.
        dDdt = 1/dH * (qlid-qbl);   % Michaut et al. Equation 19
        dTldt = (Tm-DTbl - T_last(1))/dt; %rate of change of temperature at base of lid.

        % determine the timestep - apply a courant type condition
        % to lid thickness change
        if abs(dDdt*dt) > (grid_r(2)-grid_r(1))/2
            dt = abs( (grid_r(2)-grid_r(1))/2/(dDdt) );
        end
        % apply a limiter based on mantle temperature change
        if abs(dTmdt*dt) > 1.0
            dt = abs(1.0/dTmdt);
        end
        if dt < dtmin
            dt = dtmin;
            warning('Setting dt = dtmin');
        end
        if any(failure_mask)
            dt = dtmin;
        end

        % update the mantle temperature
        delta_rb = dDdt*dt;
        z = z_last + delta_rb;
        % update the basal temperature
        Tb = Tb + dTldt*dt;
        Tm = Tm + dTmdt*dt;

        % compute the melting temperature for the new NH3 content at
        % the ocean-ice interface:
        % Tmelt = ammonia_melting(X);
        % Tb = Tmelt;
        % dzdt = delta_rb/dt;

        if (Ri-z-delta_rb <= Rc)
            % code seems to get very unstable when the ocean is too
            % thin...
            break
        end
        % calculate new ocean pressure (Manga and Wang 2007, equation 5)
        % Pex_pred = 3*K_eff*(Ri-z)^2/((Ri-z)^3-Rc^3)*( -ur_last(1) ) + K_eff*alpha_v*(Tm-Tm0); % ur_last because we don't yet know the uplift
        % Re-mesh and interpolate the solution onto the new grid.
        new_grid_r = linspace(Ri-z,Ro,nr);
        dTdr_last = (T_last(2)-T_last(1))/(grid_r(2)-grid_r(1));
        [T_last,sigma_r_last,sigma_t_last,er_last,et_last] = interpolate_solution(new_grid_r,grid_r,T_last,sigma_r_last,sigma_t_last,er_last,et_last,Tb);
        grid_r = new_grid_r; % end interpolation step

        % 2. form discrete operators and solve the heat equation
        H = zeros(nr,1);
        iscrust = grid_r>=(Ro-h_crust);
        H( iscrust ) = mars_heating((time+tstart)/seconds_in_year)*crustal_heating_factor*rhoc;
        H(~iscrust ) = mars_heating((time+tstart)/seconds_in_year)*mantle_heating_factor*rho;
        kvec = k(Tm0)*ones(nr,1);
        kvec(iscrust) = kcrust;
        Cpvec = Cp*ones(nr,1);
        rhovec = rho*ones(nr,1); rhovec(iscrust) = rhoc;

        [T,dTdotdr] = solve_temperature_shell_mars(grid_r,T_last,Tb,Ts,kvec,rhovec,Cpvec,H,dt,delta_rb);

        % 3. Nonlinear loop over pressure.
        % because the ocean pressure depends on the uplift, we make a guess
        % (above). Using this guess, we calculate stresses, strains, and
        % displacements. Then we re-calculate the pressure using the new value
        % of radial displacement. We continue until the pressure used in the
        % calculations has converged to the pressure consistent with the
        % calculated displacement;
        converged = false;
        pex_store = zeros(maxiter,1);
        pexpost_store = zeros(maxiter,1);
        for iter=1:maxiter
            if iter>100
                [tmp,ind] = unique(pex_store(1:iter-1));
                Pex = interp1(pexpost_store(ind)-pex_store(ind),pex_store(ind),0,'linear','extrap');
            elseif iter>1
                Pex = Pex + relaxation_parameter*(Pex_post-Pex);
            else
                Pex = Pex_last;
            end

            % compute height to which water would rise
            % [ptmp,Ttmp,ztmp,rhotmp] = ammonia_adiabatic_profile(X,g*rho_lith*(grid_r(end)-grid_r(1)),g);
            % rhobar = cumtrapz(ztmp,rhotmp);
            % rhobar = 1/(ztmp(end)-ztmp(1))*rhobar(end);

            % Pex_crit = (rhobar-rho_lith)*(Ro-(Ri-z))*g;

            % calculate viscosity at each node
            visc_converged = false;
            visc_iter = 100;
            ivisc = 1;
            while ~visc_converged && ivisc <= visc_iter
                % compute mu for current siiD
                if ivisc == 1
                    siiD = siiD_last;
                else
                    siiD = siiD_post;
                end


                mu_node = zeros(nr,1);
                mu_node(:) = mu(T,siiD);
                % reduce Maxwell time in region experiencing failure
                if all(failure_mask)
                    if Pex_last >= Pex_crit
                        % Calculate the volume erupted (dP)*beta*V0 + V-V0
                        pressure_contribution = (Pex_last-Pex_crit)*beta_w*(4/3*pi*((Ri-z)^3-Rc^3));
                        urelax = (Ri-z)/E*(1-2*nu)*(Pex_last-Pex_crit); % Manga and Wang (2007) equation 4
                        volume_contribution = (Ri-z)^2*urelax*4*pi; % (4*pi*R^2)*dr
                    else
                        pressure_contribution = 0;
                        volume_contribution = 0;
                    end
                    % reset stresses and uplift
                    sigma_r_last = 0*sigma_r_last;
                    sigma_t_last = 0*sigma_t_last;
                    er_last = 0*er_last;
                    et_last = 0*et_last;
                    ur_last = 0*ur_last;
                    Pex=0; % force zero pressure.
                    converged = true;
                    Ri = Ri - z;
                    z = 0;
                    z_last=0;
                    % move the inner radius effectively to the current position
                    % of the base of the ice shell. Then set the amount of
                    % freezing to zero.
                elseif( any(failure_mask) )
                    minimum_viscosity_prefactor = 0; % maximum allowable fractional reduction in viscosity
                    mu_node(failure_mask) = min(mu_node(failure_mask),max(minimum_viscosity_prefactor*mu_node(failure_mask),0.1*E*dt));  % timestep = 10x(maxwell time)
                    above_crack = find( failure_mask,1,'last');
                    above_crack_mask = false(size(failure_mask));
                    above_crack_mask( above_crack+1 : end ) = true;
                    %                 mu_node(above_crack_mask) = min( mu_node(above_crack_mask),100*E*dt ); % limit maximum viscosity to 100*maxwell time
                    %                 for i=1:3
                    %                 tmp = exp(smooth(log( mu_node )));
                    %                 mu_node(~failure_mask) = tmp(~failure_mask);
                    %                 end
                    %                 mu_node = exp(smooth(log( mu_node )));
                    if iter==1
                        Pex=0; % If failure occurs, it's better to guess that all pressure is relieved. Other choices could cause convergence problems.
                    end
                end

                % Calculate Stresses
                [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-Pex,E,nu,dt);
                siiD_post = sqrt( 0.5*(sigma_rD.^2 + 2*sigma_tD.^2) );
                norm_change = min(norm(siiD_post-siiD)/norm(siiD),norm(siiD_post-siiD));
                %                     disp([num2str(ivisc) ' change in norm of siiD:' num2str(norm_change)]);
                if isnan(norm_change)
                    keyboard
                elseif norm_change < 1e-4
                    visc_converged = true;
                end

                ivisc = ivisc+1;
            end
            if all(failure_mask)
                erupted_volume = erupted_volume + pressure_contribution + volume_contribution;
                erupted_volume_pressurechange = erupted_volume_pressurechange + pressure_contribution;
                erupted_volume_volumechange = erupted_volume_volumechange + volume_contribution;
            end

            % 5. Calculate the strains
            dT = T-T_last;
            dr1 = grid_r(2)-grid_r(1);
            dr2 = grid_r(3)-grid_r(1);
            L = [0 0 1;
                dr1^2 dr1 1;
                dr2^2 dr2 1];
            R1 = T(1:3);
            coef = L\R1;
            dTdr_b = coef(2);
            %             dTdr_b=(T(2)-Tb)/(grid_r(2)-grid_r(1));
            dT(1) = delta_rb*dTdr_b;

            dsigma_t = sigma_t - sigma_t_last;
            dsigma_r = sigma_r - sigma_r_last;
            %             mu_node(2:end-1) = exp(0.5*(log(mu_node(1:end-2))+log(mu_node(3:end))));
            de_t = 1/E*(dsigma_t-nu*(dsigma_t+dsigma_r))+alpha_l*dT + dt/2*(sigma_tD./mu_node); % change in tangential strain
            de_r = 1/E*(dsigma_r-2*nu*dsigma_t)         +alpha_l*dT + dt/2*(sigma_rD./mu_node); % HS91 equations A5-6
            er = er_last + de_r;
            et = et_last + de_t;
            ur = grid_r'.*et; %radial displacement

            ei = 2*de_t + de_r; % first invariant of strain
            de_tD = de_t - 1/3*ei;
            de_rD = de_r - 1/3*ei;
            eiiD = sqrt( 0.5*(de_rD.^2 + 2*de_tD.^2) ); % second invariant of deviatoric strain

            % re-calculate excess pressure using new uplift
            %             Pex_post = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur(1));
            Pex_post = 0*Pex_last + 3*K_eff*(Ri-z)^2/((Ri-z)^3-Rc^3)*(-(ur(1)-0*ur_last(1))) + K_eff*alpha_v*(Tm-Tm0);
            % Pex_post = Pex_last;
            % Calculate the critical excess presssure necessary to
            % erupt water onto the surface.
            %fprintf('iter %d. Pex_post %.2e Pex %.2e\n',iter,Pex_post,Pex);

            % check for convergence
            if abs( Pex_post-Pex )/abs(Pex) < 1e-3 || abs(Pex_post-Pex) < 1e2
                fprintf('dt=%.2e yr, time=%.3e Myr, Pex_post %.6e Pex %.6e, converged in %d iterations\n',dt/seconds_in_year,(time+dt)/seconds_in_year/1e6,Pex_post,Pex,iter);
                converged = true;
            elseif iter==maxiter
                error('Nonlinear loop failed to converge');
            end

            pex_store(iter) = Pex;
            pexpost_store(iter) = Pex_post;
            if converged
                break;
            end
        end%end nonlinear loop

        % 5. Determine whether tensile failure has occurred
        failure = tensile_failure_criterion(Ro-grid_r',sigma_t,rho,g,tensile_strength);
        if(any(failure)) % failure is occurring
            disp(['Failure criterion has been reached']);
            idx_shallow = find(failure,1,'last');
            idx_deep = find(failure,1,'first');
            fprintf('Shallowest, deepest failure: %f, %f\n\n',Ro-grid_r(idx_shallow),Ro-grid_r(idx_deep));
            fprintf('Failure time: %f Myr\n',time / seconds_in_year / 1e6);
            fprintf('Surface stress at failure: %f MPa\n',sigma_t(end)/1e6);

            % check to see if a crack could propagate to surface
            % 1. Find the midpoint of the crack
            % 2. Look upward - balance stresses on crack in upward
            % direction
            % 3. If crack reached surface, balance stresses on entire
            % crack. Otherwise balance stresses in downward direction.
            sigma_t_tot = sigma_t - rho*g*(Ro-grid_r');
            depth = Ro-grid_r; % depth will be in descending order, i.e. deepest first
            midpoint_depth = mean(depth([idx_shallow idx_deep]));
            [~,midpoint_ind] = max( sigma_t_tot );
            if midpoint_ind == nr
                stress_above = 0;
            else
                stress_above = cumtrapz( grid_r(midpoint_ind:end), sigma_t_tot(midpoint_ind:end) );  % integrate in upward direction
            end
            stress_below = cumtrapz( depth(midpoint_ind:-1:1), sigma_t_tot(midpoint_ind:-1:1) ); % integrate in downward direction
            if stress_above(end) >= 0
                disp('Crack reached surface');
                surface_failure = true;
                above_stress_integral = stress_above(end);
                net_tension = above_stress_integral + stress_below;
                % find depth at which crack stops
                ind = find(net_tension > 0,1,'last'); % net tension is ordered by increasing depth
                depth_tmp = depth(midpoint_ind:-1:1);
                max_depth = depth_tmp(ind); % depth at which crack stops
                min_depth = 0;
                if net_tension > 0
                    disp('Crack reaches ocean!');
                end
            else
                disp('Crack cannot reach surface');
                surface_failure = false;
                % find location where integral of stress is zero in upward
                % direction
                ind = find( stress_above > 0,1,'last');
                depth_tmp = depth(midpoint_ind:end);
                min_depth = depth_tmp(ind);
                % find depth at which crack stops
                ind = find(stress_below > 0,1,'last'); % net tension is ordered by increasing depth
                depth_tmp = depth(midpoint_ind:-1:1);
                max_depth = depth_tmp(ind);
            end
            fprintf('Relieving stresses between %e-%e m\n',min_depth,max_depth);
            results.failure_thickness(ifail) = max_depth-min_depth;
            results.failure_time(ifail) = time/seconds_in_year/1e6;
            results.failure_P(ifail) = Pex;
            results.failure_Pex_crit(ifail) = Pex_crit;

            results.failure_top(ifail) = min_depth;
            results.failure_bottom(ifail) = max_depth;
            results.failure_sigma_t{ifail} = sigma_t;
            results.failure_sigma_r{ifail} = sigma_r;
            results.failure_r{ifail} = grid_r;


            results.failure_z(ifail) = ztmp(end);

            ifail = ifail + 1;
            now_failing = depth >= min_depth & depth <= max_depth;
            failure_mask = failure_mask | now_failing;

            %                 failure_mask = false(size(sigma_r));
            %                 failure_mask(failure) = true;
            failure_time(now_failing) = time+dt;
        else
            no_longer_failing = failure_mask & (time - failure_time) >= 10*dtmin;
            if any(failure_mask(no_longer_failing))
                results.failure_dP(ifail-1) = Pex-results.failure_P(ifail-1);
            end
            if all(failure_mask) && any(failure_mask(no_longer_failing))
                %if erupted_volume > 0
                results.failure_erupted_volume(ifail-1) = erupted_volume;
                results.failure_erupted_volume_volumechange(ifail-1) = erupted_volume_volumechange;
                results.failure_erupted_volume_pressurechange(ifail-1) = erupted_volume_pressurechange;
                %end
                erupted_volume = 0;
                erupted_volume_volumechange = 0;
                erupted_volume_pressurechange = 0;
            end
            failure_mask(no_longer_failing) = false;
        end
        yielding = eiiD > (cohesion - 1/3*ei*friction); % note that compression is negative
        if any(yielding)
            keyboard
        end

        %5.75 consider resetting stresses if ice shell is
        %thinning?
        if time < no_stress_time;
            sigma_r = 0*sigma_r;
            sigma_t = 0*sigma_t;
            siiD = 0*siiD;
            er = 0*er;
            et = 0*et;
            ur = 0*ur;
            Pex = 0.0;
            Tm0 = Tm;
        end


        % 6. advance to next time step and plot (if needed)
        % eccentricity_last = eccentricity;
        sigma_r_last = sigma_r;
        sigma_t_last = sigma_t;
        siiD_last = siiD;
        T_last = T;
        er_last = er;
        et_last = et;
        z_last = z;
        ur_last = ur;
        Pex_last = Pex;
        Tb_last = Tb;
        Tm_last = Tm;

        time = time + dt;

        if (parameters.do_plots && (time >= plot_times(iplot) || time >= t_end ))
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
            h=plot(fig1a.ax(1),(Ro-grid_r)/1e3,sigma_t_last/1e6);
            plot(fig1a.ax(2),(Ro-grid_r)/1e3,T_last,'--','Color',h.Color);


            last_plot_time = time;
            drawnow();
        end
        if (time-last_store >= save_interval || time >= t_end || any(failure_mask))
            sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
            time_store(isave) = time;

            results.time(isave) = time;
            % results.eccentricity(isave) = eccentricity;
            results.thickness(isave) = grid_r(end)-grid_r(1);
            results.z(isave) = z;
            results.Ri(isave) = Ri;
            results.Tm(isave) = Tm;
            % results.qb(isave) = total_heating;
            results.sigma_t(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
            results.sigma_r(:,isave) = interp1(Ro-grid_r,sigma_r_last,save_depths);
            results.e_t(:,isave) = interp1(Ro-grid_r,et_last,save_depths);
            results.e_r(:,isave) = interp1(Ro-grid_r,er_last,save_depths);
            results.ur(:,isave) = interp1(Ro-grid_r,ur_last,save_depths);
            results.ur_base(isave) = ur_last(1);
            results.dTdr(:,isave) = interp1(Ro-grid_r,dTdotdr*dt,save_depths);
            results.T(:,isave) = interp1(Ro-grid_r,T,save_depths);
            results.Tb(isave) = Tb;
            results.Pex(isave) = Pex;

            [strtmp,indtmp] = max(sigma_t-sigma_r);
            results.maximum_differential_stress(isave) = strtmp;
            results.minimum_differential_stress(isave) = min(sigma_t-sigma_r);
            results.maximum_stress_depth(isave) = Ro-grid_r(indtmp);
            [indtmp] = find( sigma_t-sigma_r >= 0,1,'last'); %shallowest value where sigma_t > sigma_r
            
            if indtmp < nr
                ytmp = sigma_t(indtmp:indtmp+1)-sigma_r(indtmp:indtmp+1);
                dydr = diff(ytmp)/(grid_r(indtmp+1)-grid_r(indtmp));
                rtmp = grid_r(indtmp) - ytmp(1)/dydr;
            else
                rtmp = grid_r(indtmp);
            end
            if ~isempty(rtmp)         
                results.stresss_crossover_depth(isave) = Ro-rtmp;
            end
            % results.Pex_crit(isave) = Pex_crit;
            % results.XNH3(isave) = X;
            last_store = time; isave = isave+1;
        end
    end
end

% add quantities needed for plotting.
results.last_save_time=last_store;
results.last_isave = isave;
results.save_depths = save_depths;
results.Ro = Ro;
results.h_crust = h_crust;