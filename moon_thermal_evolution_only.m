clear;
% close all;

% properties of the mantle and crust
% Script to solve coupled ice shell thermal and stress evolution
% Max Rudolph, March 19, 2020
% adapted for the moon, June 2025
%
% Cite:
% Rudolph, M.L., Manga, M., Walker, M., and Rhoden, A. Cooling Crusts
% Create Concommitant Cryovolcanic Cracks. Geophysical Research Letters
% 49(5) e2021GL094421

clear;
% close all;
addpath core; % this is where the helper functions live.
addpath mars;
addpath ~/sw/matlab/crameri

nrs = [256]; % number of points used in the radial direction

for isetup = 6:6
    if isetup == 6 % Moon
        % Settings related to numerics
        label='Moon';
        seconds_in_year = 3.1556952e7;
        max_depth = 8e5; % maximum depth for saving solution values (m)
        relaxation_parameter = 1e-4;%1e-3; % used for fixed point iteration in pressure convergence loop.
        t_end = 4500e6*seconds_in_year;%  3*perturbation_period; 5e8*seconds_in_year;
        dtmax = 5e6*seconds_in_year;
        dtmin = 100*seconds_in_year;%*seconds_in_year;
        no_stress_time = 0.5e9*seconds_in_year; % time before which stresses are not allowed to increase

        % Stuff related to the Mars thermal evolution model
        arh =2.54;   % constant from Michaut equation 12
        C   = 0.5;   % Davaille and Jaupart 1993 constant for heat flux

        % Rheology
        viscosity_model = 2;    % 0 = Nimmo (2004), 1 = Goldsby and Kohlstedt (2001), 2=Arrhenius
        viscosity.d = NaN;      %1e-3; % grain size in m used to calculate the viscosity (for G-K)
        viscosity.P = NaN;      %1e5; % Pressure in MPa used to calculate the viscosity (for G-K)
        Tref = 1600;            % Reference temperature, Kelvin.
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
        thickness = 50e3;       % Initial lid thickness (m)
        Ts = 250;               % Surface temperature, K.
        tstart = 0*seconds_in_year; % starting time of the model, used for radioactive heating.
        % note - initial temperature profile is steady state, calculated later
        % Initial basal temperature
        Tm0 = 1700;             % initial mantle temperature
        DTbl = arh*dTnu(Tm0);
        Tb = Tm0-DTbl;          % Temperature at base of lid

        % Planet properties
        Ro = 1.740e6;      % outer radius of lithosphere
        Ri = Ro-thickness;  % initial inner radius of lithosphere
        Rc = 390*1000;       % core radius, m
        h_crust = 40e3;     % crust thickness (assumed constant)
        crust_heat_fraction = 0.8; % fraction of primitive mantle heat production concentrated within crust
        moon_mass = 0.07346e24;
        silicate_mass = 0.98*moon_mass; % assuming core 2%
        silicate_density = silicate_mass / (4/3*pi*(Ro^3-Rc^3));% density of bulk silicate mars
        % rho = silicate_density; % uniform density approximation
        crust_mass = rhoc*4/3*pi*(Ro^3-(Ro-h_crust)^3);     % mass of the crust
        mantle_mass = silicate_mass - crust_mass; % mass of the mantle
        mantle_density = mantle_mass/( 4/3*pi*((Ro-h_crust)^3-Rc^3));
        rho=mantle_density;
        crust_mass_fraction = crust_mass/silicate_mass;
        % compute the heating per unit mass in the crust
        % mantle heating = [h]*rho*V
        crustal_heating_factor = crust_heat_fraction/crust_mass_fraction; % this is the enrichment in volumetric heating relative to primitive mantle material
        mantle_heating_factor = (1-crust_heat_fraction)/(1-crust_mass_fraction);
        g = 1.625;           % surface gravity (m/s^2)

        kappa = k(Tb)/rho/Cp;           % thermal diffusivity m^2/s        

        % Failure criterion:
        tensile_strength = 1e99; % tensile strength, Pa
        cohesion = 1e99;  % plastic yield strength, MPa
        friction = 0.0; % friction angle for plastic yielding
    end

    if viscosity_model == 0
        label = [label '-nimmovisc'];
    elseif viscosity_model == 1
        label = [label '-goldsbykohlstedt'];
    elseif viscosity_model == 2
        label = [label '-arrhenius']
    elseif viscosity_model == 3
        label = [label '-diffdisl']
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
        results.Tp = zeros(nsave,1);
        results.z_lith = zeros(nsave,1); results.z_lith(1) = Ro-Ri;
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
        results.minimum_differential_stress = NaN*zeros(1,nsave);
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
        H( iscrust ) = moon_heating((time+tstart)/seconds_in_year)*crustal_heating_factor*rhoc;
        H(~iscrust ) = moon_heating((time+tstart)/seconds_in_year)*mantle_heating_factor*rho;
        kvec = k(Tm0)*ones(nr,1);
        kvec(iscrust) = kcrust;
        Cpvec = Cp*ones(nr,1);
        rhovec = rho*ones(nr,1); rhovec(iscrust) = rhoc;

        [T_last,dTdotdr] = solve_temperature_shell_mars(grid_r,T_last,Tb,Ts,kvec,rhovec,Cpvec,H,Inf,0.0); % call solver with infinite timestep
        %T_last = Tb+(Ts-Tb)/(grid_r(end)-grid_r(1))*(grid_r-grid_r(1))';
        % dTdotdr = zeros(size(T_last));
        Tm = Tm0; % mantle temperature

        er_last = zeros(nr,1); % strains
        et_last = zeros(nr,1);
        ur_last = zeros(nr,1); % displacement
        z_last = 0;    % total amount of thickening
        dzdt_last = 0; % thickening rate
        Pex_last = 0; %initial overpressure

        % Set up plot
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
            D = Ro-(Ri-z_last); % z is the amount by which the lid has thickened
            mantle_volume = 4/3*pi*((Ro-D)^3-Rc^3);
            Cm = rho*Cp*mantle_volume; % mantle heat capacity
            Slid = 4*pi*(Ro-D)^2;
            h_conv = mantle_heating_factor*rho*moon_heating(time/seconds_in_year); % mantle volumetric heat production
            % boundary layer heat transport into the lid:
            qbl = C*k(Tm)*(alpha_v_bl*rho*g/kappa/mu(Tm,1e9,0))^(1/3)*dTnu(Tm)^(4/3);% note assumes 1 GPa-pressure creep viscosity
            % temperature difference across the boundary layer:
            DTbl = arh*dTnu(Tm);
            delta_bl = k(Tm)*DTbl/qbl; % boundary layer thickness
            Tl = Tm-DTbl;% temp at base of conductive layer
            dH = rho*Cp*DTbl; % enthalpy change across the lid
            % rate of change of mantle temperature
            dTmdt = 1/Cm * (-Slid*qbl + h_conv*mantle_volume); % Michaut et al. Equation 18
            % rate of change of lid thickness.
            dDdt = 1/dH * (qlid-qbl);   % Michaut et al. Equation 19

            % determine the timestep - apply a courant type condition
            % to lid thickness change
            if abs(dDdt*dt) > (grid_r(2)-grid_r(1))/2
                dt = abs( (grid_r(2)-grid_r(1))/2/(dDdt) );
            end
            % apply a limiter based on mantle temperature change
            max_dTm = 0.1;
            if abs(dTmdt*dt) > max_dTm
                dt = abs(max_dTm/dTmdt);
            end
            if dt < dtmin
                dt = dtmin;
                warning('Setting dt = dtmin');
            end
            if any(failure_mask)
                dt = dtmin;
            end
            % dTldt = (Tm-DTbl - T_last(1))/dt; %rate of change of temperature at base of lid.


            % update the mantle temperature
            delta_rb = dDdt*dt;
            z = z_last + delta_rb;
            % update the basal temperature
            Tm = Tm + dTmdt*dt;
            Tb = Tm - arh*dTnu(Tm);

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
            H( iscrust ) = moon_heating((time+tstart)/seconds_in_year)*crustal_heating_factor*rhoc;
            H(~iscrust ) = moon_heating((time+tstart)/seconds_in_year)*mantle_heating_factor*rho;
            kvec = k(Tm0)*ones(nr,1);
            kvec(iscrust) = kcrust;
            Cpvec = Cp*ones(nr,1);
            rhovec = rho*ones(nr,1); rhovec(iscrust) = rhoc;

            [T,dTdotdr] = solve_temperature_shell_mars(grid_r,T_last,Tb,Ts,kvec,rhovec,Cpvec,H,dt,delta_rb);

            %5.75 consider resetting stresses if ice shell is
            %thinning?

            % 6. advance to next time step and plot (if needed)

            T_last = T;
            z_last = z;
            Tb_last = Tb;
            Tm_last = Tm;

            % compute the mantle potential temperature
            z_lith = (Ro-Ri)+z + delta_bl;% lithosphere thickness = lid thickness + boundary layer thickness
            Tp = Tm/exp(alpha_v*g*z_lith/Cp);
            if time == 0
                results.Tp(1) = Tp;
            end

            time = time + dt;

            if (time >= plot_times(iplot) || time >= t_end )
                iplot = iplot+1;

                figure(hf);
                subplot(1,4,1);
                subplot(1,4,3);
                plot(T,Ro-grid_r);


                figure(hf2);

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
                results.z_lith(isave)=z_lith;
                results.Tp(isave) = Tp;
                % results.qb(isave) = total_heating;

                results.dTdr(:,isave) = interp1(Ro-grid_r,dTdotdr*dt,save_depths);
                results.T(:,isave) = interp1(Ro-grid_r,T,save_depths);
                results.Tb(isave) = Tb;
                last_store = time; isave = isave+1;
            end
        end
        %%

        % compute (sigma_t-sigma_r)


        %% Pseudocolor stress plot

        mask = ~isnan(results.time);

       t=tiledlayout(4,1);

        % for i=1:ifail-1
        %     plot(results.failure_time(i)*1e6*[1 1],[results.failure_top(i) results.failure_bottom(i)]/1e3,'r');
        % end
        % TEMPERATURE
        nexttile
        contourf(results.time(mask)/seconds_in_year/1e6,save_depths/1000,results.T(:,mask),64,'Color','none'); %shading flat;
        hold on
        plot(results.time(mask)/seconds_in_year/1e6,((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        %         set(gca,'YLim',[0 ceil(1+max(((Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        % ax1 = gca();
        % ax1.FontSize=8;
        set(gca,'Colormap',crameri('-lajolla'))
        hcb = colorbar();
        hcb.Label.String = 'Temperature (K)';
        text(0.025,0.85,char('C'),'FontSize',12,'Units','normalized');
        % xlabel('Time (years)');
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        hold on;

        nexttile
        plot(results.time(mask)/seconds_in_year/1e6,results.ur(1,mask),'k')
        ylabel('u_r (m)')
        text(0.025,0.85,char('D'),'FontSize',12,'Units','normalized');

        nexttile
        plot(results.time(mask)/seconds_in_year/1e6,results.Pex(mask)/1e6,'k');
        ylabel('P_{ex} (MPa)');
        set(gca,'XScale',xscale);
        % ax2 = gca();
        % ax2.Position(3) = ax1.Position(3);
        % ax2.XLim = ax1.XLim;
        % ax2.FontSize=8;
        hold on
        % plot(results.failure_time(1:ifail-1)*1e6,results.failure_P(1:ifail-1)/1e6,'r.');
        % end_color = [0 0.9 0];
        % plot(results.failure_time(1:ifail-1)*1e6,(results.failure_P(1:ifail-1)+results.failure_dP(1:ifail-1))/1e6,'LineStyle','none','Color',end_color,'Marker','o','MarkerFaceColor',end_color,'MarkerSize',2);
        text(0.025,0.85,char('E'),'FontSize',12,'Units','normalized');
        % plot(results.time(mask)/seconds_in_year,results.Pex_crit(mask)/1e6,'k-');

        % xlabel('Time (years)');
        nexttile
        % for i=1:ifail-1
        %     if isnan(results.failure_erupted_volume(i))
        %         % plot nothing
        %     else
        %         if results.failure_P(i) - results.failure_Pex_crit(i) > 0
        %             plot(results.failure_time(i)*1e6*[1 1],[0 1],'b');
        %         else
        %             plot(results.failure_time(i)*1e6*[1 1],[0 1],'b--');
        %         end
        %     end
        % end
        plot(results.time(mask)/seconds_in_year/1e6,results.Tm(mask),'k-');

        ylabel('T_m (K)');
        xlabel('Time (years)');
        set(gca,'XScale',xscale);
        % ax3=gca();
        % ax3.XLim = ax1.XLim;
        % ax3.Position(3) = ax1.Position(3);
        % ax3.Box = 'on';
        % ax3.FontSize=8;
        text(0.025,0.85,char('F'),'FontSize',12,'Units','normalized');
        linkaxes(t.Children,'x');
        set(gca,'XLim',[0 4500]);


        fig = gcf();
        fig.Position(3:4) = [385   650];
        axmask = arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'),t.Children);

        set(t.Children(axmask),'XTickLabel',[])
        set(gca,'XTickLabel',get(gca,'XTick'))

        fig.Color = 'w';
        filename = sprintf('moon-thermal-evolution-zerotime-%f.pdf',no_stress_time/seconds_in_year/1e9);
        % exportgraphics(gcf,filename,'ContentType','vector');

        %% new multi-panel plot
        tga = 4.5-results.time/seconds_in_year/1e9;
        mask1 = mask & results.time>(no_stress_time+save_interval*10);

        xscale = 'linear';
        ax=[];
        figure();
        t=tiledlayout(5,1,'TileSpacing','compact','Padding','none');

        nexttile([2,1])
        contourf(tga(mask),save_depths/1000,results.differential_stress(:,mask)/1e6,64,'Color','none'); %shading flat;
        contourf(tga(mask),save_depths/1000,results.sigma_t(:,mask)/1e6,64,'Color','none'); %shading flat;
        hold on
        plot(tga(mask),((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        plot(tga(mask),results.z_lith(mask)/1000,'--','Color','k','LineWidth',1);
        hold on
        contour(tga(mask),save_depths/1000,results.T(:,mask),[1000 1000],'Color','k','LineStyle','-'); %


        set(gca,'YDir','reverse');
        hcb = colorbar();
        set(gca,'Colormap',crameri('-roma'))
        stmax = max(max(abs(results.sigma_t(:,mask)/1e6)));
        caxis([-1 1]*stmax)
        hcb.Label.String = '\sigma_t-\sigma_r (MPa)';
        text(-0.16,0.95,char('A'+0),'FontSize',12,'Units','normalized');

        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        hold on;

        nexttile
        contourf(tga(mask),save_depths/1000,results.T(:,mask),64,'Color','none'); %shading flat;
        hold on
        plot(tga(mask),((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        set(gca,'YDir','reverse');
        % ax1 = gca();
        % ax1.FontSize=8;
        set(gca,'Colormap',crameri('-lajolla'))
        hcb = colorbar();
        hcb.Label.String = 'Temperature (K)';
        text(-0.16,0.95,char('A'+1),'FontSize',12,'Units','normalized');
        % xlabel('Time (years)');
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        hold on;
        %
        % e_t
        %
        nexttile
        plot(tga(mask1),results.e_t(1,mask1)*1e3,'k');
        ylabel('\epsilon_t (10^{-3})')
        text(-0.16,0.95,char('A'+2),'FontSize',12,'Units','normalized');

        %
        % Tm
        %
        nexttile


        plot(tga(mask),results.Tm(mask),'r-','LineWidth',1);
        hold on
        plot(tga(mask),results.Tp(mask),'k','LineWidth',1);
        text(-0.16,0.95,char('A'+3),'FontSize',12,'Units','normalized');


        set(gca,'Box','on')

        ylabel('T_m (K)');
        xlabel('Time (Ga)');
        set(gca,'XScale',xscale);

        fig = gcf();
        fig.Position(3:4) = [385   650];
        axmask = arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'),t.Children);
        linkaxes(t.Children(axmask),'x');
        set(gca,'XLim',[0 4.5]);

        set(t.Children(axmask),'XTickLabel',[])
        set(gca,'XTickLabel',get(gca,'XTick'))
        set(t.Children(axmask),'XDir','reverse')

        fig.Color = 'w';
        filename = sprintf('mars-thermal-evolution-zerotime-%f.pdf',no_stress_time/seconds_in_year/1e9);
        % exportgraphics(gcf,filename,'ContentType','vector');



    end
end

