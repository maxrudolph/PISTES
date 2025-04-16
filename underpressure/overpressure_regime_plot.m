clear;
close all;

% moon = 'enceladus';

rhoi = 910;
rhow = 1000.0;
rhoc = 3000;% density of rock (core)
G = 6.67430e-11; % gravitational constant, SI

beta = 4e-10; % ocean compressibility Pa^-1
E = 5e9; % Young's modulus, Pa
nu = 0.33;
% 
% % Define data structures for each satellite
% europa.name = 'Europa';
% europa.R = 1560.8e3;
% europa.ri = europa.R - 2.4e3;  % Manga and Wang 2007 value
% europa.rc = europa.R - 1.20e5; % Manga and Wang, cites Cammarano et al. 2006
% europa.g = 1.3;
% europa.style = '--';
% 
% enceladus.name = 'Enceladus';
% enceladus.R = 2.52e5; % outer radius
% enceladus.ri = enceladus.R - 5e4; % Manga and Wang 2007 value.
% enceladus.rc = 1.61e5;
% enceladus.g = 0.113;
% enceladus.style = '-';
% 
% mimas.name = 'Mimas';
% mimas.R = 198.2e3;
% mimas.rc = mimas.R-1.266e5;
% mimas.ri = mimas.rc;
% mimas.g = 0.064;
% mimas.style = ':';
% 
% charon.name = 'Charon';
% charon.R = 606e5;
% charon.rc = 3.76e5;
% charon.ri = charon.rc;
% charon.g = 0.288;
% charon.style = '-.';

% get_water_fraction = @(planet) 

initial_frozen_fraction = NaN;%0.05; % the initial fraction of the H2O that is frozen
initial_frozen_thickness= 1e3;
% cohesion = 4e7;
tensile_strength = 3e6;

elastic_fraction = 1/2;


label = sprintf('_initial-%f_tensilestrength-%e',initial_frozen_fraction,tensile_strength);

% functions for excess pressure & tensile stress
Pex = @(z_values,ri_values,xi_values,planet) (z_values .* (1-rhoi/rhow)) ./ (beta*(ri_values.^3-planet.rc^3)./(3*ri_values.^2) + ...
    xi_values/E.*(1 + 2*nu*(1+0.5*(planet.R./xi_values).^3)./( (planet.R./xi_values).^3-1) ) );
sigma_t = @(z_values,ri_values,xi_values,planet) 3/2*Pex(z_values,ri_values,xi_values,planet)./((planet.R./xi_values).^3-1);
sigma_rr = @(Pex1,r,ri,xi,planet) Pex1/((planet.R/xi)^3-1)*(1-(planet.R./r).^3);
sigma_tt = @(Pex1,r,ri,xi,planet) Pex1/((planet.R/xi)^3-1)*(1+0.5*(planet.R./r).^3);

nf = 21;
nr = 20;

ff = linspace(0.001,0.999,nf);     % mass fraction water (+ice)
RR = linspace(1.0e5,1.8e6,nr);     % moon radius
% fractional_thickening_boil = NaN*ones(nf,nr);
% fractional_thickening_yield_constant = NaN*ones(nf,nr);
% fractional_thickening_yield_coulomb = NaN*ones(nf,nr);
fractional_thickening_fail_tension = NaN*ones(nf,nr);
% yield_constant = NaN*ones(nf,nr);
% yield_coulomb = NaN*ones(nf,nr);
fail_tension = NaN*ones(nf,nr);
% sdmax_boil = NaN*ones(nf,nr);
% sdmax_yield_constant = NaN*ones(nf,nr);
% sdmax_yield_coulomb = NaN*ones(nf,nr);
% sdmax_fail_tension = NaN*ones(nf,nr);
net_loading_at_failure = NaN*ones(nf,nr);
Pex_at_failure = NaN*ones(nf,nr);
Pexcrit_at_failure = NaN*ones(nf,nr);

h2o_thickness = zeros(nf,nr);
sdmax = NaN*zeros(nf,nr);
boil = NaN*zeros(nf,nr);

for i = 1:nf
    f = ff(i);
    for j=1:nr
        clear p
        p.R = RR(j);
        % thickness of water layer
        M = 4/3*pi*p.R^3 / ( (1-f)/rhoc + f/rhow );
        p.rc = (3*(1-f)*M/(4*pi*rhoc) )^(1/3);
        p.g = G*M/p.R^2;
        % sanity check:
        Mw = 4/3*pi*(p.R^3 - p.rc^3)*rhow;
        fw = Mw/M;
        gg(i,j) = p.g;
        h2o_thickness(i,j) = p.R-p.rc;

        nthick = 50; % number of thickness values to test
        % assume a constant initial thickness
        z_values = linspace(0,(p.R-p.rc-initial_frozen_thickness),nthick);
        ri_values = p.R-initial_frozen_thickness - z_values;

        % assume that the ice shell is X%  frozen:
        % z_values = linspace(0,(p.R-p.rc)*(1-initial_frozen_fraction)-1,nthick);  % amount of thickening (thinning negative)
        % ri_values = p.R - (p.R-p.rc)*(initial_frozen_fraction) - z_values;       % coordinate of base of ice shell

        elastic_thickness = (elastic_fraction)*(p.R-ri_values);  % elastic thickness
        xi_values = p.R - elastic_thickness;        % coordinate of base of elastic layer
        % compute the pressure
        Pex_values = Pex(z_values,ri_values,xi_values,p);
        P_values = rhoi*p.g*(p.R-ri_values) + Pex_values;
        sigma_t_values = sigma_t(z_values,ri_values,xi_values,p);
        Pex_crit_values = (p.R-ri_values)*p.g*(rhow-rhoi);

        % more correct but costly check for yielding
        % ind1 = []; % index for yielding via mohr-columb
        % ind2 = []; % index for yielding at constant strength
        ind3 = []; % index for tensile failure
        % fail1_values = zeros(nthick,1);
        % fail2_values = zeros(nthick,1);
        fail3_values = zeros(nthick,1);
        % fail1_depth = zeros(nthick,1);
        % fail2_depth = zeros(nthick,1);
        fail3_depth = zeros(nthick,1);
        % sdmax_values = zeros(nthick,1);
        % sdmax_depth = zeros(nthick,1);
        stensile_max = zeros(nthick,1);
        net_loading = zeros(nthick,1);

        for k=1:nthick % loop over thickening amounts (z_values)
            Pex1 = Pex_values(k);
            rr = linspace(xi_values(k),p.R,100);% start at base of elastic layer, go to surface.
            srr = sigma_rr(Pex1,rr,ri_values(k),xi_values(k),p);
            stt = sigma_tt(Pex1,rr,ri_values(k),xi_values(k),p);
            % Pex,r,z,ri,xi,planet
            pp = 1/3*(srr+2*stt); % elastic pressure
            plith = -rhoi*p.g*(p.R-rr); % lithostatic stress (compression negative)
            ptot = pp+plith;
            % sd = srr-stt; % differential stress
            stensile = stt+plith;

            rrp = linspace(ri_values(k),xi_values(k));
            plithp = -rhoi*p.g*(p.R-rrp);
            tmp1 = cumtrapz(rr,plith);
            tmp2 = cumtrapz(rrp,plithp);
            pressure_integral = tmp1(end)+tmp2(end);
            loading_integral = cumtrapz(rr,stt);
            loading_integral = loading_integral(end);
            net_loading(k) = loading_integral + pressure_integral;

            [stmp,itmp] = max(stensile); % this is (tensile stress)-(lithostatic pressure)
            fail3_values(k) = stmp;
            fail3_depth(k) = itmp;           

            % if isempty(ind1) && any( sd > strength )
                % ind1 = k;
            % end
            % if isempty(ind2) && any( sd > cohesion )
                % ind2 = k;
            % end
            if any( stensile > tensile_strength )
                ind3 = k;              
            end
            % if ~isempty(ind1) && ~isempty(ind2)
                % break;
            % end
           
            if ~isempty(ind3) % tensile failure has occurred
                break;
            end
        end
        
        fail3_values = fail3_values(1:k);
        fail3_depth = fail3_depth(1:k);
        net_loading = net_loading(1:k);
        
        if ~isempty(ind3)
            
            fail_tension(i,j) = interp1(fail3_values,1:k,tensile_strength); % this is the interpolated index of pressure where failure occurred
            Pex_at_failure(i,j) = interp1(1:k,Pex_values(1:k),fail_tension(i,j));
            Pexcrit_at_failure(i,j) = interp1(1:k,Pex_crit_values(1:k),fail_tension(i,j));
            net_loading_at_failure(i,j) = interp1(1:k,net_loading,fail_tension(i,j));
            fractional_thickening_fail_tension(i,j) = interp1(1:nthick,z_values,fail_tension(i,j))/(p.R-p.rc);
        end     
    end
end

%% load table of satellites
TableofSatellites = import_satellite_table('Table_of_Satellites.xlsx', "Sheet1", [1, Inf]);
% This table assumes a core density of 3500 and ice+water density of 1000
satellite_f = 1-TableofSatellites.McMp;
satellite_R = TableofSatellites.RcKm./TableofSatellites.RcRp;

%% Plotting
figure
contourf(RR/1e3,ff,gg,128);
hold on
scatter(satellite_R,satellite_f,'ko','MarkerFaceColor','k');
text(satellite_R+25,satellite_f,TableofSatellites.Satellite);
shading flat
colorbar()

figure()
% pcolor(RR,ff,h2o_thickness);
hold on
contour(RR,ff,h2o_thickness,25000:25000:1000000);
shading flat
colorbar()
title('h2o thickness')

extrusion_possible = Pex_at_failure > Pexcrit_at_failure;

f=figure(101);
% pos = get(gcf,'Position');
% f.Position(3) = pos(3)*2;
% subplot(1,2,1);
% pcolor(RR/1e3,ff,fractional_thickening_fail_tension);
contourf(RR/1e3,ff,net_loading_at_failure,128,'Color','none');
shading flat;
hcb=colorbar();
hcb.Label.String = 'Fractional Thinning to Fail';
xlabel('Planetary Body Radius (km)');
ylabel('H_20 Mass Fraction (-)')
hold on
scatter(satellite_R,satellite_f,'ko','MarkerFaceColor','k');
text(satellite_R+25,satellite_f,TableofSatellites.Satellite);
contour(RR/1e3,ff,net_loading_at_failure,[0.0 0.0],'k','LineWidth',1);
contour(RR/1e3,ff,extrusion_possible,[0.5 0.5],'r','LineWidth',1);
% contour(RR/1e3,ff,boil_region_constant,[0.5 0.5],'k--','LineWidth',1);
set(gca,'FontSize',14);
set(gca,'Layer','top');
% caxis([0 1]);
set(gca,'YLim',[0 1]);
text(0.05,0.95,'A','Units','normalized','FontSize',18)
% set(gcf,'Color','none');
% exportgraphics(gcf,['fractional_thickening' label '.pdf'],'Resolution',600);

% figure
% pcolor(RR/1e3,ff,yield_coulomb);
% shading flat;
% hold on
% % contour(RR/1e3,ff,yield,[8 40]*10^6,'k');
% hcb=colorbar();
% % set(gca,'ColorScale','log');
% hcb.Label.String = 'yield criterion?';
% set(gca,'FontSize',14)
% set(gca,'Layer','top')
%% 
figure();
% subplot(1,2,2);
max_stress = sdmax_yield_coulomb;
max_stress(boil_region_coulomb) = sdmax_boil(boil_region_coulomb);

hh=pcolor(RR/1e3,ff,max_stress);
shading flat;
hold on
hcb=colorbar();
set(gca,'ColorScale','log');
hcb.Label.String = 'Maximum differential stress (Pa)';
% contour(RR/1e3,ff,sdmax,[8 40]*10^6,'k');
contour(RR/1e3,ff,boil_region_coulomb,[0.5 0.5],'k','LineWidth',1);
% contour(RR/1e3,ff,yield_coulomb,[0.5 0.5],'k--','LineWidth',1);
set(gca,'FontSize',14)

% hcb.Label.String = 'Fractional Thinning to Boil';
xlabel('Planetary Body Radius (km)');
ylabel('H_2O Mass Fraction (-)')
hold on
scatter(satellite_R,satellite_f,'ko','MarkerFaceColor','k');
text(satellite_R+25,satellite_f,TableofSatellites.Satellite);
set(gca,'FontSize',14);
set(gca,'Layer','top');
set(gca,'YLim',[0 1]);
text(0.05,0.95,'B','Units','normalized','FontSize',18);
set(gcf,'Color','none');
exportgraphics(gcf,['differential_stress_coulomb' label '.pdf'],'Resolution',600);


% CONSTANT YIELD STRESS
figure()
max_stress = sdmax_yield_constant;
max_stress(boil_region_constant) = sdmax_boil(boil_region_constant);

hh=pcolor(RR/1e3,ff,max_stress);
shading flat;
hold on
hcb=colorbar();
set(gca,'ColorScale','log');
hcb.Label.String = 'Maximum differential stress (Pa)';
% contour(RR/1e3,ff,sdmax,[8 40]*10^6,'k');
contour(RR/1e3,ff,boil_region_constant,[0.5 0.5],'k--','LineWidth',1);
set(gca,'FontSize',14)
xlabel('Planetary Body Radius (km)');
ylabel('H_2O Mass Fraction (-)')
hold on
scatter(satellite_R,satellite_f,'ko','MarkerFaceColor','k');
text(satellite_R+25,satellite_f,TableofSatellites.Satellite);
set(gca,'FontSize',14);
set(gca,'Layer','top');
set(gca,'YLim',[0 1]);
text(0.05,0.95,'B','Units','normalized','FontSize',18)
exportgraphics(gcf,'differential_stress_constant.pdf','Resolution',600);
%% 

figure, 

% make the plot
nthick = 1000; % number of thicknesses

Pex = @(z_values,ri_values,xi_values,planet) (z_values .* (1-rhoi/rhow)) ./ (beta*(ri_values.^3-planet.rc^3)./(3*ri_values.^2) + ...
    xi_values/E.*(1 + 2*nu*(1+0.5*(planet.R./xi_values).^3)./( (planet.R./xi_values).^3-1) ) );
sigma_t = @(z_values,ri_values,xi_values,planet) 3/2*Pex(z_values,ri_values,xi_values,planet)./((planet.R./xi_values).^3-1);

% reproduce Manga and Wang 2007 results - 1/3 elastic fraction
% p = enceladus;
% z_values = logspace(1,4,nthick); % amount of thickening
% ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
% elastic_thickness = 1/3*(p.R-(ri_values));
% xi_values = p.R-elastic_thickness;
% Pex_13 = Pex(z_values,ri_values,xi_values,p);
% sigmat_13 = sigma_t(z_values,ri_values,xi_values,p);
% Pb_13 = rhoi*p.g*(p.R-ri_values) + Pex_13;
% 
% ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
% elastic_thickness = 1e3; % fixed at 1 km
% xi_values = p.R-elastic_thickness;
% Pex_1km = Pex(z_values,ri_values,xi_values,p);
% sigmat_1km = sigma_t(z_values,ri_values,xi_values,p);
% Pb_1km = rhoi*p.g*(p.R-ri_values) + Pex_1km;
% 
% figure;
% plot(z_values,Pex_13,'k--','DisplayName','Pex 1/3');
% hold on
% plot(z_values,sigmat_13,'k','DisplayName','\sigma_t 1/3');
% plot(z_values,Pex_1km,'r--','DisplayName','Pex 1km');
% plot(z_values,sigmat_1km,'r','DisplayName','\sigma_t 1km');
% 
% legend('Location','northwest');
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% title(p.name);
% 
% figure;
% plot(z_values,Pb_13,'DisplayName','Pb 1/3')
% hold on
% plot(z_values,Pb_1km,'DisplayName','Pb 1 km')
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% legend()
% 
% %% Manga and Wang - Europa
% p = europa;
% z_values = logspace(1,4,nthick); % amount of thickening
% ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
% elastic_thickness = 1/3*(p.R-(ri_values));
% xi_values = p.R-elastic_thickness;
% Pex_13 = Pex(z_values,ri_values,xi_values,p);
% sigmat_13 = sigma_t(z_values,ri_values,xi_values,p);
% Pb_13 = rhoi*p.g*(p.R-ri_values) + Pex_13;
% 
% ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
% elastic_thickness = 1e3; % fixed at 1 km
% xi_values = p.R-elastic_thickness;
% Pex_1km = Pex(z_values,ri_values,xi_values,p);
% sigmat_1km = sigma_t(z_values,ri_values,xi_values,p);
% Pb_1km = rhoi*p.g*(p.R-ri_values) + Pex_1km;
% 
% figure;
% plot(z_values,Pex_13,'k--','DisplayName','Pex 1/3');
% hold on
% plot(z_values,sigmat_13,'k','DisplayName','\sigma_t 1/3');
% plot(z_values,Pex_1km,'r--','DisplayName','Pex 1km');
% plot(z_values,sigmat_1km,'r','DisplayName','\sigma_t 1km');
% 
% legend('Location','northwest');
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% title(p.name);
% 
% figure;
% plot(z_values,Pb_13,'DisplayName','Pb 1/3')
% hold on
% plot(z_values,Pb_1km,'DisplayName','Pb 1 km')
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% legend(['total pressure ' p.name])
% 
% 
% %% thinning ice shell...
% z_values = -linspace(1,5e4,nthick); % amount of thickening
% 
% figure(1001); clf;
% figure(1002); clf;
% for p = [enceladus europa mimas charon]
%     p.ri = p.rc;
%     z_values = linspace(0,-((p.R-1e3-p.ri)),nthick); % amount of thickening - limited by ice shell thickness
% 
%     ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
%     elastic_thickness = 1/3*(p.R-(ri_values));
%     xi_values = p.R-elastic_thickness;
%     Pex_13 = Pex(z_values,ri_values,xi_values,p);
%     sigmat_13 = sigma_t(z_values,ri_values,xi_values,p);
%     Pb_13 = rhoi*p.g*(p.R-ri_values) + Pex_13;
% 
%     % ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
%     % elastic_thickness = 1e3; % fixed at 1 km
%     % xi_values = p.R-elastic_thickness;
%     % Pex_1km = Pex(z_values,ri_values,xi_values,p);
%     % sigmat_1km = sigma_t(z_values,ri_values,xi_values,p);
%     % Pb_1km = rhoi*p.g*(p.R-ri_values) + Pex_1km;
% 
%     figure(1001);
%     plot(z_values,Pex_13,'k','LineStyle',p.style,'DisplayName',[p.name ' 1/3']);
%     hold on
%     plot(z_values,sigmat_13,'k','DisplayName','\sigma_t 1/3');
%     % plot(z_values,Pex_1km,'r--','DisplayName','Pex 1km');
%     % plot(z_values,sigmat_1km,'r','DisplayName','\sigma_t 1km');
% 
%     figure(1002);
%     Pb0 = p.g*rhoi*(p.R-p.rc);
%     R0 = (p.R-p.rc);
%     mask = Pb_13 > 600;
%     plot(-z_values(mask)/R0,(Pb_13(mask)-(600))/(Pb0),'DisplayName',[p.name ' 1/3']);
%     hold on
%     % plot(z_values,Pb_1km,'DisplayName','Pb 1 km');
% end
% figure(1001);
% legend('Location','northwest');
% % set(gca,'XScale','log');
% % set(gca,'YScale','log');
% title(p.name);
% 
% 
% 
% % set(gca,'XScale','log');
% % set(gca,'YScale','log');/
% 
% figure(1002);
%     plot(get(gca,'XLim'),[0 0],'k')
%     ylabel('(P-P_{tp})/P_c (-)')
% xlabel('Fracional thinning (-)')
% legend()
% title('Basal Total Pressure')
% 
% % plot the excess pressure
% % plot elastic stresses in the ice shell
% 
% % plot the pressure profile
% 
