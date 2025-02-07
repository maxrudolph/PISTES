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

initial_frozen_fraction = 1.0; % the initial fraction of the H2O that is frozen
cohesion = 2e7;

label = sprintf('_initial-%f_strength-%e',initial_frozen_fraction,cohesion);

% functions for excess pressure & tensile stress
Pex = @(z_values,ri_values,xi_values,planet) (z_values .* (1-rhoi/rhow)) ./ (beta*(ri_values.^3-planet.rc^3)./(3*ri_values.^2) + ...
    xi_values/E.*(1 + 2*nu*(1+0.5*(planet.R./xi_values).^3)./( (planet.R./xi_values).^3-1) ) );
sigma_t = @(z_values,ri_values,xi_values,planet) 3/2*Pex(z_values,ri_values,xi_values,planet)./((planet.R./xi_values).^3-1);
sigma_rr = @(Pex1,r,ri,xi,planet) Pex1/((planet.R/xi)^3-1)*(1-(planet.R./r).^3);
sigma_tt = @(Pex1,r,ri,xi,planet) Pex1/((planet.R/xi)^3-1)*(1+0.5*(planet.R./r).^3);

nf = 501;
nr = 500;

ff = linspace(0.001,0.999,nf);     % mass fraction water (+ice)
RR = linspace(1.0e5,1.8e6,nr);     % moon radius
fractional_thickening_boil = NaN*ones(nf,nr);
fractional_thickening_yield_constant = NaN*ones(nf,nr);
fractional_thickening_yield_coulomb = NaN*ones(nf,nr);
yield_constant = NaN*ones(nf,nr);
yield_coulomb = NaN*ones(nf,nr);
sdmax_boil = NaN*ones(nf,nr);
sdmax_yield_constant = NaN*ones(nf,nr);
sdmax_yield_coulomb = NaN*ones(nf,nr);


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

        nthick = 100; % number of thickness values to test
        % assume that the ice shell is entirely frozen:
        z_values = linspace(0,-(p.R-p.rc-1),nthick);  % amount of thickening (thinning negative)
        ri_values = p.rc + z_values;                  % coordinate of base of ice shell

        % assume that the ice shell is X%  frozen:
        z_values = linspace(0,-(p.R-p.rc)*(initial_frozen_fraction)+1,nthick);  % amount of thickening (thinning negative)
        ri_values = p.rc + (p.R-p.rc)*(1-initial_frozen_fraction) + z_values;     % coordinate of base of ice shell

        elastic_thickness = (1/2)*(p.R-ri_values);  % elastic thickness
        xi_values = p.R - elastic_thickness;        % coordinate of base of elastic layer
        % compute the pressure
        Pex_values = Pex(z_values,ri_values,xi_values,p);
        P_values = rhoi*p.g*(p.R-ri_values) + Pex_values;
        sigma_t_values = sigma_t(z_values,ri_values,xi_values,p);
        
        % more correct but costly check for yielding
        ind1 = []; % index for yielding via mohr-columb
        ind2 = []; % index for yielding at constant strength
        fail1_values = zeros(nthick,1);
        fail2_values = zeros(nthick,1);
        fail1_depth = zeros(nthick,1);
        fail2_depth = zeros(nthick,1);
        sdmax_values = zeros(nthick,1);
        sdmax_depth = zeros(nthick,1);
        for k=1:nthick
            Pex1 = Pex_values(k);
            rr = linspace(xi_values(k),p.R,100);% start at base of elastic layer, go to surface.
            srr = sigma_rr(Pex1,rr,ri_values(k),xi_values(k),p);
            stt = sigma_tt(Pex1,rr,ri_values(k),xi_values(k),p);
            % Pex,r,z,ri,xi,planet
            pp = 1/3*(srr+2*stt); % elastic pressure
            plith = -rhoi*p.g*(p.R-rr); % lithostatic stress (compression negative)
            ptot = pp+plith;
            sd = srr-stt; % differential stress
            strength = -ptot*0.6 + cohesion;% 40 MPa + 0.6*P

            [stmp,itmp] = max(sd-strength);
            fail1_values(k) = stmp;
            fail1_depth(k) = itmp;

            [stmp,itmp] = max(sd-cohesion);
            fail2_values(k) = stmp;
            fail2_depth(k) = itmp;

            [stmp,itmp] = max(sd);
            sdmax_depth(k) = itmp; 
            sdmax_values(k) = stmp;

            if isempty(ind1) && any( sd > strength )
                ind1 = k;
            end
            if isempty(ind2) && any( sd > cohesion )
                ind2 = k;
            end
            if ~isempty(ind1) && ~isempty(ind2)
                break;
            end
        end
        fail1_values = fail1_values(1:k);
        fail1_depth = fail1_depth(1:k);
        fail2_values = fail2_values(1:k);
        fail2_depth = fail2_depth(1:k);
        sdmax_values = sdmax_values(1:k);
        sdmax_depth = sdmax_depth(1:k);

        [ind] = find(P_values <= 600,1,'first'); % total pressure less than boiling pressure
        if ~isempty(ind)
            boil(i,j) = interp1(P_values,1:nthick,600);% interpolated index where pressure reaches 600 Pa            
            fractional_thickening_boil(i,j) = -interp1(1:nthick,z_values,boil(i,j))/(p.R-p.rc);
            sdmax_boil(i,j) = interp1(1:k,sdmax_values,boil(i,j));
        end      
        if ~isempty(ind1)
            % yield_coulomb(i,j) = ind1;
            yield_coulomb(i,j) = interp1(fail1_values,1:k,0);
            fractional_thickening_yield_coulomb(i,j) = -interp1(1:nthick,z_values,yield_coulomb(i,j))/(p.R-p.rc);
            sdmax_yield_coulomb(i,j) = interp1(1:k,sdmax_values,yield_coulomb(i,j));
            yield_depth_coulomb(i,j) = fail1_depth(ind1);
        end
        if ~isempty(ind2)
            yield_constant(i,j) = interp1(fail2_values,1:k,0);
            fractional_thickening_yield_constant(i,j) = -interp1(1:nthick,z_values,yield_constant(i,j))/(p.R-p.rc);
            sdmax_yield_constant(i,j) = interp1(1:k,sdmax_values,yield_constant(i,j));
            yield_depth_constant(i,j) = fail2_depth(ind2);
        end
        % [ind1] = find(-sigma_t_values > 4e7,1,'first');
        % 
        % if isempty(ind) && isempty(ind1) % no boiling, no yielding
        %     fractional_thickening(i,j) = NaN;
        %     max_stress(i,j) = NaN;
        %     sdmax(i,j) = NaN;
        %     yield_coulomb(i,j) = Inf;
        %     yield_constant(i,j) = Inf;
        % 
        %     % for k=1:nthick
        %     %     Pex1 = Pex_values(k);
        %     %     rr = linspace(xi_values(k),p.R,100);
        %     %     srr = sigma_rr(Pex1,rr,ri_values(k),xi_values(k),p);
        %     %     stt = sigma_tt(Pex1,rr,ri_values(k),xi_values(k),p);
        %     %     % Pex,r,z,ri,xi,planet
        %     %     pp = 1/3*(srr+2*stt);
        %     %     plith = -rhoi*p.g*(p.R-rr); % lithostatic stress (compression negative)
        %     %     ptot = pp+plith;
        %     %     sd = srr-stt; % differential stress
        %     %     strength = -ptot*0.6 + 4e7;% 40 MPa + 0.6*P
        %     %     if any(sd > strength)
        %     %         fractional_thickening_to_yield(i,j) = -z_values(k)/(p.R-p.rc);
        %     %         break;
        %     %     end
        %     % end
        % elseif isempty(ind) && ~isempty(ind1)
        %     % yielding            
        %     fractional_thickening(i,j) = -z_values(ind1)/(p.R-p.rc);
        %     max_stress(i,j) = sigma_t_values(ind1);
        %     yield_constant(i,j) = true;
        %     sdmax(i,j) = sigma_t_values(ind1);
        % elseif ~isempty(ind) && isempty(ind1)
        %     % boiling
        %     fractional_thickening(i,j) = -z_values(ind)/(p.R-p.rc);
        %     % fractional_thickening_to_yield(i,j) = -z_values(ind)/(p.R-p.rc);
        %     max_stress(i,j) = sigma_t_values(ind);
        %     sdmax(i,j) = sigma_t_values(ind);
        % else % boiling AND yielding
        %     if ind < ind1
        %         % boiling occurs first
        %         fractional_thickening(i,j) = -z_values(ind)/(p.R-p.rc);
        %         max_stress(i,j) = sigma_t_values(ind);
        %         sdmax(i,j) = sigma_t_values(ind);
        %     else
        %         % yielding occurs first
        %         fractional_thickening(i,j) = -z_values(ind1)/(p.R-p.rc);
        %         max_stress(i,j) = sigma_t_values(ind1); 
        %         yield_constant(i,j) = true;
        %         sdmax(i,j) = sigma_t_values(ind1);
        %     end
        %     % % compute stresses in the elastic layer at the time of boiling
        %     % Pex1 = Pex_values(ind);
        %     % rr = linspace(xi_values(ind),p.R,100);
        %     % srr = sigma_rr(Pex1,rr,ri_values(ind),xi_values(ind),p);
        %     % stt = sigma_tt(Pex1,rr,ri_values(ind),xi_values(ind),p);
        %     %               % Pex,r,z,ri,xi,planet
        %     % pp = 1/3*(srr+2*stt);
        %     % plith = -rhoi*p.g*(p.R-rr); % lithostatic stress (compression negative)
        %     % ptot = pp+plith;
        %     % sd = srr-stt; % differential stress
        %     % strength = -ptot*0.6 + 4e7;% 40 MPa + 0.6*P
        %     % yield_coulomb(i,j) = any(sd > strength);
        %     % yield_constant(i,j) = any(sd > 4e7);
        %     % 
        %     % sdmax(i,j) = max(sd);
        % end        
    end
end

%% load table of satellites
TableofSatellites = import_satellite_table('Table_of_Satellites.xlsx', "Sheet1", [1, Inf]);
% This table assumes a core density of 3500 and ice+water density of 1000
satellite_f = 1-TableofSatellites.McMp;
satellite_R = TableofSatellites.RcKm./TableofSatellites.RcRp;
%% 
fractional_thickening_coulomb = min(fractional_thickening_boil,fractional_thickening_yield_coulomb);
boil_first_coulomb = fractional_thickening_boil <= fractional_thickening_yield_coulomb; % mask of boil before yield
boil_first_constant = fractional_thickening_boil <= fractional_thickening_yield_constant;
% This logic is necessary because comparisons involving NaN are always
% false
boil_region_coulomb = boil_first_coulomb | ...
    (isnan(fractional_thickening_yield_coulomb) & ~isnan(fractional_thickening_boil));
boil_region_constant = boil_first_constant | ...
    (isnan(fractional_thickening_yield_constant) & ~isnan(fractional_thickening_boil));
% yield_first_coulomb = fractional_thickening_yield_coulomb < fractional_thickening_boil;

figure
pcolor(RR,ff,gg);
shading flat
colorbar()

figure
% pcolor(RR,ff,h2o_thickness);
hold on
contour(RR,ff,h2o_thickness,25000:25000:1000000);
shading flat
colorbar()
title('h2o thickness')

f=figure(101);
% pos = get(gcf,'Position');
% f.Position(3) = pos(3)*2;
% subplot(1,2,1);
pcolor(RR/1e3,ff,fractional_thickening_coulomb);
shading flat;
hcb=colorbar();
hcb.Label.String = 'Fractional Thinning to Boil or Fail';
xlabel('Planetary Body Radius (km)');
ylabel('H_20 Mass Fraction (-)')
hold on
scatter(satellite_R,satellite_f,'ko','MarkerFaceColor','k');
text(satellite_R+25,satellite_f,TableofSatellites.Satellite);
contour(RR/1e3,ff,boil_region_coulomb,[0.5 0.5],'k','LineWidth',1);
contour(RR/1e3,ff,boil_region_constant,[0.5 0.5],'k--','LineWidth',1);
set(gca,'FontSize',14);
set(gca,'Layer','top');
caxis([0 1]);
set(gca,'YLim',[0 1]);
exportgraphics(gcf,['fractional_thinning' label '.pdf'],'Resolution',600);

figure
pcolor(RR/1e3,ff,yield_coulomb);
shading flat;
hold on
% contour(RR/1e3,ff,yield,[8 40]*10^6,'k');
hcb=colorbar();
% set(gca,'ColorScale','log');
hcb.Label.String = 'yield criterion?';
set(gca,'FontSize',14)
set(gca,'Layer','top')
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
