clear;
close all;

% moon = 'enceladus';

rhoi = 910;
rhow = 1000.0;
beta = 4e-10; % ocean compressibility Pa^-1
E = 5e9; % Young's modulus, Pa
nu = 0.33;

% Define data structures for each satellite
europa.name = 'Europa';
europa.R = 1560.8e3;
europa.ri = europa.R - 2.4e3;  % Manga and Wang 2007 value
europa.rc = europa.R - 1.20e5; % Manga and Wang, cites Cammarano et al. 2006
europa.g = 1.3;
europa.style = '--';

enceladus.name = 'Enceladus';
enceladus.R = 2.52e5; % outer radius
enceladus.ri = enceladus.R - 5e4; % Manga and Wang 2007 value.
enceladus.rc = 1.61e5;
enceladus.g = 0.113;
enceladus.style = '-';

mimas.name = 'Mimas';
mimas.R = 198.2e3;
mimas.rc = mimas.R-1.266e5;
mimas.ri = mimas.rc;
mimas.g = 0.064;
mimas.style = ':';

charon.name = 'Charon';
charon.R = 606e5;
charon.rc = 3.76e5;
charon.ri = charon.rc;
charon.g = 0.288;
charon.style = '-.';

% make the plot
nthick = 1000; % number of thicknesses

Pex = @(z_values,ri_values,xi_values,planet) (z_values .* (1-rhoi/rhow)) ./ (beta*(ri_values.^3-planet.rc^3)./(3*ri_values.^2) + ...
    xi_values/E.*(1 + 2*nu*(1+0.5*(planet.R./xi_values).^3)./( (planet.R./xi_values).^3-1) ) );
sigma_t = @(z_values,ri_values,xi_values,planet) 3/2*Pex(z_values,ri_values,xi_values,planet)./((planet.R./xi_values).^3-1);

% reproduce Manga and Wang 2007 results - 1/3 elastic fraction
p = enceladus;
z_values = logspace(1,4,nthick); % amount of thickening
ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
elastic_thickness = 1/3*(p.R-(ri_values));
xi_values = p.R-elastic_thickness;
Pex_13 = Pex(z_values,ri_values,xi_values,p);
sigmat_13 = sigma_t(z_values,ri_values,xi_values,p);
Pb_13 = rhoi*p.g*(p.R-ri_values) + Pex_13;

ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
elastic_thickness = 1e3; % fixed at 1 km
xi_values = p.R-elastic_thickness;
Pex_1km = Pex(z_values,ri_values,xi_values,p);
sigmat_1km = sigma_t(z_values,ri_values,xi_values,p);
Pb_1km = rhoi*p.g*(p.R-ri_values) + Pex_1km;

figure;
plot(z_values,Pex_13,'k--','DisplayName','Pex 1/3');
hold on
plot(z_values,sigmat_13,'k','DisplayName','\sigma_t 1/3');
plot(z_values,Pex_1km,'r--','DisplayName','Pex 1km');
plot(z_values,sigmat_1km,'r','DisplayName','\sigma_t 1km');

legend('Location','northwest');
set(gca,'XScale','log');
set(gca,'YScale','log');
title(p.name);

figure;
plot(z_values,Pb_13,'DisplayName','Pb 1/3')
hold on
plot(z_values,Pb_1km,'DisplayName','Pb 1 km')
set(gca,'XScale','log');
set(gca,'YScale','log');
legend()

%% Manga and Wang - Europa
p = europa;
z_values = logspace(1,4,nthick); % amount of thickening
ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
elastic_thickness = 1/3*(p.R-(ri_values));
xi_values = p.R-elastic_thickness;
Pex_13 = Pex(z_values,ri_values,xi_values,p);
sigmat_13 = sigma_t(z_values,ri_values,xi_values,p);
Pb_13 = rhoi*p.g*(p.R-ri_values) + Pex_13;

ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
elastic_thickness = 1e3; % fixed at 1 km
xi_values = p.R-elastic_thickness;
Pex_1km = Pex(z_values,ri_values,xi_values,p);
sigmat_1km = sigma_t(z_values,ri_values,xi_values,p);
Pb_1km = rhoi*p.g*(p.R-ri_values) + Pex_1km;

figure;
plot(z_values,Pex_13,'k--','DisplayName','Pex 1/3');
hold on
plot(z_values,sigmat_13,'k','DisplayName','\sigma_t 1/3');
plot(z_values,Pex_1km,'r--','DisplayName','Pex 1km');
plot(z_values,sigmat_1km,'r','DisplayName','\sigma_t 1km');

legend('Location','northwest');
set(gca,'XScale','log');
set(gca,'YScale','log');
title(p.name);

figure;
plot(z_values,Pb_13,'DisplayName','Pb 1/3')
hold on
plot(z_values,Pb_1km,'DisplayName','Pb 1 km')
set(gca,'XScale','log');
set(gca,'YScale','log');
legend(['total pressure ' p.name])


%% thinning ice shell...
z_values = -linspace(1,5e4,nthick); % amount of thickening

figure(1001); clf;
figure(1002); clf;
for p = [enceladus europa mimas charon]
    p.ri = p.rc;
    z_values = linspace(0,-((p.R-1e3-p.ri)),nthick); % amount of thickening - limited by ice shell thickness
    
    ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
    elastic_thickness = 1/3*(p.R-(ri_values));
    xi_values = p.R-elastic_thickness;
    Pex_13 = Pex(z_values,ri_values,xi_values,p);
    sigmat_13 = sigma_t(z_values,ri_values,xi_values,p);
    Pb_13 = rhoi*p.g*(p.R-ri_values) + Pex_13;

    % ri_values = p.ri-z_values; % ri should evolve with the changing ice shell thickness.
    % elastic_thickness = 1e3; % fixed at 1 km
    % xi_values = p.R-elastic_thickness;
    % Pex_1km = Pex(z_values,ri_values,xi_values,p);
    % sigmat_1km = sigma_t(z_values,ri_values,xi_values,p);
    % Pb_1km = rhoi*p.g*(p.R-ri_values) + Pex_1km;

    figure(1001);
    plot(z_values,Pex_13,'k','LineStyle',p.style,'DisplayName',[p.name ' 1/3']);
    hold on
    plot(z_values,sigmat_13,'k','DisplayName','\sigma_t 1/3');
    % plot(z_values,Pex_1km,'r--','DisplayName','Pex 1km');
    % plot(z_values,sigmat_1km,'r','DisplayName','\sigma_t 1km');

    figure(1002);
    Pb0 = p.g*rhoi*(p.R-p.rc);
    R0 = (p.R-p.rc);
    mask = Pb_13 > 600;
    plot(-z_values(mask)/R0,(Pb_13(mask)-(600))/(Pb0),'DisplayName',[p.name ' 1/3']);
    hold on
    % plot(z_values,Pb_1km,'DisplayName','Pb 1 km');
end
figure(1001);
legend('Location','northwest');
% set(gca,'XScale','log');
% set(gca,'YScale','log');
title(p.name);



% set(gca,'XScale','log');
% set(gca,'YScale','log');/

figure(1002);
    plot(get(gca,'XLim'),[0 0],'k')
    ylabel('(P-P_{tp})/P_c (-)')
xlabel('Fracional thinning (-)')
legend()
title('Basal Total Pressure')

% plot the excess pressure
% plot elastic stresses in the ice shell

% plot the pressure profile

