clear;
% close all;


% mars-like parameters
seconds_in_year = 3.15e7;
p = struct();
p.rhom = 3500;
p.Cpm  = 1142;
p.Ro = 3.3895e6;      % outer radius of lithosphere
p.Rc = 1.830e6;       % core radius, m (Samuel et al., 2023)
p.alpha = 2.5e-5;
p.g = 3.73;           % surface gravity (m/s^2)
p.k = 4.4;
p.kappa = p.k/p.rhom/p.Cpm;
p.Ts = 270;            % Surface temperature, K.
p.crustal_heating_factor = 10.7108;
p.mantle_heating_factor  = 0.4273;
p.hc = 60e3; % crustal thickness (m)
mub = 1e21;
Q = 300;
R = 8.314e-3;
Tref=1600;

eta = @(T,stress) mub*exp(Q/R*(1./T - 1./Tref)); % Michaut et al. 2025 - Arrhenius form
dTnu = @(T) R/Q*T^2; % rheological temperature scale (positive sign??)
p.eta = eta;
p.dTnu = dTnu;

Tm0 = 1700;
D0  = 60000;

lid_rhs(0.0,[Tm0 D0],p)

options = odeset();
options.InitialStep = seconds_in_year*1e5;
[t,y] = ode45(@(t,y) lid_rhs(t,y,p),[0 4.5e9*seconds_in_year],[Tm0,D0],options);
figure();
subplot(2,1,1);
plot(t/seconds_in_year,y(:,1));
subplot(2,1,2);
plot(t/seconds_in_year,y(:,2));

function dydt = lid_rhs(t,y,p)
        seconds_in_year = 3.15e7;

        Tm = y(1);  % temperature at base of lid
        D  = y(2);  % lid thickness
        arh =2.0;   % constant from Michaut equation 12 
        C = 0.47;   % Davaille and Jaupart 1993
        mantle_volume = 4/3*pi*((p.Ro-D)^3-p.Rc^3);
        Cm = p.rhom*p.Cpm*mantle_volume;
        Slid = 4*pi*(p.Ro-D)^2;      
        
        h_conv = p.mantle_heating_factor*p.rhom*mars_heating(t/seconds_in_year); % mantle volumetric heat production
        
        qbl = C*p.k*(p.alpha*p.rhom*p.g/p.kappa/p.eta(Tm,0))^(1/3)*p.dTnu(Tm)^(4/3);
        DTbl = arh*p.dTnu(Tm);
        % for now make a really naive assumption about heat conduction
        % through the lid... revise later!
        % qlid = p.k*(Tm-DTbl)/D; % Tm-dTbl is the temperature at the base of the conductive lid.
        
        % steady state conduction solution for a lid with heating
        % only above depth D (ignore heating below crust)
        Tl = Tm-DTbl;% temp at base of conductive layer
        Hc = p.crustal_heating_factor*mars_heating(t/seconds_in_year)*p.rhom;
        Hm = p.mantle_heating_factor*mars_heating(t/seconds_in_year)*p.rhom;
        qlid = steady_conduction_sphere_heating(Tl,Hc,Hm,D,p);
        % a=1/D*(Tl-p.Ts+Hc/2/p.k*hc^2+Hc*hc/p.k*(D-hc));
        % dTdz = a;
        % qlid = p.k*dTdz;

        dH = p.rhom*p.Cpm*DTbl; % enthalpy change across the lid
        dTmdt = 1/Cm * (-Slid*qlid + h_conv*mantle_volume); % Michaut et al. Equation 18
        dDdt = 1/dH * (qlid-qbl);   % Michaut et al. Equation 19
        dydt=[dTmdt dDdt]';
end