function [qbot] = steady_conduction_sphere_heating(Tl,Hu,Hl,D,p)
% solve the steady heat equation subject to internal heating that varies
% within two layers
% Tl is temperature at base of layer (lid)
% Hu is crust (upper layer) heating
% Hl is mantle (lower layer) heating

% 
% Tl = 1700;
% p.hc = 60e3;
% p.D = 100e3;
% p.Ro = 3.3895e6;      % outer radius of lithosphere
% p.Ts = 250;            % Surface temperature, K.
% p.k = 4.4;
Ri = p.Ro - D;
r1 = p.Ro-p.hc;
k=p.k;
% 
% Hu = 1e-6;
% Hl = 3e-6*0;
% The solution is piecewise constant and has the form
% T(r) = -H/k r^2/6+C0/r+C1
% lower layer:
% set up equations such that
% [L]*[C0u C1u C0l C1l] = [R]
L = zeros(4,4); R = zeros(4,1);
% T(R0) = Ts:
L(1,:) = [-1/p.Ro 1 0 0]; R(1) = p.Ts + Hu/k*p.Ro^2/6;
% T(R0-D) = Tl:
L(2,:) = [0 0 -1/Ri 1];   R(2) = Tl + Hl/k*Ri^2/6;
% equal heat fluxes across hc:
L(3,:) = [1/r1^2 0 -1/r1^2 0]; R(3) = -Hl*r1/k/3 + Hu*r1/k/3;
% equal temperature at h
L(4,:) = [-1/r1 1 1/r1 -1]; R(4) = -Hl/k*r1^2/6 + Hu/k*r1^2/6;

S = L\R;

qbot = -k*(-Hl/k*Ri/3   + S(3)/Ri^2)


doplot=0;
if doplot
    nr = 501;
    r = linspace(Ri,p.Ro,nr);
    T = zeros(size(r));
    mask = r>r1;
    T(mask)  = -Hu/k*r(mask).^2/6 - S(1)./r(mask) + S(2);
    T(~mask) = -Hl/k*r(~mask).^2/6 - S(3)./r(~mask) + S(4);
    figure()
    plot(r,T);
end
% Debugging stuff: 
% Basal heat flux
% qbot = -k*(-Hl/k*Ri/3   + S(3)/Ri^2)
% qtop = -k*(-Hu/k*p.Ro/3 + S(1)/p.Ro^2)
% check energy conservation
% Htop = 4/3*pi*(p.Ro^3-r1^3)*Hu;
% Hbtm = 4/3*pi*(r1^3-Ri^3)*Hl;
% Atop = 4*pi*p.Ro^2;
% Abtm = 4*pi*Ri^2;
% qbot*Abtm + Htop+Hbtm - qtop*Atop
return