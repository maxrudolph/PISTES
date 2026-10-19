function [etaeff] = hirth_kohlstedt(s,T,d,P)
% s is the second invariant of the deviatoric stress, sqrt(0.5*tau_ij*tau_ij), in Pascal
% T is temperature in K
% d is grain size in m
% P is pressure in Pascal
% etaeff is effective viscosity in Pa-s, defined by tau_ij = 2*etaeff*edot_ij
%
% Hirth & Kohlstedt (2003) give axial strain rate as a function of
% differential stress. For uniaxial loading, s = sigma_d/sqrt(3) and
% edot_II = (sqrt(3)/2)*edot_axial, so each flow law is evaluated at
% sigma_d = sqrt(3)*s and converted back below (equivalent to the F2 factor
% of Gerya, Introduction to Numerical Geodynamic Modelling, for triaxial
% compression experiments).
assert(all(size(s) == size(T))); % ensure that sizes of s and T are compatible.
sd = sqrt(3)*s; % equivalent differential stress (Pa)
%
%
R = 8.3144621; % ideal gas constant, J/mol/K */
%   /* diffusion creep */
if( P < 1.0 )
    P = 1.0; %/* prevent spurious pressures from affecting viscosity */
end
% diffusion creep - dry olivine
A = 1.5e9; % grain size in micron, stress in MPa
n = 1;
p=3;
E=375000; %J/mol
V=5e-6;   %2-10 quoted in paper
ediff =  A*(sd/1e6).^n .* (d*1e6).^-p .* exp(-(E+P*V)/R./T);
ediff0 = A*(1/1e6)    .* (d*1e6).^-p .* exp(-(E+P*V)/R./T);

% dislocation creep - dry olivine
A = 1.1e5;
n = 3.5;
p = 0;
E = 530000;
V = 14e-6;
edisl = A*(sd/1e6).^n .* (d*1e6).^-p .* exp(-(E+P*V)/R./T);

etot = ediff + edisl ;%+ 1.0./(1.0./ebs + 1.0./egbs) ; % axial strain rate at sigma_d
%
etaeff = s./(sqrt(3)*etot); % = s/(2*edot_II), with edot_II = (sqrt(3)/2)*etot

mask = s < 1e0; %{/* if stress is really small (i.e. for first timestep), we cannot compute dislocation creep or gbs strain rates */
etaeff(mask) = 1.0./(3.0*ediff0(mask)); % Newtonian limit of the line above
