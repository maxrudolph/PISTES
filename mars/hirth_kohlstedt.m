function [etaeff] = hirth_kohlstedt(s,T,d,P)
% s is stress in Pascal
% T is temperature in K
% d is grain size in m
% P is pressure in Pascal
% etaeff is effective viscosity in Pa-s
assert(all(size(s) == size(T))); % ensure that sizes of s and T are compatible.
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
ediff =  A*(s/1e6).^n .* (d*1e6).^-p .* exp(-(E+P*V)/R./T);
ediff0 = A*(1/1e6)    .* (d*1e6).^-p .* exp(-(E+P*V)/R./T);

% dislocation creep - dry olivine
A = 1.1e5;
n = 3.5;
p = 0;
E = 530000;
V = 14e-6;
edisl = A*(s/1e6).^n .* (d*1e6).^-p .* exp(-(E+P*V)/R./T);

etot = ediff + edisl ;%+ 1.0./(1.0./ebs + 1.0./egbs) ;
%
etaeff = s./(2.0*etot);

mask = s < 1e0; %{/* if stress is really small (i.e. for first timestep), we cannot compute dislocation creep or gbs strain rates */
etaeff(mask) = 0.5./ediff0(mask);
