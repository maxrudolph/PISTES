function [etaeff] = goldsby_kohlstedt(s,T,d,P)
% s is the second invariant of the deviatoric stress, sqrt(0.5*tau_ij*tau_ij), in Pascal
% T is temperature in K
% d is grain size in m
% P is pressure in Pascal
% etaeff is effective viscosity in Pa-s, defined by tau_ij = 2*etaeff*edot_ij
%
% Goldsby & Kohlstedt (2001) give axial strain rate as a function of
% differential stress (uniaxial experiments). For uniaxial loading,
% s = sigma_d/sqrt(3) and edot_II = (sqrt(3)/2)*edot_axial, so each flow law
% is evaluated at sigma_d = sqrt(3)*s and converted back below
% (equivalent to the F2 factor of Gerya, Introduction to Numerical
% Geodynamic Modelling, for triaxial compression experiments).
assert(all(size(s) == size(T))); % ensure that sizes of s and T are compatible.
sd = sqrt(3)*s; % equivalent differential stress (Pa)
%
%
% PetscScalar goldsbyKohlstedt(PetscScalar s, PetscScalar T, PetscScalar d, PetscScalar P){
%   /* returns strain rate predicted by goldsby and kohlstedt 2001 in units of s^-1*/
R = 8.3144621; % ideal gas constant, J/mol/K */
%   /* diffusion creep */

P = max(P,1.0); %/* prevent spurious pressures from affecting viscosity (element-wise if P is an array) */
%   {
%     /* quantities from table 6: Diffusion Creep Parameters */
b = 4.52e-10; %/* burgers vector */
Vm = 1.97e-5; %/* m^3 */
Qv = 59.4*1000.0;%/* Activation energy for volume diffusion, J/mol */
delta = 2.0*b; %/* grain boundary width, m */
Qb = 49*1000.0;
D0v = 9.10e-4; %/* m^2/s */
D0b = 5.8e-4; %/* m^2/s - this is lower of two upper bounde presented in text */
%     /* diffusivities */
Dv = D0v*exp(-Qv./(R*T));
Db = D0b*exp(-Qb./(R*T));
%     /* note that all units above are SI, so strain rate will be in 1/s for stress in Pa */
ediff0 = 42.0*Vm./(R*T*d*d).*(Dv + pi*delta/d*Db );

ediff = ediff0 .* sd;
%   {
%     PetscScalar A,n,Q;
egbs = zeros(size(T));
V = -13e-6; %/* m^3/mol, activation volume from text */
p = 1.4; %/* grain size dependence */
mask = T<255;

A=3.9e-3;% /* MPa^-1.8 m^1.4 s^-1 */
n=1.8;
Q=49*1000.0; %/* J/mol */
egbs(mask) = A*(sd(mask)/1e6).^n/d^p.*exp( (-Q + P*V)./(R*T(mask)) );
%     }else{
A = 3.0e26; %/* MPa^-1.8 m^1.4 S^-1 */
n = 1.8;
Q = 192*1000.0;
egbs(~mask) = A*(sd(~mask)/1e6).^n/d^p.*exp( (-Q + P*V)./(R*T(~mask)) );

%   PetscScalar edisl;
%   {
%     PetscScalar A,n,Q;
% T is a vector, so choose parameters node by node (an 'if' on an array is
% true only if every element is true).
warm = T >= 258;
A = 4.0e5*ones(size(T));% /* MPa^-4.0 s^-1, T < 258 K */
Q = 60*1000.0*ones(size(T)); %/* J/mol, T < 258 K */
A(warm) = 6.0e28; %/* MPa^-4.0 S^-1 */
Q(warm) = 180*1000.0; %/* Note that there is a mistake in G&K table 5. Should be 180 kJ/mol per text */
n = 4.0;
V = -13e-6; %/* m^3/mol, activation volume from text */
edisl = A.*(sd/1e6).^n.*exp( (-Q + P*V)./(R*T) );
%   }
%   PetscScalar ebs;
%   {
A = 5.5e7; %/* MPa^-2.4 s^-1 */
n = 2.4;
Q = 60.0*1000.0; %/* J/mol */
V = -13e-6; %/* m^3/mol, activation volume from text */
ebs = A*(sd/1e6).^n.*exp( (-Q + P*V)./(R*T) );
%   }
etot = ediff + 1.0./(1.0./ebs + 1.0./egbs) + edisl; % axial strain rate at sigma_d
%
etaeff = s./(sqrt(3)*etot); % = s/(2*edot_II), with edot_II = (sqrt(3)/2)*etot

mask = s < 1e0; %{/* if stress is really small (i.e. for first timestep), we cannot compute dislocation creep or gbs strain rates */
etaeff(mask) = 1.0./(3.0*ediff0(mask)); % Newtonian limit of the line above
