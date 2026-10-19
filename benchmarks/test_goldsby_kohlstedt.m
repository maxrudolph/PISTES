% test the goldsby kohlstedt rheology.
% I used this code to reproduce the curves in Figure 7
% of Goldsby and Kohlstedt (2001).
clear
close all

addpath ../core/

sigma = logspace(-4,2)*1e6;
T = 258* ones(size(sigma));
d=1e-3;
P=1.3*1000*1e4;

% sigma is the differential stress of a uniaxial experiment. goldsby_kohlstedt
% takes the deviatoric stress invariant, sigma/sqrt(3), and the axial strain
% rate is then (2/3*sigma)/(2*eta) = sigma/(3*eta).
eta = goldsby_kohlstedt(sigma/sqrt(3),T,d,P);

figure()
plot(sigma/1e6,sigma./(3*eta));
xlabel('stress (MPa)');
ylabel('strain rate (1/s)');
set(gca,'YLim',[1e-12 1e0]);
set(gca,'XLim',[1e-4 1e2]);
set(gca,'XScale','log');
set(gca,'YScale','log');

% mixed-temperature check: a profile that straddles the 258 K
% dislocation-creep transition, evaluated in one call, must match
% evaluating each node on its own.
Tprof = linspace(100,273,60);
for s0 = [0.1 1 5]*1e6
    sprof = s0*ones(size(Tprof));
    eta_vec = goldsby_kohlstedt(sprof,Tprof,d,P);
    eta_node = arrayfun(@(i) goldsby_kohlstedt(sprof(i),Tprof(i),d,P),1:numel(Tprof));
    assert(max(abs(eta_vec-eta_node)./eta_node) < 1e-12, ...
        'goldsby_kohlstedt: vector evaluation differs from per-node evaluation at %g MPa',s0/1e6);
end
disp('goldsby_kohlstedt mixed-temperature check passed');

