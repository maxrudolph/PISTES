function [z,Hprob,Dprob,Tprob]=ProbaQuakeDepth()
% clear all
% close all
run("Catalogue_Depth.m")

%Drilleau
SMAX=size(S);
sum=@(z)0;

for i=1:SMAX(1)
h=@(z)1./S(i,8)/(2*pi)^0.5.*exp(-(z-S(i,7)).^2./(2*S(i,8)^2));
hint=integral(h,0,200);
sum=@(z)sum(z)+h(z)./hint;
end
Hsum=integral(sum,0,200);
z=0:0.1:200;
% plot(sum(z)/Hsum,-z,'b')
% hold on
Hprob = sum(z)/Hsum;

%Duran
DMAX=size(D);
dsum=@(z)0;

for i=1:DMAX(1)
hd=@(z)1./D(i,4)/(2*pi)^0.5.*exp(-(z-D(i,3)).^2./(2*D(i,4)^2));
hdint=integral(hd,0,200);
dsum=@(z)dsum(z)+hd(z)./hdint;
end
HDsum=integral(dsum,0,200);
z=0:0.1:200;
% plot(dsum(z)/HDsum,-z,'k--')
Dprob = dsum(z)/HDsum;

% Stahler
TMAX=size(T);
tsum=@(z)0;

for i=1:TMAX(1)
ht=@(z)1./T(i,2)/(2*pi)^0.5.*exp(-(z-T(i,1)).^2./(2*T(i,2)^2));
htint=integral(ht,0,200);
tsum=@(z)tsum(z)+ht(z)./htint;
end
HTsum=integral(tsum,0,200);
z=0:0.1:200;
% plot(tsum(z)/HTsum,-z,'g-.')
Tprob = tsum(z)/HTsum;

% ylabel('Depth (km)')
% xlabel('Density Probability for Quake Depths')
