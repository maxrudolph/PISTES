function [H] = mars_heating(t)
if t>4.6e9
    warning('Time is greater than 4.5 Ga');
end
% time should be specified in years
% the heating model implements Table 1 of Thiriet et al. 2019, PEPI.
% this is a model for the primitive mantle of mars.
H0 = 23e-12; %pW/kg
% H4p5 = 4;%pW/kg
% assume H = H0*exp(-lam_eff*t)
% H(4.5) = H0*exp(-lam_eff*4.5)
% log(H4p5/H0) = -lam_eff*4.5
% lam_eff = -log(H4p5)/log(H0)/4.5e9;
lam_eff = 3.8871e-10; % year^-1
H = H0*exp(-lam_eff*t);
end

