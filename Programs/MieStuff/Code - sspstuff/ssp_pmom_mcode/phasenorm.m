function f = phasenorm(mu,P_mu,debug)
% Normalizes phase function

Normfactor = -trapz(mu,P_mu)/2;
if (debug == 1) disp(['Normfactor = ' num2str(Normfactor)]); end
f = P_mu/Normfactor;

return