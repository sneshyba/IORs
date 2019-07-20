% runreadsspfile

% Get the data
[NANG,Niors,Nreff,ppa,Pnrm_psd,iops_psd,reff,iors] = readsspfile('ssp.txt','reff.txt','iors_used.txt');
Nreff = length(reff);
reff_cm = reff/1e4;
Niors = length(iors);
m_iors = complex(iors(:,2),iors(:,3));
nu_iors = iops_psd(2,:,1)';

% Testing ... compare to fig. 3.12 of Goody
for ilambda = 1:Niors;
x=2*pi*reff_cm*nu_iors(ilambda);
rho=2*x*abs(m_iors(ilambda)-1);

iwant = 3;
if (iwant == 1)
    semilogx(...
        rho,reshape(iops_psd(9, ilambda,:),Nreff,1),...
        rho,reshape(iops_psd(10,ilambda,:),Nreff,1),...
        rho,reshape(iops_psd(11,ilambda,:),Nreff,1));
    xlabel('\rho');
    xlim([0.1 1e4]);
    ylim([0 3.5]);
    grid;legend('e','s','a');
elseif (iwant ==2)
    semilogx(...
        reff,reshape(iops_psd(9, ilambda,:),Nreff,1),...
        reff,reshape(iops_psd(10,ilambda,:),Nreff,1),...
        reff,reshape(iops_psd(11,ilambda,:),Nreff,1));
    xlabel('r (\mum');
    xlim([min(reff) max(reff)]);
    ylim([0 3.5]);
    grid;legend('e','s','a');
elseif (iwant ==3)
    semilogx(...
        reff,reshape(iops_psd(8,ilambda,:),Nreff,1));
    xlabel('r (\mum)');
    xlim([min(reff) max(reff)]);
    ylim([0 1]);
end
title(num2str(nu_iors(ilambda)))
pause(.1)
end

