% Converts ssp file to moments file


% Read the sspfile
%sspfileroot = 'ssp_db.mie_wat.gamma_sigma_0p100';
sspfileroot = '../../../../Data/BL_hybrid_CRIs/hybrid_BL_273/ssp_getpsd_T273_S331';
sspfile = [sspfileroot '.txt']
pmomarrayfile = 'pmom.mat'
Npmomarrayfile = 'Npmom.mat'
ssp = readsspfileblind(sspfile)

% Some parameters
lmax = 70;
lmin = 8;
lbackoff = 4;
mu = cos(ssp.ppa/180*pi);

% One way to store the data
pmomarray = zeros(lmax,ssp.NIORS,ssp.Nreff);
Npmomarray = zeros(ssp.NIORS,ssp.Nreff);

% Loop over all the phase functions
for ireff = 1:ssp.Nreff
%for ireff = 10:11
    
    display(['working on reff = ', num2str(ssp.reff(ireff))]);
    
    for iIORS = 1:ssp.NIORS

        % Normalize the phase functions
        P_mu = phasenorm(mu,ssp.Pnrm_psd(:,ireff,iIORS),0);


        % Get the Legendre moments
        chilist = phase2pmom4(mu,P_mu,lmax,lmin,lbackoff);

        % Do some QC
%         figure(1)
%         subplot(1,2,1)
%         F = pmom2phase(chilist,mu);
%         plot(mu,P_mu,F.mu,F.P_mu,mu,F.P_mu-P_mu)
%         title('phase functions')
%         grid
% 
%         subplot(1,2,2)
%         semilogy(chilist,'-o')
%         title('moments')
%         xlim([0 lmax]);
%         ylim([.000001 1.2]);
%         grid
% 
%         pause(.02)

        if(chilist(end) > chilist(end-1))
            display(['Up-swing on ' num2str(ireff) ' ' num2str(iIORS)]);
        elseif(chilist(end) < 0)
            display(['Unphysical on ' num2str(ireff) ' ' num2str(iIORS)]);
        end
                   
        % Store
        pmomarray(1:length(chilist),iIORS,ireff) = chilist;
        Npmomarray(iIORS,ireff) = length(chilist);
        
        
    end
end

save ('pmomarrayfile') %,'pmomarray') syntax error;
save ('Npmomarrayfile') %, 'Npmomarray')