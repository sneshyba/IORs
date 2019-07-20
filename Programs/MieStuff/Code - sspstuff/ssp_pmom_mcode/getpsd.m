function f = getpsd(...
    sspfile, radiusfile, refffile, iorsfile, ppafile, ...
    LOGSG, Nforssp, ...
    NANG, NIORS, Pnrm, iops)

%   getpsd.m
%
%   Description:
%   getpsd calculates scattering parameters for a 
%   distribution of particle sizes
%
%   Input: 
%
%   Nstuff.txt      Number of phase function values, radii, etc
%   Pnrm.txt        Phase functions, Mie (indiv. radii)
%   iops.txt        inherent optical properties, Mie (qext, etc)
%   radius.txt      radii, Mie
%   reff.txt        radii, effective (desired)
%   iors_used.txt   frequencies and indices of refraction
%   ppa.txt         phase function angles
%
%
%   Output:
%
%   getpsd.log      A log of information on each distribution (QC)
%   ssp.txt         ssp file (phase functions, iops for range of r & nu)
%
%
%   Author: S. Neshyba, April 2001
%   Adapted to Matlab April 2009
%

% Parameters related to the width
nstddeviations = 4; % This says how far out to look for Mie candidates
%SG = exp(LOGSG);
kappa = exp(2.5*(LOGSG)^2);
LOGXDIFFMAX = LOGSG*nstddeviations; 

% Load list of un-averaged radii
radius_um = load(radiusfile);
radius_cm = radius_um/1e4;
Nradius = length(radius_um);

% Load list of desired effective radii & convert to geometric means
reff_um = load(refffile);
reff_cm = reff_um/1e4;
rbar_cm = reff_cm/kappa;
%rbar_um = rbar_cm*1e4;
Nreff = length(reff_um);

% Load list of indices of refraction; Niors had better equal NIORS
iors = load (iorsfile); 
nu_iors = iors(:,1); 
%m_iors = complex(iors(:,2),iors(:,3)); 
lambda_iors = 1e4./nu_iors;
%Niors = length(nu_iors);

% Load list of phase function angles
ppa = load(ppafile);

% Preallocate output arrays
Pnrm_psd = zeros(NANG,NIORS,Nreff); % zeros(NANG,1,1); % 
iops_psd = zeros(Nforssp,NIORS,Nreff); % zeros(Nforssp,1,1); %

% Open up a file to log some answers for quality control
fid = fopen('getpsd.log','w');

% Loop over frequencies
for i_iors = 1:NIORS
    
    % Report
    %['working on ... ' num2str(nu_iors(i_iors))]
    
    % Get size parameters for the radii
    XVAL = 2*pi*radius_cm*nu_iors(i_iors);
    AREA = pi*radius_um.^2;
    %VOL  = 4/3*pi*radius_um.^3;
    %XEFF = 2*pi*  reff_cm*nu_iors(i_iors);
    XBAR = 2*pi*  rbar_cm*nu_iors(i_iors);

    % Need some inherent optical properties too
    QEXT = iops(1,i_iors,:); QEXT = reshape(QEXT,Nradius,1);
    QSCA = iops(2,i_iors,:); QSCA = reshape(QSCA,Nradius,1);
    QABS = iops(3,i_iors,:); QABS = reshape(QABS,Nradius,1);
    %QBA  = iops(4,i_iors,:); QBA  = reshape(QBA, Nradius,1);
    ASYM = iops(5,i_iors,:); ASYM = reshape(ASYM,Nradius,1);
    %QRAT = iops(6,i_iors,:); QRAT = reshape(QRAT,Nradius,1);
    P11  = Pnrm(:,i_iors,:); P11  = reshape(P11,NANG,Nradius);
    CEXT = QEXT.*AREA;
    CSCA = QSCA.*AREA;
    CABS = QABS.*AREA;
    
    % Loop over desired radii
    for i_reff = 1:Nreff
        
        % Do the averaging
        [XBAR_PSD,X2_PSD,X3_PSD,XEFF_PSD,...
         CEXT_PSD,CSCA_PSD,CABS_PSD, ...
         QEXT_PSD,QSCA_PSD,QABS_PSD, ...
         ASYM_PSD, ...
         P11_PSD, ...
         ERROR,ICOUNT] = ...
                 psd(...
                     XBAR(i_reff),LOGSG,LOGXDIFFMAX,XVAL, ...
                     CEXT,CSCA,CABS, ...
                     QEXT,QSCA,QABS, ...
                     ASYM, ...
                     P11);

        
        % Some moments of r & associated means
        kappa_OUT = XEFF_PSD/XBAR_PSD;
        %R1 = XBAR_PSD/(2*pi)*lambda_iors(i_iors);
        R2 = X2_PSD^.5/(2*pi)*lambda_iors(i_iors);
        %R3 = X3_PSD^(1/3)/(2*pi)*lambda_iors(i_iors);
        AREA_PSD = pi*R2^2;
        %VOL_PSD  = 4/3*pi*R3^3;
        VOL_PSD = AREA_PSD *reff_um(i_reff)*4/3; % This should equal previous line
        %test1 = R2^2/R1^2;
        %test2 = R2^2/R3^(2/3);
        %est3 = R2^2/reff_um(i_reff)^2;
        test4 = (CEXT_PSD - QEXT_PSD*AREA_PSD)/CEXT_PSD*100;

        
        % Quality control
        fprintf(fid, ...
            '%8.3f %8.3f %5i %8.3f %8.3f %8.3f %8.3f \n', ...
            [nu_iors(i_iors) reff_um(i_reff) ICOUNT ERROR*1000 kappa kappa_OUT test4]);
        %plot(ppa,P11SUM,'o',ppa,P11)


        % Single scattering albedo
        SSA_PSD = QSCA_PSD/QEXT_PSD;


        % Repackage for saving
        Pnrm_psd(:,i_iors,i_reff)  = P11_PSD; % phase function
        iops_psd(1,i_iors,i_reff)  = lambda_iors(i_iors); % wavelength, um
        iops_psd(2,i_iors,i_reff)  = nu_iors(i_iors); % wavenumber, cm-1
        iops_psd(3,i_iors,i_reff)  = reff_um(i_reff); % effective radius, um
        iops_psd(4,i_iors,i_reff)  = CEXT_PSD; % extinction cross section
        iops_psd(5,i_iors,i_reff)  = CSCA_PSD; % scattering cross section
        iops_psd(6,i_iors,i_reff)  = CABS_PSD; % absorption cross section
        iops_psd(7,i_iors,i_reff)  = SSA_PSD; % single scattering albedo
        iops_psd(8,i_iors,i_reff)  = ASYM_PSD; % asymetry parameter
        iops_psd(9,i_iors,i_reff)  = QEXT_PSD; % extinction efficiency
        iops_psd(10,i_iors,i_reff) = QABS_PSD; % absorption efficiency
        iops_psd(11,i_iors,i_reff) = QSCA_PSD; % scattering efficiency
        iops_psd(12,i_iors,i_reff) = VOL_PSD; % volume, um^3
        iops_psd(13,i_iors,i_reff) = AREA_PSD; % projected area, um^2
    end
end
% Save
f = writesspfile (sspfile, NANG, NIORS, Nreff, ppa, Pnrm_psd, iops_psd);

% Close the log file
fclose(fid)





     