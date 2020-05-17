% makemieset
%
% Makes a set of mie parameters for a given range of wavelengths and
% particle sizes
%
% Author: Steven Neshyba, April 2009
%
% Input: 
%   iors.txt (indices of refraction)
%   ppa.txt (angles for phase function, degrees)
%
% Output:
%   radius.txt (desired radii, microns)
%   Nstuff (#angles in phase function, #iors, #radii, #in columns in iops)
%   Pnrm (normalized phase functions, ordered Nppa, Niors, Nradius) 
%   iops = (inherent optical properties,  ordered NfromMie, Niors, Nradius)
%   
%   
%   

% Report out where we are
!pwd

% Decide on radii
LOGSG = 0.1; % Assuming this is what getpsd uses
logsg_mieset = LOGSG/3; % This should be smaller than LOGSG used in getpsd
radius1 = .1;
radius2 = 200;
logradius_um = log10(radius1):logsg_mieset:log10(radius2);
radius_um = 10.^logradius_um';
Nradius = length(radius_um)


% Load up indices of refraction (iors) & repackage
%iors = load ('iors.txt');
iors = load ('Water_BL_298.txt');

% This line here to skip parts of the original ior dataset
%ior_interval = 20;
ior_interval = 1;
%range = 1:ior_interval:length(iors);
range = 1:ior_interval:4017;
iors_used = iors(range,:);
iors = iors_used;

% Make a record of the iors actually used
save iors_used.txt iors -ascii -double; 

% Carry on
m_iors = complex(iors(:,2),iors(:,3)); 
nu_iors = iors(:,1);
Niors = length(nu_iors);


% This is the number Mie.m returns
NfromMie = 6; 


% Load up the desired angles for phase function
ppa = load('ppa.txt');
u = cos(ppa*pi/180);
Nppa = length(u);
Nstuff = [Nppa Niors Nradius NfromMie]; 


% Save stuff
save 'radius.txt' radius_um -ascii -double;
save 'Nstuff.txt' Nstuff -ascii;


% Preallocate
Pnrm = zeros(Nppa,Niors,Nradius);
iops = zeros(NfromMie,Niors,Nradius);


% Loop over radii
for i_radius = 1:Nradius

    % Get the current radius and size parameter
    tic
    ['Working on r = ...' num2str(radius_um(i_radius)) ' (which is ' num2str(i_radius) ' out of ' num2str(Nradius) ')'] 
    radius_cm = radius_um(i_radius)/1e4; % microns converted to cm
    x_iors = 2*pi*radius_cm*nu_iors;

    % Do the Mie calculation for this set of iors
    for i_iors = 1:Niors
        %['Working on ior % = ...' num2str(i_iors)] 
        ftemp = Mie(m_iors(i_iors),x_iors(i_iors));
        qsca = ftemp(2);
        Ptemp = zeros(Nppa,1);
        for i_ppa=1:Nppa
            S12 = Mie_S12(m_iors(i_iors), x_iors(i_iors), u(i_ppa));
            Ptemp(i_ppa) = (abs(S12(1))^2 + abs(S12(2))^2);  
        end
        nrm = 2/(x_iors(i_iors)^2*qsca); % Normalization factor
        Ptemp_nrm = Ptemp*nrm;
        %trapz(u,Ptemp_nrm)
        Pnrm(:,i_iors,i_radius) = Ptemp_nrm;
        iops(:,i_iors,i_radius) = ftemp;
    end
    toc
end


% Save stuff
Pnrm_straight = reshape(Pnrm,Nppa*Niors*Nradius,1); 
save 'Pnrm.txt' Pnrm_straight -ascii -double;

iops_straight = reshape(iops,NfromMie*Niors*Nradius,1);
save 'iops.txt' iops_straight -ascii -double;



% Test the phase functions as saved and reconstructed
% [Nppa_test, Niors_test, Nradius_test, NfromMie_test, Pnrm_test, iops_test] = readmieset(...
%     'Nstuff.txt',...
%     'Pnrm.txt',...
%     'iops.txt');
%max(max(max(abs((iops_test-iops)./iops))))*100 %Should be zeros
%max(max(max(abs((Pnrm_test-Pnrm)./Pnrm))))*100 %Should be zeros




return


%     qext = zeros(Niors,1);
%     qsca = zeros(Niors,1);
%     qabs = zeros(Niors,1);
%     qb = zeros(Niors,1);
%     asy = zeros(Niors,1);
%     qratio = zeros(Niors,1);


% Plot the asymmetry parameter
% figure(1)
% plot(nu_iors,asy)
% xlabel('wavenumber (cm^{-1})')
% ylabel('g')

% Plot the phase function
% figure(2)
% plot(ppa,P)
% xlabel('Angle (degrees)')
% ylabel('Phase function')

