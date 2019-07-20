% Decide on a radius
radius_um = 10; % microns converted to cm
radius_cm = radius_um/1e4; % microns converted to cm

% Load some indices of refraction (iors) & repackage
water_xxx = load ('water_269.txt');
m_xxx = complex(water_xxx(:,2),water_xxx(:,2)); 
nu_xxx = water_xxx(:,1);
x_xxx = 2*pi*radius_cm*nu_xxx;
Nx = length(x_xxx);

% Do the Mie calculation for this set of iors
qext = zeros(Nx,1);
qsca = zeros(Nx,1);
qabs = zeros(Nx,1);
qb = zeros(Nx,1);
asy = zeros(Nx,1);
qratio = zeros(Nx,1);
for i = 1:Nx
    fi = Mie(m_xxx(i),x_xxx(i));
    qext(i) = fi(1);
    qsca(i) = fi(2);
    qabs(i) = fi(3);
    qba(i) = fi(4);
    asy(i) = fi(5);
    qratio(i) = fi(6);
end

% Plot the asymmetry parameter
figure(1)
plot(nu_xxx,asy)
xlabel('wavenumber (cm^{-1})')
ylabel('g')

% Choose a given wavelength and size parameter
m_index = 1;


% Load up the desired angles for phase function
load ppa.txt
u = cos(ppa*pi/180);
Nu = length(u);
P = zeros(Nu,1);
for i=1:Nu
    S12 = Mie_S12(m_xxx(m_index), x_xxx(m_index), u(i));
    P(i) = (abs(S12(1))^2 + abs(S12(2))^2);  
end
Nrm = 2/(x_xxx(m_index)^2*qsca(m_index)); % Normalization factor; Liou says 4/XX2QSCA
PNrm = P*Nrm;
trapz(u,PNrm)

% Plot the phase function
figure(2)
plot(ppa,P)
xlabel('Angle (degrees)')
ylabel('Phase function')

% Dump out an entry in Turner's format (approximately)
lambda_xxx = 1e4./nu_xxx; % in microns
reff = radius_um; % Setting effective radius to actual radius
geometric_area = pi*reff^2;
ext = qext*geometric_area;
scat = qsca*geometric_area;
absx = qabs*geometric_area;
w0 = qsca./qext;
Vol = 4/3*pi*reff^3;
Projarea = pi*reff^2;

fid = fopen('ssp_test.txt','w');
fprintf(fid,'Single scattering properties for wat spheres at such-and-such T\n');
fprintf(fid,'The size distribution is a delta-function (for now) \n');
fprintf(fid,'      ????  Number of data lines \n');
fprintf(fid,'       281  Number of angles in the phase function \n');
fprintf(fid,'phase function angles:\n');
fprintf(fid,'%13.5e',ppa);
fprintf(fid,'\n');
fprintf(fid,'wave [um]   wave [cm-1]     reff [um]    ext [um^2]   scat [um^2]    absx [um^2]        w0             g          Qext          Qabs          Qsca       Vol [um3]   Proj_area [um2]...and phase function values...\n');
fline1 = [...
        lambda_xxx(m_index) ...
        nu_xxx(m_index) ...
        reff(m_index) ...
        ext(m_index) ...
        scat(m_index) ...
        absx(m_index)]; 
fline2 = [... 
        w0(m_index) ...
        asy(m_index) ...
        qext(m_index) ...
        qabs(m_index) ...
        qsca(m_index)];
fline3 = [...
        Vol ...
        Projarea];

fprintf(fid,'%14.4f',fline1);
fprintf(fid,'%14.6f',fline2);
fprintf(fid,'%14.6e',fline3);
fprintf(fid,'%13.6e',PNrm);
fprintf(fid,'%\n');
fclose(fid);

