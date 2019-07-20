function f = writesspfile (sspfile, NANG, NIORS, Nreff, ppa, Pnrm, iops)
% Writes ssp-formatted file

% Calculate the number of lines
Nlines = round(NIORS*Nreff);

% Open the file and write a header
fid = fopen(sspfile,'w');
fprintf(fid,'Single scattering properties for water spheres\n');
fprintf(fid,'Size distribution is log-normal\n');
fprintf(fid,['     ' num2str(Nlines) ' data lines\n']);
fprintf(fid,['     ' num2str(NANG)   ' phase function angles\n']);
fprintf(fid,'phase function angles:\n');
fprintf(fid,'%13.5e',ppa);
fprintf(fid,'\n');
fprintf(fid,'wave [um]   wave [cm-1]     reff [um]    ext [um^2]   scat [um^2]    absx [um^2]        w0             g          Qext          Qabs          Qsca       Vol [um3]   Proj_area [um2]...and phase function values...\n');

% Decide on order to write out data as a function of r and nu
if NIORS > 1
    if iops(2,2,1) > iops(2,1,1) 
        i_iors_range = 1:NIORS;
    else
        i_iors_range = NIORS:-1:1;
    end
else
    i_iors_range = 1;
end

if Nreff > 1
    if iops(3,1,2) > iops(3,1,1)
        i_reff_range = 1:Nreff;
    else
        i_reff_range = Nreff:-1:1;
    end
else
    i_reff_range = 1;
end

% Write the phase function data  
for i_iors=i_iors_range
    for i_reff=i_reff_range
        
        % First set of data
        fline1 = iops(1:6,  i_iors,i_reff);
        fline2 = iops(7:11, i_iors,i_reff); 
        fline3 = iops(12:13,i_iors,i_reff);
        fprintf(fid,'%14.4f',fline1);
        fprintf(fid,'%14.6f',fline2);
        fprintf(fid,'%14.6e',fline3);
        fprintf(fid,'%13.6e',Pnrm(:,i_iors,i_reff));
        fprintf(fid,'\n');
    end
end
 
% Close and get out        
fclose(fid);
f=0;
return

