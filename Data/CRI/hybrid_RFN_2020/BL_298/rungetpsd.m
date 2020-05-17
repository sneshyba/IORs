% Rungetpsd makes an ssp file


% Get the un-averaged Mie dataset
alreadyread = 0
if alreadyread == 0
    'Reading the mieset ...'
    tic
    [NANG, NIORS, NXVAL, NfromMie, Pnrm, iops] = readmieset(...
    'Nstuff.txt',...
    'Pnrm.txt',...
    'iops.txt');
    Nforssp = 13;
    radiusfile = 'radius.txt';
    radius_um = load(radiusfile);
    dlogsg = log(radius_um(2)/radius_um(1))
    toc
end

% Get the distributions
f = getpsd(...
    'ssp_getpsd_BL298_S487.txt', radiusfile, 'reff.txt', 'iors_used.txt', 'ppa.txt', ...
    0.487, Nforssp, NANG, NIORS, Pnrm, iops)
f = getpsd(...
    'ssp_getpsd_BL298_S331.txt', radiusfile, 'reff.txt', 'iors_used.txt', 'ppa.txt', ...
    0.331, Nforssp, NANG, NIORS, Pnrm, iops)
