% runreadsspfile3 ... need to be in the Collaboration - Walden folder to
% run this

% This is Dave's data ... 
ftest1 = readsspfile2('iors_water_Turner/ssp_db.mie_wat.gamma_sigma_0p100.txt');

% Get the data at 273
ftest4 = readsspfile2('iors_water_Turner/ssp_water_273_lognormal_sigma_p487.txt');

% Get the data at 240
ftest5 = readsspfile2('iors_water_240/ssp_water_240_lognormal_sigma_p487.txt');

% Get the data at 240
ftest6 = readsspfile2('iors_water_240/ssp_vonsrpdf.txt');



iwantkabs = 0;
if (iwantkabs == 1)
    plot(...
        squeeze(ftest4.iops_psd(2,1,:)), squeeze(ftest4.iops_psd(6,1,:) ./ftest4.iops_psd(12,1,:)),'--o', ...
        squeeze(ftest5.iops_psd(2,1,:)), squeeze(ftest5.iops_psd(6,1,:) ./ftest5.iops_psd(12,1,:)),'--x', ...
        squeeze(ftest6.iops_psd(2,1,:)), squeeze(ftest6.iops_psd(6,1,:) ./ftest6.iops_psd(12,1,:)),'--+', ...
        squeeze(ftest1.iops_psd(2,14,:)),squeeze(ftest1.iops_psd(6,14,:)./ftest1.iops_psd(12,14,:)), ...
        'linewidth',1)
    set(gca,'fontsize',14)
    xlabel('wavenumber (cm^{-1})');
    ylabel('k_{abs}');
    legend('\sigma=0.487,T=273','\sigma=0.487,T=240','Von','Turner, T=273')
    grid
else

    % Index to show ...
    %indextoshow = 4; ylab='C_{ext} (\mum)';
    indextoshow = 9; ylab='Q_{ext}';
    %indextoshow = 6; ylab='C_{abs} (\mum)';
    %indextoshow = 12; ylab='Volume (\mum^3)';
    %indextoshow = 13; ylab='area (\mum^2)';

    plot(...
        squeeze(ftest4.iops_psd(2,1,:)), squeeze(ftest4.iops_psd(indextoshow,1,:)),'--*', ...
        squeeze(ftest5.iops_psd(2,1,:)), squeeze(ftest5.iops_psd(indextoshow,1,:)),'--*', ...
        squeeze(ftest6.iops_psd(2,1,:)), squeeze(ftest6.iops_psd(indextoshow,1,:)),'--*', ...
        squeeze(ftest1.iops_psd(2,14,:)),squeeze(ftest1.iops_psd(indextoshow,14,:)), ...
        'linewidth',1)
    set(gca,'fontsize',14)
    xlabel('wavenumber (cm^{-1}');
    ylabel(ylab);
    legend('\sigma=0.487,T=273','\sigma=0.487,T=240','Von','Turner, T=273')
end

