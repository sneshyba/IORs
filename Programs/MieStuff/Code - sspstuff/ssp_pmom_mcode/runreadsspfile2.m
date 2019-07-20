% runreadsspfile2

% This is Dave's data ... 
ftest1 = readsspfile2('ssp_db.mie_wat.gamma_sigma_0p100.txt');

% Get the data
ftest2 = readsspfile2('ssp_water_273_lognormal_sigma_p100.txt');

% Get the data
ftest3 = readsspfile2('ssp_water_273_lognormal_sigma_p331.txt');

% Get the data
ftest4 = readsspfile2('ssp_water_273_lognormal_sigma_p487.txt');

iwantkabs = 1;
if (iwantkabs == 1)
    plot(...
        squeeze(ftest2.iops_psd(2,1,:)), squeeze(ftest2.iops_psd(6,1,:) ./ftest2.iops_psd(12,1,:)),'--o', ...
        squeeze(ftest3.iops_psd(2,1,:)), squeeze(ftest3.iops_psd(6,1,:) ./ftest3.iops_psd(12,1,:)),'--x', ...
        squeeze(ftest4.iops_psd(2,1,:)), squeeze(ftest4.iops_psd(6,1,:) ./ftest4.iops_psd(12,1,:)),'--*', ...
        squeeze(ftest1.iops_psd(2,14,:)),squeeze(ftest1.iops_psd(6,14,:)./ftest1.iops_psd(12,14,:)), ...
        'linewidth',1)
    set(gca,'fontsize',14)
    xlabel('wavenumber (cm^{-1})');
    ylabel('k_{abs}');
    legend('\sigma=0.1','\sigma=0.331','\sigma=0.487','Turner')

else

    % Index to show ...
    %indextoshow = 4; ylab='C_{ext} (\mum)';
    %indextoshow = 9; ylab='Q_{ext}';
    indextoshow = 6; ylab='C_{abs} (\mum)';
    %indextoshow = 12; ylab='Volume (\mum^3)';
    %indextoshow = 13; ylab='area (\mum^2)';

    plot(...
        squeeze(ftest2.iops_psd(2,1,:)), squeeze(ftest2.iops_psd(indextoshow,1,:)),'--o', ...
        squeeze(ftest3.iops_psd(2,1,:)), squeeze(ftest3.iops_psd(indextoshow,1,:)),'--x', ...
        squeeze(ftest4.iops_psd(2,1,:)), squeeze(ftest4.iops_psd(indextoshow,1,:)),'--*', ...
        squeeze(ftest1.iops_psd(2,14,:)),squeeze(ftest1.iops_psd(indextoshow,14,:)), ...
        'linewidth',1)
    set(gca,'fontsize',14)
    xlabel('wavenumber (cm^{-1}');
    ylabel(ylab);
    legend('\sigma=0.1','\sigma=0.331','\sigma=0.487','Turner')
end

