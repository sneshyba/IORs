% runreadsspfile4 ... 

% Get the data at 240
fvonxx = readsspfile2('ssp_vonsrpdf.txt');

% Get the data at 240
ftest2 = readsspfile2('ssp_method2.txt');

% Get the data at 240
ftest4 = readsspfile2('ssp_method4.txt');

% Get the data at 240
ftest5 = readsspfile2('ssp_method5.txt');


ftest6 = readsspfile2('ssp_water_240_lognormal_sigma_p331.txt');


iwantkabs = 1;
if (iwantkabs == 0)
    plot(...
        squeeze(fvonxx.iops_psd(2,1,:)),squeeze(fvonxx.iops_psd(6,1,:)./fvonxx.iops_psd(12,1,:)), ...
        squeeze(ftest2.iops_psd(2,1,:)),squeeze(ftest2.iops_psd(6,1,:)./ftest2.iops_psd(12,1,:)),'--*', ...
        squeeze(ftest4.iops_psd(2,1,:)),squeeze(ftest4.iops_psd(6,1,:)./ftest4.iops_psd(12,1,:)),'--+', ...
        squeeze(ftest5.iops_psd(2,1,:)),squeeze(ftest5.iops_psd(6,1,:)./ftest5.iops_psd(12,1,:)),'--o', ...
        squeeze(ftest6.iops_psd(2,1,:)),squeeze(ftest6.iops_psd(6,1,:)./ftest6.iops_psd(12,1,:)),'--o', ...
        'linewidth',1)
    set(gca,'fontsize',14)
    xlabel('wavenumber (cm^{-1})');
    ylabel('k_{abs}');
    legend('vonsrpdf','method 2: 0.487','method 4: 0.1','method 5: 0.331','getpsd: 0.331')
    grid
else

    % Index to show ...
    %indextoshow = 4; ylab='C_{ext} (\mum)';
    indextoshow = 9; ylab='Q_{ext}';
    %indextoshow = 6; ylab='C_{abs} (\mum)';
    %indextoshow = 12; ylab='Volume (\mum^3)';
    %indextoshow = 13; ylab='area (\mum^2)';

    plot(...
        squeeze(fvonxx.iops_psd(2,1,:)), squeeze(fvonxx.iops_psd(indextoshow,1,:)), ...
        squeeze(ftest2.iops_psd(2,1,:)), squeeze(ftest2.iops_psd(indextoshow,1,:)),'--*', ...
        squeeze(ftest4.iops_psd(2,1,:)), squeeze(ftest4.iops_psd(indextoshow,1,:)),'--+', ...
        squeeze(ftest5.iops_psd(2,1,:)), squeeze(ftest5.iops_psd(indextoshow,1,:)),'--o', ...
        squeeze(ftest6.iops_psd(2,1,:)), squeeze(ftest6.iops_psd(indextoshow,1,:)),'--o', ...
        'linewidth',1)
    set(gca,'fontsize',14)
    xlabel('wavenumber (cm^{-1}');
    ylabel(ylab);
    legend('vonsrpdf','method 2','method 4','method 5','getpsd')
end

