function [Nppa_test, Niors_test, Nradius_test, NfromMie_test, Pnrm_test, iops_test] = readmieset(f1,f2,f3)
% Gets iops from files

% Load the phase functions etc, and reconstruct
Nstuff = load(f1);
Nppa_test = Nstuff(1);
Niors_test = Nstuff(2);
Nradius_test = Nstuff(3);
NfromMie_test = Nstuff(4);
Pnrm_straight = load(f2);
iops_straight = load(f3);
Pnrm_test = reshape(Pnrm_straight,Nppa_test,Niors_test,Nradius_test);
iops_test = reshape(iops_straight,NfromMie_test,Niors_test,Nradius_test);
return
end


