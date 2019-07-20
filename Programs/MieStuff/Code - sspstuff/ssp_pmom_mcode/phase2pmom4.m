function f = phase2pmom4(mu,P_mu,lmax,lmin,lbackoff)
% Converts phase function to legendre moments
nrm = abs(trapz(mu,P_mu));
lasterr = nrm;


chilist = zeros(lmax+1,1); betalist = chilist; 
for l = 0:lmax
    test = -.5*trapz(mu,mylegendreP(l,0,mu).*P_mu);
    chilist(l+1) = test;
    Ftest = pmom2phase(chilist(1:l+1),mu);
    error = std(Ftest.P_mu-P_mu)/nrm;
    %[error lasterr]
    if (l > lmin & error > lasterr)
        break
    end
    lasterr = error;
end
indx = max(l-lbackoff,lmin);
f = chilist(1:1+indx,1);
return
