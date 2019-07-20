function f = pmom2phase(chilist,mu)
% Converts legendre moments to phase function

lmax = length(chilist)-1;
P_mu = zeros(size(mu));
for l = 0:lmax
    P_mu = P_mu + (2*l+1)*chilist(l+1)*mylegendreP(l,0,mu); 
end
f.mu = mu;
f.P_mu = P_mu;
return
