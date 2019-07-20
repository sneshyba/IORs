function ...
[XVALSUM,X2SUM,X3SUM,XCSUM,...
CEXTSUM,CSCASUM,CABSSUM, ...
QEXTSUM,QSCASUM,QABSSUM, ...
ASYMSUM, ...
P11SUM, ...
ERROR,ICOUNT] = ...
         psd(...
             XBAR,LOGSG,LOGXDIFFMAX,XVAL, ...
             CEXT,CSCA,CABS, ...
             QEXT,QSCA,QABS, ...
             ASYM, ...
             P11)
%
%   Description:
%   psd calculates scattering parameters for a 
%   distribution of particle sizes
%
%
%   Author: S. Neshyba, April 2001
%   Converted to Matlab April 2009
%

% Parameters frequently used locals
SMALL = .008;
IDEBUG = 2;
R2PI = (2*pi)^.5;
NXVAL = length(XVAL);
NANG = length(P11);
%IFLAG = 0;

% Log of x-interval stuff
dLOGX = log(XVAL(2)/XVAL(1));
LOGXBAR = log(XBAR);

% Initialize sums...
ICOUNT = 0;
LOGXVALSUM = 0;
XVALSUM = 0;
SUMN  = 0;
SUMNS = 0;
QEXTSUM = 0;
QSCASUM = 0;
QABSSUM = 0;
CEXTSUM = 0;
CSCASUM = 0;
CABSSUM = 0;
ASYMSUM = 0;
XCSUM = 0;
X2SUM = 0;
X3SUM = 0;
P11SUM = zeros(NANG,1);
ERROR = 0;


% Loop over the input size parameters around xbar
for I = 1:NXVAL
    XVALI = XVAL(I);
    LOGXVALI = log(XVALI);
    LOGXDIFF = LOGXVALI-LOGXBAR;
    if (abs(LOGXDIFF) < LOGXDIFFMAX)

        
    %   Get the weights and the volume element
    %   Note: FNdN is Grenfell's PROD01, FNXSdN is Grenfell's PROD1...
        ICOUNT = ICOUNT+1;
        ARGFN = LOGXDIFF^2/(2*LOGSG^2);
        EXPARGFN = exp(-ARGFN);
        PREFACTOR = 1/(LOGSG*XVALI*R2PI);
        FN = PREFACTOR * EXPARGFN;
        dN = XVALI * dLOGX;
        FNdN   = FN * dN;
        FNXSdN = FN * dN * XVALI^2;


    %   Normalization sums...
        SUMN   = SUMN  + FNdN;
        SUMNS  = SUMNS + FNXSdN;


    %   Accumulate first moments...
        LOGXVALSUM = LOGXVALSUM + FNdN*LOGXVALI;
        CEXTSUM    = CEXTSUM    + FNdN*CEXT(I);
        CSCASUM    = CSCASUM    + FNdN*CSCA(I);
        CABSSUM    = CABSSUM    + FNdN*CABS(I);


    %   Accumulate third moments...
        QEXTSUM = QEXTSUM + FNXSdN *QEXT(I);  
        QSCASUM = QSCASUM + FNXSdN *QSCA(I);
        QABSSUM = QABSSUM + FNXSdN *QABS(I);
        XCSUM   = XCSUM   + FNXSdN *XVALI;


    %   Accumulate Qsca - weighted third moments...
        ASYMSUM = ASYMSUM + FNXSdN *QSCA(I) *ASYM(I);
        P11SUM = P11SUM + FNXSdN *QSCA(I) *P11(:,I);
    end
end

% Clean up & do a check...
if (ICOUNT > 0)
    LOGXVALSUM = LOGXVALSUM/SUMN;
    XVALSUM = exp(LOGXVALSUM);
    ERROR = (XVALSUM/XBAR - 1);
    CEXTSUM    = CEXTSUM/SUMN;
    CSCASUM    = CSCASUM/SUMN;
    CABSSUM    = CABSSUM/SUMN;
    QEXTSUM = QEXTSUM/SUMNS;
    QSCASUM = QSCASUM/SUMNS;
    QABSSUM = QABSSUM/SUMNS;
    X2SUM   = SUMNS/SUMN;
    X3SUM   = XCSUM/SUMN;
    XCSUM   = XCSUM/SUMNS;
    ASYMSUM = (ASYMSUM/SUMNS)/QSCASUM;
    P11SUM = (P11SUM/SUMNS)/QSCASUM;
else
    XVALSUM = 9999;
end


% Get out ...
return

