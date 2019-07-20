function p=mylegendreP(l,m,x)

%disp ('here in legendreP')

% Calculate the Legendre polynomials

% Calculate coef of maximum degree in x from the explicit analytical
% formula
cl=(-1)^m*factorial(2*l)/((2^l)*factorial(l)*factorial(l-m));
px=l-m;

% Power of x changes from one term to the next by 2. Also needed for
% sqrt(1-x^2).
x2=x.*x;


% Calculate efficiently P_l^m (x)/sqrt(1-x^2)^(m/2) - that is, only the
% polynomial part. At least one coefficient is guaranteed to exist - there
% is no null Legendre polynomial.
p=cl*ones(size(x));

for j=l-1:-1:0,
    % Check the exponent of x for current coefficient, px. If it is 0 or 1,
    % just exit the loop
    if(px<2), break; end;
    % If current exponent is >=2, there is a "next" coefficient; multiply p
    % by x2 and add it. Calculate the current coefficient
    cl=-(j+j+2-l-m)*(j+j+1-l-m)/(2*(j+j+1)*(l-j))*cl;
    % ...and add to the polynomial
    p=p.*x2 + cl;
    % Decrease the exponent of x - this is the exponent of x corresponding
    % to the newly added coefficient
    px=px-2;
end;
% Now we're done adding coefficients. However, if the exponent of x
% corresponding to the last added coefficient is 1 (polynomial is odd),
% multiply the polynomial by x 
if(px==1), p=p.*x; end;

% All that's left is to multiply the whole thing with sqrt(1-x^2)^(m/2). No
% further calculations are needed if m=0.
if(m==0), return; end;

x2=1-x2;
%First, multiply by the integer part of m/2
for j=1:floor(m/2), p=p.*x2; end;
%If m is odd, there is an additional factor sqrt(1-x^2)
if(m~=2*floor(m/2)), p=p.*sqrt(x2); end;
% Done.


