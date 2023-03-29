%% Truncated multivariate Gaussian distribution entropy (numerical)
% lambda, k
Mu = [1.4598, 0.5416];
ml = Mu(1);
mk = Mu(2);
Sigma = [0.0383, -0.0019; -0.0019, 0.0182];
% std dev of each variable
sl = sqrt(Sigma(1,1));
sk = sqrt(Sigma(2,2));
rho = Sigma(1,2)/(sk*sl);

% bounds
lmax = 2;
lmin = 1;
kmax = 1;
kmin = 0;

% l = lambda = 'x', k = 'y'
% syms l k
pdf = @(l,k)exp(-(1/(2*(1-rho^2)))*(((l-ml)/sl)^2)-2*rho*((l-ml)/sl)*((k-mk)/sk)+((k-mk)/sk)^2)./(2*pi*sl*sk*sqrt(1-rho^2))

% denominator of truncated Gaussian is cdf evaluated at bounds and
% subtracted --> bounded integral of pdf 
denom = integral2(pdf,lmin,lmax,kmin,kmax)

% Truncated distribution is just original pdf divided by this integral (I
% think)
p_trunc = pdf/denom

% Entropy is integral of -plog(p)
H = integral2(p_trunc*log(p_trunc),lmin,lmax,kmin,kmax)

