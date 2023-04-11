function H = multivar_trunc(Mu, Sigma)
    %% Truncated multivariate Gaussian distribution entropy (numerical)
    % lambda, k
    % Mu = [1.4598, 0.5416]; % For amp = 1 of wavy sweep
    % Mu = [1.4499, 0.5405]; % For amp = 0.6 of wavy sweep
    ml = Mu(1);
    mk = Mu(2);
    % Sigma = [0.0383, -0.0019; -0.0019, 0.0182]; % For amp = 1 of wavy sweep
    % Sigma = [0.0156, -0.0016; -0.0016, 0.0122]; % For amp = 0.6 of wavy sweep
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
    pdf = @(l,k)exp(-(1/(2.*(1-rho^2))).*(((l-ml)./sl).^2)-2.*rho.*((l-ml)./sl).*((k-mk)./sk)+((k-mk)/sk).^2)./(2*pi*sl*sk*sqrt(1-rho.^2))
    
    % denominator of truncated Gaussian is cdf evaluated at bounds and
    % subtracted --> bounded integral of pdf 
    denom = integral2(pdf,lmin,lmax,kmin,kmax)
    
    % Truncated distribution is just original pdf divided by this integral (I
    % think)
    p_trunc = @(l,k)(1/denom).*exp(-(1/(2.*(1-rho^2))).*(((l-ml)./sl).^2)-2.*rho.*((l-ml)./sl).*((k-mk)./sk)+((k-mk)/sk).^2)./(2*pi*sl*sk*sqrt(1-rho.^2));
    entr = @(l,k)(1/denom).*exp(-(1/(2.*(1-rho^2))).*(((l-ml)./sl).^2)-2.*rho.*((l-ml)./sl).*((k-mk)./sk)+((k-mk)/sk).^2)./(2*pi*sl*sk*sqrt(1-rho.^2)).*log((1/denom).*exp(-(1/(2.*(1-rho^2))).*(((l-ml)./sl).^2)-2.*rho.*((l-ml)./sl).*((k-mk)./sk)+((k-mk)/sk).^2)./(2*pi*sl*sk*sqrt(1-rho.^2)));
    % Entropy is integral of -plog(p)
    % syms L K
    H = integral2(entr,lmin,lmax,kmin,kmax)
end
