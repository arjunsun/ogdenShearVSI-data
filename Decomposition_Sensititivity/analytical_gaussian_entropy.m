%% Analytical expected entropy
mu_k = 0.5;
mu_lam = 1.4;

syms l_ k_;
% analytical expression for entropy of a Gaussian distribution
% entropy = -plog(p)

sigma_range = (0.1:0.1:5);
ea = zeros(length(sigma_range));
ent_closedform = ea;
for ii = 1:length(sigma_range)
    for jj = 1:length(sigma_range)
        sigma_k = sigma_range(ii);
        sigma_lam = sigma_range(jj);
        entr_ = @(l_,k_)(-(1/(2*pi*sigma_k*sigma_lam))*exp(-((mu_lam-l_).^2)./(2*sigma_lam.^2) - ((mu_k-k_).^2)./(2.*sigma_k.^2))).*log((1./(2*pi.*sigma_k.*sigma_lam)).*exp(-((mu_lam-l_).^2)./(2.*sigma_lam.^2) - ((mu_k-k_).^2)./(2.*sigma_k.^2)));
        k_min = 0;% mu_k-20*sigma_k; % 0
        k_max = 1;%mu_k+20*sigma_k; % 1
        lam_min = 1;%mu_lam-20*sigma_lam; % 1
        lam_max = 1.8;%mu_lam+20*sigma_lam; % 1.8
        entropy_analytical = integral2(entr_,lam_min,lam_max,k_min,k_max);
        ea(ii,jj) = entropy_analytical;
        ent_closedform(ii,jj) = log(sigma_k*sigma_lam) + (1+log(2*pi));
    end
end    
figure()
surf(sigma_range,sigma_range,ea)
hold on
surf(sigma_range,sigma_range,ent_closedform)
xlabel('\sigma_k')
ylabel('\sigma_{\lambda}')
zlabel('entropy')
legend('integrated (bounded)', 'closed form equation (full distribution)')

%% 1 dimensional Gaussian
mu = 0;%.5;

sigma_range = (1:1:100);
ea = zeros(size(sigma_range));
ent_closedform = ea;
syms x
for ii = 1:length(sigma_range)
    sigma = sigma_range(ii);
    x_min = -1000;%mu-20*sigma;
    x_max = 1000;%mu+20*sigma;
    p_ = @(x)(1/(2*pi*sigma))*exp(-(1/2).*((x-mu)./sigma).^2);
    P_ = @(x)(1/2)*(1+erf(x/sqrt(2)));
    p2_ = @(x)p_(x)/(sigma*(P_((x_max-mu)/sigma) - P_((x_min-mu)/sigma)));
    entr_ = @(x)-p2_(x).*log(p2_(x));
    
    entropy = integral(entr_,x_min,x_max);
    ea(ii) = entropy;
    ent_closedform(ii) = 0.5 + 0.5*log(2*pi*sigma*sigma);
end
figure()
plot(sigma_range,ea,sigma_range,ent_closedform)
legend('integrated', 'closed form equation')
xlabel('\sigma')
ylabel('entropy')

%% comparing peak and flat distribution
% flat distribution on -b to b --> entropy = log(2b)
syms b_
% exponential w/ peak at 0 going done on both sides
A_ = 1/(2*(1-exp(-b_)));
ent_exp = @(b_)-2.*A_.*(-exp(-b_).*(log(A_.*exp(-b_))-1) + log(A_) - 1)
b = logspace(-2,2);
A = 1./(2.*(1-exp(-b)));
entexp_vals = -2.*A.*(-exp(-b).*(log(A.*exp(-b))-1) + log(A) - 1);
semilogx(b,log(2*b),b,entexp_vals)
legend('flat distribution','exponential peak')