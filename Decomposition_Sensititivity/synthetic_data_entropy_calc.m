%% Generate Gaussian distribution for k and lambda
% Parameters for Gaussian
% N = 10000; % num data points
k_min = 0;
k_max = 1;
% k_range = linspace(k_min,k_max,N);
lam_min = 1;
lam_max = 1.8;
% lam_range = linspace(lam_min,lam_max,N);
% averages and std deviations
mu_k = 0.5;
mu_lam = 1.4;
sigma_k = 0.5;
sigma_lam = 0.5;
% max number of data points
n = 8;

%% Analytical expected entropy
% mu_k = 0.5;
% mu_lam = 1.4;
% 
% syms l_ k_;
% % analytical expression for entropy of a Gaussian distribution
% % entropy = -plog(p)
% 
% sigma_range = (0.1:0.1:5);
% ea = zeros(length(sigma_range));
% for ii = 1:length(sigma_range)
%     for jj = 1:length(sigma_range)
%         sigma_k = sigma_range(ii);
%         sigma_lam = sigma_range(jj);
%         entr_ = @(l_,k_)(-(1/(2*pi*sigma_k*sigma_lam))*exp(-((mu_lam-l_) ...
%             .^2)./(2*sigma_lam.^2) - ((mu_k-k_).^2)./(2.*sigma_k.^2))) ...
%             .*log((1./(2*pi.*sigma_k.*sigma_lam)).*exp(-((mu_lam-l_).^2) ...
%             ./(2.*sigma_lam.^2) - ((mu_k-k_).^2)./(2.*sigma_k.^2)));
%         entropy_analytical = integral2(entr_,mu_lam-20*sigma_lam, ...
%             mu_lam+20*sigma_lam,mu_k-20*sigma_k,mu_k+20*sigma_k);
%         ea(ii,jj) = entropy_analytical;
%     end
% end    
% surf(sigma_range,sigma_range,ea)
% xlabel('\sigma_k')
% ylabel('\sigma_{\lambda}')
% zlabel('entropy')

%% make Gaussians
entropy_all = ones(1,n-2);
p2ent_all = ones(1,n-2);

for i=3:n
    % Number of data points
    N = 10^i;
    % Gaussian distribution itself
    % k = (1/sigma_k*sqrt(2*pi)).*exp(-(1/2).*((k_range-mu_k)./sigma_k).^2);
    % lambda = (1/sigma_lam*sqrt(2*pi)).*exp(-(1/2).*((lam_range-mu_lam) ...
    % ./sigma_lam).^2);
    k = normrnd(mu_k,sigma_k,[1,N]);
    lambda = normrnd(mu_lam,sigma_lam,[1,N]);
    % rescale to fit bounds
    k = k_min + (k_max-k_min).*((k-min(k))./(max(k)-min(k)));
    lambda = lam_min + (lam_max-lam_min).*((lambda-min(lambda))./(max ...
        (lambda)-min(lambda)));
    %     lambda = linspace(lam_min,lam_max,N);
    
    %% Plot 1D hist as sanity check
%     figure()
%     subplot(1,2,1)
%     histogram(k,10)
%     title('k distribution')
%     subplot(1,2,2)
%     histogram(lambda,10)
%     title('\lambda distribution')
    
    %% 2D Hist (figure command only)
    figure()
    sgtitle(strcat('N = ',num2str(N)))
    %% 2D Hist
    subplot(2,2,1)
    [~,~,h] = ndhist(lambda,k,'axis',[1,0;1.8,1]-0.005);
    imagesc(h);
    colormap turbo
    set(gca,'YDir','normal')
    axis square
    ylabel('k')
    xlabel('$\lambda$','interpreter','latex')
    title('Count (h)')
    colormap turbo
    
    subplot(2,2,2)
    imagesc(log10(h./sum(h,'all')))
    colormap turbo
    clim([-5 0])
    set(gca,'YDir','normal')
    axis square
    ylabel('k')
    xlabel('$\lambda$','interpreter','latex')
    title('Probability (p), log scale')
    colormap turbo
    
    subplot(2,2,3)
    [~,~,h] = ndhist(lambda,k,'axis',[1,0;1.8,1]-0.005);
    p = h./sum(h,'all');
    imagesc(log10(-p.*log10(p)));
    clim([-5 0]);
    set(gca,'YDir','normal')
    axis square
    ylabel('k')
    xlabel('$\lambda$','interpreter','latex')
    title('$plog(p)$, log scale','interpreter','latex')
    colormap turbo
    
    subplot(2,2,4)
    [lam_bin,k_bin,h] = ndhist(lambda,k,'axis',[1,0;1.8,1]-0.005);
    p = h./sum(h,'all');
    imagesc(log10(-p.*p.*log10(p)));
    clim([-8 0]);
    set(gca,'YDir','normal')
    axis square
    ylabel('k')
    xlabel('$\lambda$','interpreter','latex')
    title('$p^2log(p)$, log scale','interpreter','latex')
    colormap turbo
    
    %% Calculate total entropy
    entropy = sum(-p.*log10(p),"all","omitnan");
    p2ent = sum(-p.*p.*log10(p),'all','omitnan');
    entropy_all(i-2) = entropy;
    p2ent_all(i-2) = p2ent;
end
%% Plot N vs entropy
figure()
subplot(1,2,1)
semilogx(logspace(3,8,6),entropy_all)
xlabel('Number of data points')
ylabel('entropy (plogp)')

subplot(1,2,2)
semilogx(logspace(3,8,6),p2ent_all)
xlabel('Number of data points')
ylabel('entropy (p^2logp)')

%% Generate violin plots
% Generate 10(+) entropy values for each value of N, columns are different
% N
% clc;
num_samples = 25;
n_max = 4;
entropy_all = ones(num_samples,n_max-1);
ent_anl_all = entropy_all;
figure()
for i = 2:n_max
    N = 10^i;
    for j = 1:num_samples
        % using Gaussian parameters from above generate Gaussian
        k = normrnd(mu_k,sigma_k,[1,N]);
        lambda = normrnd(mu_lam,sigma_lam,[1,N]);
        % rescale to fit bounds
%         k = k_min + (k_max-k_min).*((k-min(k))./(max(k)-min(k)));
%         lambda = lam_min + (lam_max-lam_min).*((lambda-min(lambda)) ...
%             ./(max(lambda)-min(lambda)));
        
        [~,~,h] = ndhist(lambda,k,'axis',[1,0;1.8,1]-0.005);
        p = h./sum(h,'all');
        entropy = sum(-p.*log10(p),"all","omitnan");
        entropy_all(j,i-1) = entropy;
        [LAM,K] = meshgrid(lambda,k); 
        ent_anl_disc = sum(entr_(LAM,K),'all');
        ent_anl_all(j,i-1) = ent_anl_disc;
    end
end
violin(entropy_all);%, 'xlabel',{'2','3','4','5','6','7'});
xlabel('log10(N)-1')
ylabel('entropy');
ylim([1 4]);
