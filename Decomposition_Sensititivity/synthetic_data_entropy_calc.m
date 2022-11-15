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
sigma_k = 0.05;
sigma_lam = 0.3;
% max number of data points
n = 8;

entropy_all = ones(1,n-2);
p2ent_all = ones(1,n-2);

for i=3:n
    % Number of data points
    N = 10^i;
    % Gaussian distribution itself
    % k = (1/sigma_k*sqrt(2*pi)).*exp(-(1/2).*((k_range-mu_k)./sigma_k).^2);
    % lambda = (1/sigma_lam*sqrt(2*pi)).*exp(-(1/2).*((lam_range-mu_lam)./sigma_lam).^2);
    k = normrnd(mu_k,sigma_k,[1,N]);
    lambda = normrnd(mu_lam,sigma_lam,[1,N]);
    % rescale to fit bounds
    k = k_min + (k_max-k_min).*((k-min(k))./(max(k)-min(k)));
    lambda = lam_min + (lam_max-lam_min).*((lambda-min(lambda))./(max(lambda)-min(lambda)));
    
    %% Plot 1D hist as sanity check
    figure()
    subplot(1,2,1)
    histogram(k,10)
    title('k distribution')
    subplot(1,2,2)
    histogram(lambda,10)
    title('\lambda distribution')
    
    %% 2D Hist (figure command only)
    figure()
    sgtitle(strcat('N = ',num2str(N)))
    %% 2D Hist
    subplot(2,2,1)
    [lam_bin,k_bin,h] = ndhist(lambda,k,'axis',[1,0;1.8,1]-0.005);
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