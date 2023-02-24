%% Fit Data to Gaussian distribution
clear;
clc;
close all;
%% Generate data (read in from wavy sweep sim)
% X = randn(1000,2);
load('22-1215-Wavy_Sweep\sensitivity.mat');
% loop through all iterations of sweep (11 total iterations)
amplitudes = linspace(0,2,11);
all_diff_ent = zeros(1,11);
all_pwise_ent = all_diff_ent;
all_outlier_ent = all_diff_ent;
all_outlier_contrib = all_diff_ent;
for ii=1:11
    % X1 = lam, X2 = k
    X = [lam{ii}' k{ii}'];
    %% visualize raw data - scatter plot
%     figure()
%     scatter(X(:,1),X(:,2),'.')
%     xlabel('\lambda')
%     ylabel('k')
%     title(['Amplitude = ',num2str(amplitudes(ii))])
%     xlim([1,2])
    
    %% generate histogram
    figure()
    [lam_bin,k_bin,h] = ndhist(X,'axis',[1,0;2,1]);
    colormap turbo
    hold on

    % Calculate pointwise (binwise) entropy
    p = h'./sum(h(:));
    H_pwise = -sum(p.*log(p), 'all','omitnan')
    all_pwise_ent(ii) = H_pwise;
    %% use fitgmdist to fit to gaussian distribution
    GMModel = fitgmdist(X,1);
    
    % visualize data with Gaussian fit
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    g = gca;
    fcontour(gmPDF, [g.XLim,g.YLim], 'LineWidth',3, 'LineColor','flat')
    title(['Amplitude = ',num2str(amplitudes(ii))])
    hold off
    
    % calculate entropy of Gaussian
    D = 2; % 2 independent variables --> dim(X) = 2
    Sigma = GMModel.Sigma
    Mu = GMModel.mu
    diff_entropy = (1/2)*(log(det(Sigma)) + (D/2)*(1 + log(2*pi)))
    all_diff_ent(ii) = diff_entropy;
%     Sigma1 = GMModel.Sigma(:,:,1);
%     Sigma2 = GMModel.Sigma(:,:,2);
%     diff_entropy = (1/2)*(log(det(Sigma1)) + (D/2)*(1 + log(2*pi))) + (1/2)*(log(det(Sigma2)) + (D/2)*(1 + log(2*pi)))
    
    %% Visualize surface (Gaussian and hist)
    figure()
    k_ = linspace(0,1,length(k_bin));
    lam_ = linspace(1,2,length(lam_bin));
    [K,LAM] = meshgrid(k_,lam_);
    surf(LAM,K,gmPDF(LAM,K))
    hold on
    % scale up p to make it visible
    p_scaled = p.*(max(max(gmPDF(LAM,K))/max(max(p))));
    surf(LAM,K,p_scaled)
    legend('Gaussian fit','histogram probability')
    colormap turbo
    title('Gaussian fit surface vs scaled histogram probability')
    xlabel('\lambda')
    ylabel('k')
    zlabel('probability')
    hold off
    % error between Gaussian fit and scaled hist probability
    figure()
    surf(LAM,K,p_scaled-gmPDF(LAM,K))
    xlabel('\lambda')
    ylabel('k')
    zlabel('error')
    title('scaled hist - Gaussian fit')

    %% Look at histogram points "outside" of the Gaussian and contribution
    % to entropy
    % Look at locations where Gaussian is <10% of max value (chosen arbitrarily)
    cutoff = max(max(gmPDF(LAM,K)))*.1;
    H_outlier = -sum(p(gmPDF(LAM,K)<cutoff).*log(p(gmPDF(LAM,K)<cutoff)), 'all','omitnan')
    all_outlier_ent(ii) = H_outlier;
    outlier_contrib = H_outlier/H_pwise
    all_outlier_contrib(ii) = outlier_contrib;
    % Visualise this outlier data
    % "turn off" the region inside the Gaussian
    p_scaled(gmPDF(LAM,K)>=cutoff)=0;
    figure()
    surf(LAM,K,p_scaled)
    colormap turbo
    title('points outside Gaussian')
    xlabel('\lambda')
    ylabel('k')
    zlabel('probability')
    zlim([0 max(max(gmPDF(LAM,K)))])

end
%% Plot all entropy values
figure()
plot(amplitudes,all_diff_ent,'LineWidth',2)
hold on
plot(amplitudes,all_pwise_ent,'LineWidth',2)
plot(amplitudes,all_outlier_ent,'LineWidth',2)
hold off
xlabel('Wave amplitude')
ylabel('simulation entropy')
legend('Gaussian entropy','pointwise entropy','outlier entropy','Location','northwest')
% xlim([0,0.8])
figure()
plot(amplitudes,all_outlier_contrib*100,'LineWidth',2)
title('relative entropy contribution of outlier data')
xlabel('Wave amplitude')
ylabel('percentage contribution')