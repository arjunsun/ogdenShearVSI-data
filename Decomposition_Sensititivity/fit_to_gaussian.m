%% Fit Data to Gaussian distribution
clear;
clc;
close all;
%% Flags to turn stuff on and off
% 1 if you want x, 0 if not
% histogram and gaussian fit always on 
plot_gauss = 1;
surfs_on = 0;%1;
outside_gauss = 0;%1;
all_ent_plot = 0;%1;
check_bins = 0;
gauss_data = 0;
sig_effect = 1;%0;
%% Generate data (read in from wavy sweep sim)
% X = randn(1000,2);
load('22-1215-Wavy_Sweep\sensitivity.mat');
% loop through all iterations of sweep (11 total iterations)
amplitudes = linspace(0,2,11);
all_diff_ent = zeros(1,11);
all_pwise_ent = all_diff_ent;
all_outlier_ent = all_diff_ent;
all_outlier_contrib = all_diff_ent;

for ii=6%1:11
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
    binfactor = 10;
    [lam_bin,k_bin,h] = ndhist(X,'axis',[1,0;2,1],'bins',binfactor);
    colormap turbo
    hold on
    xlabel('\lambda')
    ylabel('k')

    % Calculate pointwise (binwise) entropy
    p = h'./sum(h(:));
    H_pwise = -sum(p.*log(p), 'all','omitnan')
    all_pwise_ent(ii) = H_pwise;
    %% use fitgmdist to fit to gaussian distribution
    GMModel = fitgmdist(X,1);
    
    % visualize data with Gaussian fit
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    if plot_gauss
        g = gca;
        fcontour(gmPDF, [g.XLim,g.YLim], 'LineWidth',3, 'LineColor','flat')
        title(['Amplitude = ',num2str(amplitudes(ii))])
    end
    hold off
    
    % calculate entropy of Gaussian
    D = 2; % 2 independent variables --> dim(X) = 2
    Sigma = GMModel.Sigma;
    Mu = GMModel.mu;
    diff_entropy = (1/2)*(log(det(Sigma)) + (D/2)*(1 + log(2*pi)))
    all_diff_ent(ii) = diff_entropy;
%     Sigma1 = GMModel.Sigma(:,:,1);
%     Sigma2 = GMModel.Sigma(:,:,2);
%     diff_entropy = (1/2)*(log(det(Sigma1)) + (D/2)*(1 + log(2*pi))) + (1/2)*(log(det(Sigma2)) + (D/2)*(1 + log(2*pi)))
    
    %% Visualize surface (Gaussian and hist)
    if surfs_on
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
    end
    %% Look at histogram points "outside" of the Gaussian and contribution
    % to entropy
    if outside_gauss
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
end
%% Plot all entropy values
if all_ent_plot
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
end

%% Look at histogram binning effects for data from Gaussian
if gauss_data
    X = mvnrnd(Mu,Sigma,10000);
    figure()
    scatter(X(:,1),X(:,2),'.')
    xlabel('\lambda')
    ylabel('k')
    xlim([1,2])
end
%% Effect of bin size/number
if check_bins
    bf_range = logspace(-0.5,1.5,20);
    nbins = zeros(size(bf_range));
    H_bins = nbins;
    for ii=1:length(bf_range)
        figure()
        binfactor = bf_range(ii);
        [lam_bin,k_bin,h] = ndhist(X,'axis',[1,0;2,1],'bins',binfactor);
        colormap turbo
        title(['binfactor = ',num2str(binfactor)])
        
        % calculate number of bins
        nbins(ii) = numel(h);
    
        % Calculate pointwise (binwise) entropy
        p = h'./sum(h(:));
        H_pwise = -sum(p.*log(p), 'all','omitnan')
        H_bins(ii) = H_pwise;
    end
    
    % Baseline values (for bin size = 0.01, hard coded, change with dataset)
    bins_baseline = 10201;
    H_baseline = 7.9139;
    
    % Theoretical max - if each point is unique and in it's own bin
    N = numel(X)/2;
    H_max_theoretical = log(N);
    
    % Plot entropy vs bins
    figure()
    semilogx(nbins,H_bins,'LineWidth',2);
    hold on
    
    plot(nbins, H_max_theoretical.*ones(size(nbins)),'--k')
    xlabel('number of bins')
    ylabel('entropy')
    xlim([min(nbins) max(nbins)])
    if gauss_data
        plot(nbins,diff_entropy.*ones(size(nbins)),'--m')
        plot(nbins, log(nbins),'--r')
        legend('vary with bin size','max possible','value from Gaussian','log(nbins)','Location','east')
        hold off
    else
        plot(bins_baseline,H_baseline,'r*','MarkerSize',10)
        legend('vary with bin size','max possible','baseline (0.01 bin width)','Location','southeast')
        hold off
    end
end
%% Effect of sigma
if sig_effect
    % Keep mu constant, scale down sigma by just multiplying by a factor <1
    scale = linspace(0.1,1,5);
    Sigma_base = Sigma;

    % Look at sigma's effect on both binwise and the analytical form of
    % entropy
    sig_diff_ent = zeros(size(scale));

    % Will be plotting all bin vs entropy plots on the same figure to see
    % impact of sigma, arbitrarily choose to plot on figure 99 in case
    % other parts of code have generated a lot of plots
    figure(99)
    legend_labels = cell(1,length(scale)+2);
    for jj = 1:length(scale)
        % Scale down sigma and then generate data
        Sigma = scale(jj).*Sigma_base;
        % X = mvnrnd(Mu,Sigma,10000);
        X = mvnrnd_trn([1,0],[2,1],Mu,Sigma,10000);
        
        % Look at various bin sizes (same process as above section)
        % Again vary bin size/number by varying bin factor in ndhist
        bf_range = logspace(-0.5,1.5,10);
        nbins = zeros(size(bf_range));
        H_bins = nbins;
        for ii=1:length(bf_range)
            % Don't actually care about each histogram, so just have them
            % plot over each other in figure 100
            figure(100)
            binfactor = bf_range(ii);
            [lam_bin,k_bin,h] = ndhist(X,'axis',[1,0;2,1],'bins',binfactor);
            colormap turbo
            title(['binfactor = ',num2str(binfactor)])
            
            % calculate number of bins
            nbins(ii) = numel(h);
        
            % Calculate pointwise (binwise) entropy
            p = h'./sum(h(:));
            H_pwise = -sum(p.*log(p), 'all','omitnan')
            H_bins(ii) = H_pwise;
        end
        % Plot entropy vs bins
        figure(99)
        semilogx(nbins,H_bins,'LineWidth',2);
        hold on
        % Get legend entry (add all entries at the end)
        legend_labels{jj} = strcat(num2str(scale(jj)),'*sigma');

        % Also calculate entropy of Gaussian dist (analytical form)
        D = 2;
        sig_diff_ent(jj) = (1/2)*(log(det(Sigma)) + (D/2)*(1 + log(2*pi)));
    end
    
    figure(99)
    % also plot theoretical max entropy (log(N)) for both overall number of
    % data points and per num of bins
    b = logspace(2,7,20);
    N = numel(X)/2;
    plot(b,log(N).*ones(size(b)),'--')
    plot(b,log(b),'--')
    
    % also add these to the legend
    legend_labels{jj+1} = 'overall max possible';
    legend_labels{jj+2} = 'max for num bins';
    
    % figure labeling etc.
    xlim([10^2 10^7]);
    ylim([4 10]);
    xlabel('number of bins')
    ylabel('binwise entropy')
    title('binwise entropy, effect of sigma')
    legend(legend_labels,'Location','southeast')
    hold off

    % Plot result of analytical form
    figure()
    plot(scale,sig_diff_ent)
    xlabel('scale factor on sigma')
    ylabel('entropy of gaussian')
    title('differential entropy, effect of sigma')
end