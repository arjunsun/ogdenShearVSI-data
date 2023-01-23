%% Fit Data to Gaussian distribution
clear;
clc;
close all;
%% Generate data (read in from wavy sweep sim)
% X = randn(1000,2);
load('22-1215-Wavy_Sweep\sensitivity.mat');
% loop through all iterations of sweep (11 total iterations)
amplitudes = linspace(0,2,11);
all_ent = zeros(1,11);
for ii=1:11
    % X1 = lam, X2 = k
    X = [lam{ii}' k{ii}'];
    %% visualize raw data
    figure()
    scatter(X(:,1),X(:,2),'.')
    xlabel('\lambda')
    ylabel('k')
    title(['Amplitude = ',num2str(amplitudes(ii))])
    xlim([1,2.5])
    hold on
    
    %% use fitgmdist to fit to gaussian distribution
    GMModel = fitgmdist(X,1);
    
    % visualize data with Gaussian fit
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    g = gca;
    fcontour(gmPDF, [g.XLim,g.YLim], 'LineWidth',3)
    hold off
    
    % calculate entropy of Gaussian
    D = 2; % 2 independent variables --> dim(X) = 2
    Sigma = GMModel.Sigma
    Mu = GMModel.mu
    entropy = (1/2)*(log(det(Sigma)) + (D/2)*(1 + log(2*pi)))
    all_ent(ii) = entropy;
end
% Plot all entropy values
figure()
plot(amplitudes,all_ent)
xlabel('Wave amplitude')
ylabel('simulation entropy')
xlim([0,0.8])