%% Run sensitivity for S_alpha functions
sensitivity_metric_scrap
close all

type = {'Sim','Exp'};
fn = 'cax-5to0LOG';
clim = [-5 0];%[0 0.05];

ii = 2;

%% k vs. lam Simulations
load('22-0325-Ogden_Uniaxial/Simulations/sensitivity.mat')
figure();
[lam_bin,k_bin,h] = ndhist(lam{1}(:),k{1}(:),'axis',[1,0;1.8,1]-0.005);
LAM_BIN_sim = lam_bin; K_BIN_sim = k_bin; H_sim = h';
axis square
title(['Uniaxial for $U_d=$7 mm'],'interpreter','latex')
ylabel('k')
xlabel('$\lambda$','interpreter','latex')
imagesc(log10(h/sum(h(:)))); caxis(clim); 
colormap turbo;
set(gca,'YDir','normal'); 
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/Sim/histo_uniaxial.png']);
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/Sim/histo_uniaxial.pdf']);


%% k vs. lam Experiment
load('22-0325-Ogden_Uniaxial/Experimental/sensitivity.mat')
load('22-0325-Ogden_Uniaxial/Experimental/MRI-3Ddefs_RectPrism_190919_222hfilt.mat')

for i = 1:10
    figure();
    [lam_bin,k_bin,h] = ndhist(lam{i}(:),k{i}(:),'axis',[1,0;1.8,1]-0.005);
    LAM_BIN{i} = lam_bin; K_BIN{i} = k_bin; H{i} = h';
    axis square
    title(['Uniaxial for $U_d=$' num2str(prescribedU(i)) ' mm'],'interpreter','latex')
    ylabel('k')
    xlabel('$\lambda$','interpreter','latex')
    imagesc(log10(h/sum(h(:)))); caxis(clim); 
    colormap turbo;
    set(gca,'YDir','normal'); 
    saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/histo_uniaxial' num2str(prescribedU(i)) '.png']);
    saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/histo_uniaxial' num2str(prescribedU(i)) '.pdf']);
end

%% k vs. lam Sensitivity Plots

figure; % imagesc(squeeze(Salpha_all(:,1,:))') % a = -5
[LAM,K] = ndgrid(lam_bin,k_bin);
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,1,:)),500);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=-5$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_alpha-5.png']);
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_alpha-5.pdf']);

figure; % imagesc(squeeze(Salpha_all(:,1001,:))') % a = 5
[~,hand] = contourf(LAM,K,squeeze(Salpha_all(:,1001,:)),500);
set(hand,'LineColor','none');
axis square;
title([type{ii} ': $\alpha=5$'],'interpreter','latex')
ylabel('k'); xlabel('$\lambda$','interpreter','latex')
axis square;
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_alpha5.png']);
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_alpha5.pdf']);

%% S_alpha{test}(metric;alpha) Plots Experimental
% Tests: 
% Uniaxial - [0, 0.5, 1, 2, 3, 4, 4.5, 5, 6, 7];
figure();
for test = 1:length(prescribedU)
    for i = 1:length(alpha)
        mat_temp{i} = squeeze(Salpha_all(:,i,:));
        % mat_temp{i} = mat_temp{i}./sum(mat_temp{i}(:)); % Normalized
        S_lookup{test}(i,1) = alpha(i);
        S_lookup{test}(i,2) = sum(H{test}(:).*mat_temp{i}(:));
    end
    if test == length(prescribedU)
        S_hat_lookup = S_lookup;
        for j = 1:length(prescribedU)
            nrmlz(j) = sum(H{j}(:));
            S_hat_lookup{j}(:,2) = S_lookup{j}(:,2)/sum(H{j}(:)); % Normalized
            plot(S_hat_lookup{j}(:,1),S_hat_lookup{j}(:,2),'Color',[j/length(prescribedU) 0 0])
            hold on
        end
    end
end
H_ = H; S_hat_ = S_hat_lookup; S_ = S_lookup; Salpha_all_ = mat_temp;
title(['$S_{metric}$(test,$\alpha$) for ' type{ii}],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_metric.png'])
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_metric.pdf'])

save(['k_lam_SensitivityPlots_uniaxial/S_metrics_histos_' type{ii} '.mat'],'H','S_hat_lookup','S_lookup');

%% S_alpha{test}(metric;alpha) Plots Simulation
figure();
for i = 1:length(alpha)
    mat_temp{i} = squeeze(Salpha_all(:,i,:));
    % mat_temp{i} = mat_temp{i}./sum(mat_temp{i}(:)); % Normalized
    S_lookup_sim(i,1) = alpha(i);
    S_lookup_sim(i,2) = sum(H_sim(:).*mat_temp{i}(:));
end
S_hat_lookup_sim = S_lookup_sim;
nrmlz(11) = sum(H_sim(:));
S_hat_lookup_sim(:,2) = S_lookup_sim(:,2)/sum(H_sim(:)); % Normalized
plot(S_hat_lookup_sim(:,1),S_hat_lookup_sim(:,2),'Color',[j/length(prescribedU) 0 0])
title(['$S_{metric}$(test,$\alpha$) for Sim'],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_metric.png'])
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/' type{ii} '/sens_metric.pdf'])

Salpha_ALL = mat_temp;

H_{11} = H_sim; S_hat_{11} = S_hat_lookup_sim; S_{11} = S_lookup_sim;

%% S_alpha{test}(metric;alpha) Plots (Both sim and exp)
save(['k_lam_SensitivityPlots_uniaxial/S_metrics_histos_both.mat'],'H_','S_hat_','S_','Salpha_all_');
figure()
plot(S_hat_{11}(:,1),S_hat_{11}(:,2))
hold on
plot(S_hat_{4}(:,1),S_hat_{4}(:,2))
plot(S_hat_{8}(:,1),S_hat_{8}(:,2))
plot(S_hat_{10}(:,1),S_hat_{10}(:,2))
title(['$S_{metric}$(test,$\alpha$)'],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend('Sim; Uniaxial; $U_d = 7 mm$','Exp; Uniaxial; $U_d = 7 mm$','Exp; Uniaxial; $U_d = 5 mm$',...
    'Exp; Uniaxial; $U_d = 2 mm$','interpreter','latex')

saveas(gcf,['k_lam_SensitivityPlots_uniaxial/sens_metric_both.png'])
saveas(gcf,['k_lam_SensitivityPlots_uniaxial/sens_metric_both.pdf'])

%% Plot both uniaxial and hists

load('k_lam_SensitivityPlots/S_metrics_histos_both.mat')
for ii = 1:length(H_)
    if ii == 1
        for jj = 1:length(H_{ii})
            plot(S_hat_{ii}{jj}(:,1),S_hat_{ii}{jj}(:,2))
        end
    else
        for jj = 1:length(H_{ii})
            plot(S_hat_{ii}{jj}(:,1),S_hat_{ii}{jj}(:,2))
        end
    end
end
title(['$S_{metric}$(test,$\alpha$)'],'interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
ylabel('$S_{metric}=\sum_{i,j}\hat{hist}\left(k_i,\lambda_j,test\right)\;.*\;\hat{S}_\alpha\left(k_i,\lambda_j,\alpha\right)$',...
    'interpreter','latex')
legend('Sim; Uniaxial; $U_d = 7 mm$','Exp; Uniaxial; $U_d = 7 mm$','Exp; Uniaxial; $U_d = 5 mm$',...
    'Exp; Uniaxial; $U_d = 2 mm$',...
    'Sim; Holes; $U_d = 2.5 mm$','Sim; Holes; $U_d = 5 mm$','Sim; No Holes; $U_d = 2.5 mm$',...
    'Sim; No Holes; $U_d = 5 mm$', 'Sim; No Holes; $U_d = 7 mm$','Exp; Holes; $U_d = 2.5 mm$',...
    'Exp; Holes; $U_d = 5 mm$','Exp; No Holes; $U_d = 2.5 mm$',...
    'Exp; No Holes; $U_d = 5 mm$', 'Exp; No Holes; $U_d = 7 mm$','interpreter','latex')
saveas(gcf,'sens_metric_all.pdf')
saveas(gcf,'sens_metric_all.png')

%% Other

% a = [-4.94,-4.22,-3.87,-2,-1.84,-1.33,1,1.45,1.5,1.65,1.89,1.91,2,2.86,3.03,3.3];
% 
% for i = 1:length(a)
%     idx = round((a(i) + 5.01)*100);
%     figure; imagesc(squeeze(Salpha_all(:,idx,:))') % a = 2
%     set(gca,'YDir','normal'); axis image;
%     title(['$\alpha=$',num2str(a(i))],'interpreter','latex')
%     ylabel('k'); xlabel('$\lambda$','interpreter','latex')
%     axis square;
%     saveas(gcf,['k_lam_SensitivityPlots_uniaxial/sens_' num2str(i) '.png']);
% end

% figure; imagesc(squeeze(Salpha_all(:,501,:))') % a = 2
% set(gca,'YDir','normal'); axis image;
% axis square;
% 
% 
% figure; imagesc(squeeze(Salpha_all(:,1,:))') % a = -3
% set(gca,'YDir','normal'); axis image;
% axis square;