%% Script to analyze the cell responses to optogentic stimulation of GABAergic feedback on spontaneous activity
% Please refere to Mazo et al., Nat Comm 2022
% srcipt used to analyze and plot data presented in Fig. 6d and
% supplementary figure 8d
% written by Camille Mazo

%% load the data
% load 'MCexp','MCctrl','MCshutterCtl','GLexp','GLctrl','TCctrl','TCex'

%% Windows for response and baseline, in frames (acquisition at 30 fps)
tAna = 135:150;
tBase = 90:105;
% z-score = (mu,resp - mu,baseline) / (sqrt [ (sigma,resp)^{2}/n + (sigma,baseline)^{2}/n ]).

%% Calculate the z-scores
% - MC
MC_cells = data.CL.light;
Exp_resp_mu = squeeze(nanmean(nanmean(MC_cells(:,tAna,:),1),2));
Exp_base_mu = squeeze(nanmean(nanmean(MC_cells(:,tBase,:),1),2));
Exp_base_sd = squeeze(std(nanmean(MC_cells(:,tBase,:),2),[],1,'omitnan'));
Exp_resp_sd = squeeze(std(nanmean(MC_cells(:,tAna,:),2),[],1,'omitnan'));
z_exp_CL = (Exp_resp_mu - Exp_base_mu)./sqrt(Exp_resp_sd.^2/size(MC_cells,1)+Exp_base_sd.^2/size(MC_cells,1));

MC_cells = data.Pulsed.light;
Exp_resp_mu = squeeze(nanmean(nanmean(MC_cells(:,tAna,:),1),2));
Exp_base_mu = squeeze(nanmean(nanmean(MC_cells(:,tBase,:),1),2));
Exp_base_sd = squeeze(std(nanmean(MC_cells(:,tBase,:),2),[],1,'omitnan'));
Exp_resp_sd = squeeze(std(nanmean(MC_cells(:,tAna,:),2),[],1,'omitnan'));
z_exp_Pulsed = (Exp_resp_mu - Exp_base_mu)./sqrt(Exp_resp_sd.^2/size(MC_cells,1)+Exp_base_sd.^2/size(MC_cells,1));

% - light controls
LightCtl = cat(3,MClightControl{1,:});
LightCtrl_resp_mu = squeeze(nanmean(nanmean(LightCtl(:,tAna,:),2),1));
LightCtrl_base_mu = squeeze(nanmean(nanmean(LightCtl(:,tBase,:),2),1));
LightCtrl_resp_sigma = squeeze(std(nanmean(LightCtl(:,tAna,:),2),[],1,'omitnan'));
LightCtrl_base_sigma = squeeze(std(nanmean(LightCtl(:,tBase,:),2),[],1,'omitnan'));
z_LightCtrl_CL = (LightCtrl_resp_mu - LightCtrl_base_mu)./sqrt(LightCtrl_resp_sigma.^2/size(LightCtl,1)+LightCtrl_base_sigma.^2/size(LightCtl,1));

LightCtl = cat(3,MClightControl{2,:});
LightCtrl_resp_mu = squeeze(nanmean(nanmean(LightCtl(:,tAna,:),2),1));
LightCtrl_base_mu = squeeze(nanmean(nanmean(LightCtl(:,tBase,:),2),1));
LightCtrl_resp_sigma = squeeze(std(nanmean(LightCtl(:,tAna,:),2),[],1,'omitnan'));
LightCtrl_base_sigma = squeeze(std(nanmean(LightCtl(:,tBase,:),2),[],1,'omitnan'));
z_LightCtrl_Pulsed = (LightCtrl_resp_mu - LightCtrl_base_mu)./sqrt(LightCtrl_resp_sigma.^2/size(LightCtl,1)+LightCtrl_base_sigma.^2/size(LightCtl,1));

% Shutter control
Shutter = cat(3,MC_ShutterControl{:});
meanTrial_Shutter_resp = squeeze(nansum(Shutter(:,tAna,:),2));
meanTrial_Shutter_base = squeeze(nansum(Shutter(:,tBase,:),2));
Shutter_mean = squeeze(nanmean(nanmean(Shutter(:,tAna,:),2),1));
Shutter_base = squeeze(nanmean(nanmean(Shutter(:,tBase,:),2),1));
Shutter_base_sigma = std(meanTrial_Shutter_base,[],1,'omitnan');
Shutter_resp_sigma = std(meanTrial_Shutter_resp,[],1,'omitnan');
z_shutter = (Shutter_mean - Shutter_base)./sqrt(Shutter_base_sigma'.^2/size(Shutter,1)+Shutter_resp_sigma'.^2/size(Shutter,1));

% - TC 
TC_cells = TCexp.CL.light;
Exp_resp_mu = squeeze(nanmean(nanmean(TC_cells(:,tAna,:),1),2));
Exp_base_mu = squeeze(nanmean(nanmean(TC_cells(:,tBase,:),1),2));
Exp_base_sd = squeeze(std(nanmean(TC_cells(:,tBase,:),2),[],1,'omitnan'));
Exp_resp_sd = squeeze(std(nanmean(TC_cells(:,tAna,:),2),[],1,'omitnan'));
z_TC_exp_CL = (Exp_resp_mu - Exp_base_mu)./sqrt(Exp_resp_sd.^2/size(TC_cells,1)+Exp_base_sd.^2/size(TC_cells,1));

TC_cells = TCexp.Pulsed.light;
Exp_resp_mu = squeeze(nanmean(nanmean(TC_cells(:,tAna,:),1),2));
Exp_base_mu = squeeze(nanmean(nanmean(TC_cells(:,tBase,:),1),2));
Exp_base_sd = squeeze(std(nanmean(TC_cells(:,tBase,:),2),[],1,'omitnan'));
Exp_resp_sd = squeeze(std(nanmean(TC_cells(:,tAna,:),2),[],1,'omitnan'));
z_TC_exp_Pulsed = (Exp_resp_mu - Exp_base_mu)./sqrt(Exp_resp_sd.^2/size(TC_cells,1)+Exp_base_sd.^2/size(TC_cells,1));

% - TC light control
LightCtl = cat(3,TClightControl{1,:});
LightCtrl_resp_mu = squeeze(nanmean(nanmean(LightCtl(:,tAna,:),2),1));
LightCtrl_base_mu = squeeze(nanmean(nanmean(LightCtl(:,tBase,:),2),1));
LightCtrl_resp_sigma = squeeze(std(nanmean(LightCtl(:,tAna,:),2),[],1,'omitnan'));
LightCtrl_base_sigma = squeeze(std(nanmean(LightCtl(:,tBase,:),2),[],1,'omitnan'));
z_TC_LightCtrl_CL = (LightCtrl_resp_mu - LightCtrl_base_mu)./sqrt(LightCtrl_resp_sigma.^2/size(LightCtl,1)+LightCtrl_base_sigma.^2/size(LightCtl,1));

LightCtl = cat(3,TClightControl{2,:});
LightCtrl_resp_mu = squeeze(nanmean(nanmean(LightCtl(:,tAna,:),2),1));
LightCtrl_base_mu = squeeze(nanmean(nanmean(LightCtl(:,tBase,:),2),1));
LightCtrl_resp_sigma = squeeze(std(nanmean(LightCtl(:,tAna,:),2),[],1,'omitnan'));
LightCtrl_base_sigma = squeeze(std(nanmean(LightCtl(:,tBase,:),2),[],1,'omitnan'));
z_TC_LightCtrl_Pulsed = (LightCtrl_resp_mu - LightCtrl_base_mu)./sqrt(LightCtrl_resp_sigma.^2/size(LightCtl,1)+LightCtrl_base_sigma.^2/size(LightCtl,1));


%% Main plot (Fig. 6d)
figure;hold on
plot([0 10],[0 0],'k:')
scatter(0.8*ones(size(z_shutter)),z_shutter,'filled','markerfacealpha',0.2)
scatter(1.8*ones(size(z_exp_Pulsed)),z_exp_Pulsed,'filled','markerfacealpha',0.2)
scatter(3.8*ones(size(z_exp_CL)),z_exp_CL,'filled','markerfacealpha',0.2)
scatter(2.8*ones(size(z_LightCtrl1)),z_LightCtrl1,'filled','markerfacealpha',0.2)
scatter(4.8*ones(size(z_LightCtrl2)),z_LightCtrl2,'filled','markerfacealpha',0.2)
scatter(5.8*ones(size(z_TC_exp_Pulsed)),z_TC_exp_Pulsed,'filled','markerfacealpha',0.2)
scatter(7.8*ones(size(z_TC_exp_CL)),z_TC_exp_CL,'filled','markerfacealpha',0.2)
scatter(6.8*ones(size(z_TC_LightCtrl1)),z_TC_LightCtrl1,'filled','markerfacealpha',0.2)
scatter(8.8*ones(size(z_TC_LightCtrl2)),z_TC_LightCtrl2,'filled','markerfacealpha',0.2)

size_max = 650;
z_shutter2 = [z_shutter ;NaN(size_max-length(z_shutter),1)];
z_exp_Pulsed = [z_exp_Pulsed;NaN(size_max-length(z_exp_Pulsed),1)];
z_LightCtrl1 = [z_LightCtrl_CL;NaN(size_max-length(z_LightCtrl_CL),1)];
z_LightCtrl2 = [z_LightCtrl_Pulsed;NaN(size_max-length(z_LightCtrl_Pulsed),1)];
z_TC_exp_Pulsed = [z_TC_exp_Pulsed;NaN(size_max-length(z_TC_exp_Pulsed),1)];
z_TC_exp_CL = [z_TC_exp_CL;NaN(size_max-length(z_TC_exp_CL),1)];
z_TC_LightCtrl1 = [z_TC_LightCtrl_CL;NaN(size_max-length(z_TC_LightCtrl_CL),1)];
z_TC_LightCtrl2= [z_TC_LightCtrl_Pulsed;NaN(size_max-length(z_TC_LightCtrl_Pulsed),1)];

violin_plot([z_shutter2 z_exp_Pulsed z_LightCtrl1 z_TC_exp_CL z_LightCtrl1 z_TC_exp_Pulsed z_TC_LightCtrl1 z_TC_exp_CL z_TC_LightCtrl2]);

ylabel('light responses (z-score)')
xticks([1:9]);
xticklabels({'ShutterCtrl','MC 33Hz', 'MC 33Hz Ctrl', 'MC CL',...
    'MC CL Ctrl', 'TC 33Hz', 'TC 33Hz Ctrl', 'TC CL', 'TC CL ctrl'});
xtickangle(45)

%% GL
GL_cells = GLexp.CL.light;
Exp_resp_mu = squeeze(nanmean(nanmean(GL_cells(:,tAna,:),1),2));
Exp_base_mu = squeeze(nanmean(nanmean(GL_cells(:,tBase,:),1),2));
Exp_base_sd = squeeze(std(nanmean(GL_cells(:,tBase,:),2),[],1,'omitnan'));
Exp_resp_sd = squeeze(std(nanmean(GL_cells(:,tAna,:),2),[],1,'omitnan'));
z_GL_exp = (Exp_resp_mu - Exp_base_mu)./sqrt(Exp_resp_sd.^2/size(GL_cells,1)+Exp_base_sd.^2/size(GL_cells,1));

LightCtl = cat(3,GLlightControl{2,:});
LightCtrl_resp_mu = squeeze(nanmean(nanmean(LightCtl(:,tAna,:),2),1));
LightCtrl_base_mu = squeeze(nanmean(nanmean(LightCtl(:,tBase,:),2),1));
LightCtrl_resp_sigma = squeeze(std(nanmean(LightCtl(:,tAna,:),2),[],1,'omitnan'));
LightCtrl_base_sigma = squeeze(std(nanmean(LightCtl(:,tBase,:),2),[],1,'omitnan'));
z_GL_LightCtrl = (LightCtrl_resp_mu - LightCtrl_base_mu)./sqrt(LightCtrl_resp_sigma.^2/size(LightCtl,1)+LightCtrl_base_sigma.^2/size(LightCtl,1));

z_GL_exp = [z_GL_exp ;NaN(size_max-length(z_GL_exp),1)];
z_GL_LightCtrl = [z_GL_LightCtrl ;NaN(size_max-length(z_GL_LightCtrl),1)];

% - supplementary figure 8d
figure; hold on
plot([0 3],[0 0],'k:')
scatter(0.8*ones(size(z_GL_exp)),z_GL_exp,'filled','markerfacealpha',0.2)
scatter(1.8*ones(size(z_GL_LightCtrl)),z_GL_LightCtrl,'filled','markerfacealpha',0.2)
violin_plot([z_GL_exp z_GL_LightCtrl]);
ylabel('light responses (z-score)')
xticks([1:2]);
xticklabels({'GL CL','GL CL ctrl'});
xtickangle(45)

%% one way ANOVA and post-hoc on all data
[p_all,~,stats_all]=anova1([z_shutter2 z_exp_Pulsed z_LightCtrl1 z_exp_CL z_LightCtrl2,...
    z_TC_exp_Pulsed z_TC_LightCtrl1 z_TC_exp_CL z_TC_LightCtrl2 ,...
    z_GL_exp z_GL_LightCtrl]);
c_all = multcompare(stats_all);