%% Script to categorize the cell responses to odor stimulation
% Please refere to Mazo et al., Nat Comm 2022
% srcipt used to analyze and plot data presented in Fig. 6e (MC and TC) and
% Supplementary Fig 8d (JuxtaG)
% written by Camille Mazo

clearvars -except data TCexp GLexp
zscore_th = 1.96; % i.e, alpha = 0.05;
tAna = 135:150;
tBase = 90:105;

MCodor = cat(3,data.CL.odorA,data.Pulsed.odorA);
MCodor(:,:,:,2) = cat(3,data.CL.odorB,data.Pulsed.odorB);

TCodor = cat(3,TCexp.CL.odorA,TCexp.Pulsed.odorA);
TCodor(:,:,:,2) = cat(3,TCexp.CL.odorB,TCexp.Pulsed.odorB);

GLodor = GLexp.CL.odorA;
GLodor(:,:,:,2) = GLexp.CL.odorB;

for cellType = 1:3
    if cellType == 1
        data_toAnalyze = MCodor;
    elseif cellType == 2
        data_toAnalyze = TCodor;
    else
        data_toAnalyze = GLodor;
    end
    nCells = size(data_toAnalyze,3);
    
    %% z-score
    mu_resp = squeeze(nanmean(nanmean(data_toAnalyze(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(data_toAnalyze(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(data_toAnalyze(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(data_toAnalyze(:,tBase,:,:),2),[],1,'omitnan'));
    
    TrialsRESP = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(data_toAnalyze,1) + sigma_base.^2/size(data_toAnalyze,1));
    TrialsBASE = (mu_base - mu_base)./sqrt(sigma_base.^2/size(data_toAnalyze,1) + sigma_base.^2/size(data_toAnalyze,1));
    
    %% test
    h_odor = false(2,nCells);
    neg = false(2,nCells);
    for i = 1:nCells
        for ii = 1:2
            respDist = TrialsRESP(i,ii);
            baseDist = TrialsBASE(i,ii);
            if abs(respDist) > zscore_th
                h_odor(ii,i) = true;
                if respDist<0
                    neg(ii,i) = true;
                end
            end
        end
    end
    
    resp_all{cellType} = cat(1,TrialsRESP(:,1),TrialsRESP(:,2));
    h_all{cellType} = cat(2,h_odor(1,:),h_odor(2,:));
    neg_all{cellType} = cat(2,neg(1,:),neg(2,:));
end

%% stacked histo
figure; hold on
MCrespCat = [sum(h_all{1}&neg_all{1}) sum(~h_all{1}) sum(h_all{1}&~neg_all{1})]./length(h_all{1});
TCrespCat = [sum(h_all{2}&neg_all{2}) sum(~h_all{2}) sum(h_all{2}&~neg_all{2})]./length(h_all{2});
b = bar([MCrespCat;TCrespCat],'stacked');

figure;
GLrespCat = [sum(h_all{3}&neg_all{3}) sum(~h_all{3}) sum(h_all{3}&~neg_all{3})]./length(h_all{3});
b = bar([GLrespCat;NaN(1,3)],'stacked');

%% violin plot
figure;
for cellType = 1:2
    nInhib = sum(h_all{cellType}&neg_all{cellType}) ;
    nExc = sum(h_all{cellType}&~neg_all{cellType});
    subplot(1,2,cellType);hold on
    scatter(0.8*ones(nInhib,1),resp_all{cellType}(h_all{cellType}&neg_all{cellType},1),'filled','markerfacealpha',0.2);
    scatter(1.8*ones(nExc,1),resp_all{cellType}(h_all{cellType}&~neg_all{cellType},1),'filled','markerfacealpha',0.2);
    violin_data{:,1} = resp_all{cellType}(h_all{cellType}&neg_all{cellType});
    violin_data{:,2} = resp_all{cellType}(h_all{cellType}&~neg_all{cellType});
    violin_plot(violin_data);
    ylim([-20 30])
    xticks([1:2]);xticklabels({'inhibited','excited'});xtickangle(45)
    if cellType == 1
        ylabel('odor responses (z-score)')
        title('MC')
    else
        title('TC')
    end
end

figure; hold on
nInhib = sum(h_all{3}&neg_all{3}) ;
nExc = sum(h_all{3}&~neg_all{3});
scatter(0.8*ones(nInhib,1),resp_all{3}(h_all{3}&neg_all{3},1),'filled','markerfacealpha',0.2);
scatter(1.8*ones(nExc,1),resp_all{3}(h_all{3}&~neg_all{3},1),'filled','markerfacealpha',0.2);
violin_data{:,1} = resp_all{3}(h_all{3}&neg_all{3});
violin_data{:,2} = resp_all{3}(h_all{3}&~neg_all{3});
violin_plot(violin_data);
ylim([-10 15])
xticks([1:2]);xticklabels({'inhibited','excited'});xtickangle(45)
ylabel('odor responses (z-score)')
title('Juxta')