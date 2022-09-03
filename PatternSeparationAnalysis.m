%% Script to calculate the euclidean distance 
% Please refere to Mazo et al., Nat Comm 2022
% srcipt used to analyze and plot data presented in Fig. 6g
% written by Camille Mazo

% load the data with the session info: MC_allSessions and TC_allSessions

% clearvars -except MCCELLS_Oct19 TCCELLS
%% Preambule
tAna = 135:150;
tBase = 90:105;

alpha_odor = 0.01;

selec_odor = 1; % 1: OR, 2: AND, 0: none

%% Stim Time-averaged and Trial-averaged response for all cells, odor-by-odor
popMC_odorA = cell(1,length(MCCELLS_Oct19.odor));
popMC_odorB = cell(1,length(MCCELLS_Oct19.odor));
popMC_odorA_light=  cell(1,length(MCCELLS_Oct19.odor));
popMC_odorB_light = cell(1,length(MCCELLS_Oct19.odor));
popTC_odorA = cell(1,length(TCCELLS.odor));
popTC_odorB = cell(1,length(TCCELLS.odor));
popTC_odorA_light = cell(1,length(TCCELLS.odor));
popTC_odorB_light = cell(1,length(TCCELLS.odor));

for j = 1:length(MCCELLS_Oct19.odor)
    mu_resp = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor(j).OdorA_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor(j).OdorA_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(MCCELLS_Oct19.odor(j).OdorA_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(MCCELLS_Oct19.odor(j).OdorA_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popMC_odorA{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(MCCELLS_Oct19.odor(j).OdorA_allTrials,1) + sigma_base.^2/size(MCCELLS_Oct19.odor(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
    mu_resp = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor(j).OdorB_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor(j).OdorB_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(MCCELLS_Oct19.odor(j).OdorB_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(MCCELLS_Oct19.odor(j).OdorB_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popMC_odorB{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(MCCELLS_Oct19.odor(j).OdorA_allTrials,1) + sigma_base.^2/size(MCCELLS_Oct19.odor(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
    mu_resp = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor_light(j).OdorA_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor_light(j).OdorA_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(MCCELLS_Oct19.odor_light(j).OdorA_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(MCCELLS_Oct19.odor_light(j).OdorA_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popMC_odorA_light{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(MCCELLS_Oct19.odor_light(j).OdorA_allTrials,1) + sigma_base.^2/size(MCCELLS_Oct19.odor_light(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
    mu_resp = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor_light(j).OdorB_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(MCCELLS_Oct19.odor_light(j).OdorB_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(MCCELLS_Oct19.odor_light(j).OdorB_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(MCCELLS_Oct19.odor_light(j).OdorB_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popMC_odorB_light{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(MCCELLS_Oct19.odor_light(j).OdorA_allTrials,1) + sigma_base.^2/size(MCCELLS_Oct19.odor_light(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
end
for j = 1:length(TCCELLS.odor)
    mu_resp = squeeze(nanmean(nanmean(TCCELLS.odor(j).OdorA_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(TCCELLS.odor(j).OdorA_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(TCCELLS.odor(j).OdorA_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(TCCELLS.odor(j).OdorA_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popTC_odorA{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(TCCELLS.odor(j).OdorA_allTrials,1) + sigma_base.^2/size(TCCELLS.odor(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
    mu_resp = squeeze(nanmean(nanmean(TCCELLS.odor(j).OdorB_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(TCCELLS.odor(j).OdorB_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(TCCELLS.odor(j).OdorB_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(TCCELLS.odor(j).OdorB_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popTC_odorB{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(TCCELLS.odor(j).OdorA_allTrials,1) + sigma_base.^2/size(TCCELLS.odor(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
    mu_resp = squeeze(nanmean(nanmean(TCCELLS.odor_light(j).OdorA_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(TCCELLS.odor_light(j).OdorA_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(TCCELLS.odor_light(j).OdorA_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(TCCELLS.odor_light(j).OdorA_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popTC_odorA_light{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(TCCELLS.odor_light(j).OdorA_allTrials,1) + sigma_base.^2/size(TCCELLS.odor_light(j).OdorA_allTrials,1));
    clear mu_resp mu_base sigma_resp sigma_base
    mu_resp = squeeze(nanmean(nanmean(TCCELLS.odor_light(j).OdorB_allTrials(:,tAna,:,:),1),2));
    mu_base = squeeze(nanmean(nanmean(TCCELLS.odor_light(j).OdorB_allTrials(:,tBase,:,:),1),2));
    sigma_resp = squeeze(std(nanmean(TCCELLS.odor_light(j).OdorB_allTrials(:,tAna,:,:),2),[],1,'omitnan'));
    sigma_base = squeeze(std(nanmean(TCCELLS.odor_light(j).OdorB_allTrials(:,tBase,:,:),2),[],1,'omitnan'));
    popTC_odorB_light{j} = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(TCCELLS.odor_light(j).OdorA_allTrials,1) + sigma_base.^2/size(TCCELLS.odor_light(j).OdorA_allTrials,1));
end


%% Select odor responsive cells only
if alpha_odor == 0.05
    zscore_th = 1.96;
elseif alpha_odor == 0.01
    zscore_th = 2.58;
end

h_odorA_MC = cell(size(popMC_odorA));
h_odorB_MC = cell(size(popMC_odorB));
for j = 1:size(popMC_odorA,2)
    h_odorA_MC{1,j} = false(size(popMC_odorA{1,j}));
    h_odorB_MC{1,j} = false(size(popMC_odorB{1,j}));
    for i = 1:length(popMC_odorA{j})
        if popMC_odorA{1,j}(i)<-zscore_th || popMC_odorA{1,j}(i) > zscore_th
            h_odorA_MC{1,j}(i) = true;
        end
        if popMC_odorB{1,j}(i)<-zscore_th || popMC_odorB{1,j}(i) > zscore_th
            h_odorB_MC{1,j}(i) = true;
        end
    end
end
h_odorA_TC = cell(size(popTC_odorA));
h_odorB_TC = cell(size(popTC_odorB));
for j = 1:size(popTC_odorA,2)
    h_odorA_TC{1,j} = false(size(popTC_odorA{1,j}));
    h_odorB_TC{1,j} = false(size(popTC_odorB{1,j}));
    for i = 1:length(popTC_odorA{j})
        if popTC_odorA{1,j}(i)<-zscore_th || popTC_odorA{1,j}(i) > zscore_th
            h_odorA_TC{1,j}(i) = true;
        end
        if popTC_odorB{1,j}(i)<-zscore_th || popTC_odorB{1,j}(i) > zscore_th
            h_odorB_TC{1,j}(i) = true;
        end
    end
end

%% Population correlations
% Population correlation, between odors
% - MC
for j = 1:length(popMC_odorA)
    selec_odorA = h_odorA_MC{j};
    selec_odorB = h_odorB_MC{j};
    selec = false(size(selec_odorA));
    for i = 1:length(selec_odorA)
        if selec_odor == 1
            if selec_odorA(i) || selec_odorB(i)
                selec(i) = true;
            end
        elseif selec_odor == 2
            if selec_odorA(i) && selec_odorB(i)
                selec(i) = true;
            end
        elseif selec_odor == 0
            selec(i) = true;
        end
    end
    
    % - euclidean distance
    if sum(selec)>5
        temp_dist = pdist2(popMC_odorA{j}(selec)',popMC_odorB{j}(selec)','euclidean');
        PopCorr_MC_odor(j) = mean2(temp_dist);
        temp_dist2 = pdist2(popMC_odorA_light{j}(selec)',popMC_odorB_light{j}(selec)','euclidean');
        PopCorr_MC_odor_light(j) = mean2(temp_dist2);
    end
    
end

% - TC
for j = 1:length(TCCELLS.odor)
    selec_odorA = h_odorA_TC{j};
    selec_odorB = h_odorB_TC{j};
    selec = false(size(selec_odorA));
    for i = 1:length(selec_odorA)
        if selec_odor == 1
            if selec_odorA(i) || selec_odorB(i)
                selec(i) = true;
            end
        elseif selec_odor == 2
            if selec_odorA(i) && selec_odorB(i)
                selec(i) = true;
            end
        elseif selec_odor == 0
            selec(i) = true;
        end
    end
    
    % - euclidean distance
    if sum(selec)>5
        temp = pdist2(popTC_odorA{j}(selec)',popTC_odorB{j}(selec)','euclidean');
        PopCorr_TC_odor(j) = mean2(temp);
        temp2 = pdist2(popTC_odorA_light{j}(selec)',popTC_odorB_light{j}(selec)','euclidean');
        PopCorr_TC_odor_light(j) = mean(temp2);
        
    end
end

% plot
NaN_MCodor = sum(isnan(PopCorr_MC_odor));
NaN_MCodor_light = sum(isnan(PopCorr_MC_odor_light));
PopCorr_TC_odor(find(~PopCorr_TC_odor))=NaN;
NaN_TCodor = sum(isnan(PopCorr_TC_odor));
PopCorr_TC_odor_light(find(~PopCorr_TC_odor_light))=NaN;
NaN_TCodor_light = sum(isnan(PopCorr_TC_odor_light));
sem_MCodor = std(PopCorr_MC_odor,'omitnan')/sqrt(size(PopCorr_MC_odor,2)-NaN_MCodor);
sem_MCodor_light = std(PopCorr_MC_odor_light,'omitnan')/sqrt(size(PopCorr_MC_odor_light,2)-NaN_MCodor_light);
sem_TCodor = std(PopCorr_TC_odor,'omitnan')/sqrt(size(PopCorr_TC_odor,2)-NaN_TCodor);
sem_TCodor_light = std(PopCorr_TC_odor_light,'omitnan')/sqrt(size(PopCorr_TC_odor_light,2)-NaN_TCodor_light);

figure;
subplot(2,2,1); hold on
bar([1 2],[nanmean(PopCorr_MC_odor) nanmean(PopCorr_MC_odor_light)])
plot([1 1],[nanmean(PopCorr_MC_odor)-sem_MCodor nanmean(PopCorr_MC_odor)+sem_MCodor],'k-')
plot([2 2],[nanmean(PopCorr_MC_odor_light)-sem_MCodor_light nanmean(PopCorr_MC_odor_light)+sem_MCodor_light],'k-')
for i = 1:length(PopCorr_MC_odor)
    plot([0.8 2.2],[PopCorr_MC_odor(i) PopCorr_MC_odor_light(i)],'k-')
end
ylabel('euclidean distance')

xticks([1 2])
xticklabels({'Odor only','Odor+light'})
xtickangle(45)
title('MCs')
subplot(2,2,2); hold on
bar([1 2],[nanmean(PopCorr_TC_odor) nanmean(PopCorr_TC_odor_light)])
plot([1 1],[nanmean(PopCorr_TC_odor)-sem_TCodor nanmean(PopCorr_TC_odor)+sem_TCodor],'k-')
plot([2 2],[nanmean(PopCorr_TC_odor_light)-sem_TCodor_light nanmean(PopCorr_TC_odor_light)+sem_TCodor_light ],'k-')
for i = 1:length(PopCorr_TC_odor)
    plot([0.8 2.2],[PopCorr_TC_odor(i) PopCorr_TC_odor_light(i)],'k-')
end
title('TCs')
xticks([1 2])
xticklabels({'Odor only','Odor+light'})
xtickangle(45)

% stats
if ~adtest(PopCorr_MC_odor) && ~adtest(PopCorr_MC_odor_light)
    [~,p2]=ttest(PopCorr_MC_odor,PopCorr_MC_odor_light);
    testName = 'ttest';
else
    [p2,~]=signrank(PopCorr_MC_odor,PopCorr_MC_odor_light);
    testName = 'Wilcoxon';
end
if ~adtest(PopCorr_TC_odor) && ~adtest(PopCorr_TC_odor_light)
    [~,p3]=ttest(PopCorr_TC_odor,PopCorr_TC_odor_light);
    testName = 'ttest';
else
    [p3,~]=signrank(PopCorr_TC_odor,PopCorr_TC_odor_light);
    testName = 'Wilcoxon';
end

% type info
subplot(2,2,3);
title([testName ', p = ', num2str(p2)]);
xl = xlim; yl = ylim;
text(xl(1),yl(2),[{'mean +- sem'},...
    {['odor: ', num2str(mean(PopCorr_MC_odor),3),'+-',num2str(sem_MCodor,3)]},...
    {['odor+light: ', num2str(mean(PopCorr_MC_odor_light),3),'+-',num2str(sem_MCodor_light,3)]},...
    {['n = ',num2str(length(PopCorr_MC_odor))]}],...
    'HorizontalAlignment','left','VerticalAlignment','top');
axis off
subplot(2,2,4);
xl = xlim; yl = ylim;
text(xl(1),yl(2),[{'mean +- sem'},...
    {['odor: ', num2str(nanmean(PopCorr_TC_odor),3),'+-',num2str(sem_TCodor,3)]},...
    {['odor+light: ', num2str(nanmean(PopCorr_TC_odor_light),3),'+-',num2str(sem_TCodor_light,3)]},...
    {['n = ',num2str(length(PopCorr_TC_odor)-NaN_TCodor)]}],...
    'HorizontalAlignment','left','VerticalAlignment','top');
axis off
title([testName ', p = ', num2str(p3)]);
text(xl(1),yl(1),['selection: ',num2str(selec_odor)],...
    'HorizontalAlignment','right','VerticalAlignment','top');


%% Pairwise correlations, average across all cells for a given odor, for every odor
figure;
for celltype = 1:2 %1=MC, 2=TC
    if celltype == 1
        AllResp_odor = cell(2,length(MCCELLS_Oct19.odor)); AllResp_odorlight = cell(2,length(MCCELLS_Oct19.odor));
    elseif celltype == 2
        AllResp_odor = cell(2,length(TCCELLS.odor)); AllResp_odorlight = cell(2,length(TCCELLS.odor));
    end
    RinterCells_odor = cell(size(AllResp_odor)); RinterCells_odorlight = cell(size(AllResp_odor));
    clear  mean_corr_odor mean_corr_odorlight GrandMean_corr sem
    if celltype == 1
        data_odor = MCCELLS_Oct19.odor;
        data_odorlight = MCCELLS_Oct19.odor_light;
    else
        data_odor = TCCELLS.odor;
        data_odorlight = TCCELLS.odor_light;
    end
    
    for j = 1:length(data_odor)
        if selec_odor == 1 || selec_odor == 2 % odor by odor, not pairwise
            if celltype == 1
                selec_odorA = h_odorA_MC{j};
                selec_odorB = h_odorB_MC{j};
            else
                selec_odorA = h_odorA_TC{j};
                selec_odorB = h_odorB_TC{j};
            end
        elseif selec_odor == 0
            if celltype == 1
                selec_odorA = true(size(h_odorA_MC{j}));
                selec_odorB = true(size(h_odorB_MC{j}));
            elseif celltype == 2
                selec_odorA = true(size(h_odorA_TC{j}));
                selec_odorB = true(size(h_odorB_TC{j}));
            end
        end
        
        if sum(selec_odorA)>5
            AllResp_odor{1,j} = squeeze(nanmean(data_odor(j).OdorA_allTrials(:,tAna,selec_odorA),2));
            AllResp_odorlight{1,j} = squeeze(nanmean(data_odorlight(j).OdorA_allTrials(:,tAna,selec_odorA),2));
        end
        
        if sum(selec_odorB)>5
            AllResp_odor{2,j} = squeeze(nanmean(data_odor(j).OdorB_allTrials(:,tAna,selec_odorB),2));
            AllResp_odorlight{2,j} = squeeze(nanmean(data_odorlight(j).OdorB_allTrials(:,tAna,selec_odorB),2));
            
        end
    end
    
    for j = 1:length(data_odor)
        for jj = 1:2
            if ~isempty(AllResp_odorlight{jj,j})
                for i = 30:-1:1
                    if isnan(AllResp_odorlight{jj,j}(i,1))
                        AllResp_odorlight{jj,j}(i,:)=[];
                    end
                end
            end
            if ~isempty(AllResp_odor{jj,j})
                for i = 30:-1:1
                    if isnan(AllResp_odor{jj,j}(i,1))
                        AllResp_odor{jj,j}(i,:)=[];
                    end
                end
            end
        end
    end
    % euclidean distance
    for j = 1:size(AllResp_odor,1)
        for jj = 1:size(AllResp_odor,2)
            RinterCells_odor{j,jj}      = pdist2(AllResp_odor{j,jj}',AllResp_odor{j,jj}','euclidean');
            RinterCells_odorlight{j,jj} = pdist2(AllResp_odorlight{j,jj}',AllResp_odorlight{j,jj}','euclidean');            
        end
    end
    for j = 1:size(AllResp_odor,1)
        for jj = 1:size(AllResp_odor,2)
            if ~isempty(RinterCells_odor{j,jj}) && size(RinterCells_odor{j,jj},2)>5
                for i = 1:size(AllResp_odor{j,jj},2)
                    RinterCells_odor{j,jj}(i,i) = NaN;
                    RinterCells_odorlight{j,jj}(i,i) = NaN;
                    mean_corr_odor(j,jj) = nanmean(nanmean(RinterCells_odor{j,jj},1),2);
                    mean_corr_odorlight(j,jj) = nanmean(nanmean(RinterCells_odorlight{j,jj},1),2);
                end
            end
        end
    end
    
    mean_corr_odor = cat(2,mean_corr_odor(1,:),mean_corr_odor(2,:));
    mean_corr_odor(find(~mean_corr_odor))= NaN;
    mean_corr_odorlight = cat(2,mean_corr_odorlight(1,:),mean_corr_odorlight(2,:));
    mean_corr_odorlight(find(~mean_corr_odorlight))= NaN;
    GrandMean_corr(celltype,1) = nanmean(mean_corr_odor');
    GrandMean_corr(celltype,2) = nanmean(mean_corr_odorlight');
    sem(celltype,1)= std(mean_corr_odor','omitnan')./sqrt(size(mean_corr_odor,2)-sum(isnan(mean_corr_odor)));
    sem(celltype,2)= std(mean_corr_odorlight','omitnan')./sqrt(size(mean_corr_odorlight,2)-sum(isnan(mean_corr_odorlight)));
    
    % stats
    if ~adtest(mean_corr_odor) && ~adtest(mean_corr_odorlight)
        [~,p5]=ttest(mean_corr_odor,mean_corr_odorlight);
        testName = 'ttest';
    else
        [p5,~]=ranksum(mean_corr_odor,mean_corr_odorlight);
        testName = 'Wilcoxon';
    end
    
    % plot
    subplot(2,2,celltype); hold on
    bar(GrandMean_corr(celltype,:))
    plot([1 1],[GrandMean_corr(celltype,1)+sem(celltype,1) GrandMean_corr(celltype,1)-sem(celltype,1)],'k-')
    plot([2 2],[GrandMean_corr(celltype,2)+sem(celltype,2) GrandMean_corr(celltype,2)-sem(celltype,2)],'k-')
    
    % finalize plotting
    if celltype == 1
        title('MCs')
    else
        title('TCs')
    end
    for i = 1:length(mean_corr_odor)
        plot([0.8 2.2],[mean_corr_odor(i) mean_corr_odorlight(i)],'k-')
    end
    xticks([1 2])
    xticklabels({'odor','odor+light'})
    xtickangle(45)
    ylabel('euclidean distance')
    
    subplot(2,2,celltype+2)
    title([testName, ', p = ' num2str(p5)])
    xl=xlim; yl=ylim;
    text(xl(1),yl(2),[{'mean +- sem'},...
        {['odor: ',num2str(GrandMean_corr(celltype,1),3),'+-',num2str(sem(celltype,1),3)]},...
        {['odor+light: ',num2str(GrandMean_corr(celltype,2),3),'+-',num2str(sem(celltype,2),3)]},...
        {['n = ',num2str(size(mean_corr_odorlight,2)-sum(isnan(mean_corr_odor))),' pairwise pop. average-odor pairs']}],...
        'HorizontalAlignment','left','VerticalAlignment','top');
    axis off
end
