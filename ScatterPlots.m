%% Script to analyze the cell responses to optogentic stimulation of GABAergic feedback on odor-evoked activity
% Please refere to Mazo et al., Nat Comm 2022
% srcipt used to analyze and plot data presented in Fig. 6f (MC and TC) and
% supplementary figure 8f (JuxtaG)
% written by Camille Mazo

%%
% Inputs:
% data is the data to analyze: MC, TC or JuxtaG (=GL)
% metric is either area under the curve (=1), dF/F (=2) or z-score (=3)
% tAna and tBase are the frames to calculate the response and baseline
% selec_light is to select only light-responsive cells (=1) or not (=0).
% selec_odor, same for odor-responsive cells
% alpha_odor is the criterion for odor-responsive
% saveFig = 1 will save the fig to saveFolder

% for the analyzes in the paper, we used metric=3, tAna=[135:150], tBase =
% [90:105], selec_light = 0, selec_odor = 1, alpha_odor = 0.05
% a typical call would be: 
% ScatterPlots(data,3,2,[135:150],[90:105],0,1,0.05,0,'YourFolder')

function ScatterPlots(data,metric, LightPattern,tAna,tBase,selec_light,selec_odor,alpha_odor,saveFig,saveFolder)

%% Preambule
if LightPattern == 1
    Data = data.CL;
    ExpType = 'CL';
elseif LightPattern == 2
    Data = data.Pulsed;
    ExpType = '33Hz';
end
ExpType = [data.name '_' ExpType];

if selec_light
    selec = Data.lightSelection.selection;
    FractionPassing = sum(Data.lightSelection.h_exp)/size(selec,2);
    disp(['Select cells based on light response being significantly greater than light control | Fraction Passing = ', num2str(FractionPassing,3)])
    ExpType = [ExpType '_select'];
elseif selec_light == 0
    selec = ones(1,size(Data.light,3));
    ExpType = [ExpType '_all'];
end

saveName_full = [saveFolder, filesep ,ExpType];
disp(saveName_full)
field = fieldnames(Data);

%% metric calculation
TrialsAVG = NaN(size(Data.(field{1}),3),numel(field)-1);
TrialsBASE = NaN(size(Data.(field{1}),3),numel(field)-1);
for i = 1:numel(field)-1
      Traces = Data.(field{i});
    if metric == 2
        TrialsAVG(:,i) = squeeze(nanmean(nanmean(Traces(:,tAna,:),1),2));
        TrialsBASE(:,i)= squeeze(nanmean(nanmean(Traces(:,tBase,:),1),2));
    elseif metric == 1
        TrialsBASE(:,i)= squeeze(nanmean(nansum(data.CL.(field{i})(:,tBase,:),2),1));
        TrialsAVG(:,i) = squeeze(nanmean(nansum(data.CL.(field{i})(:,tAna,:),2),1));
    elseif metric == 3
        mu_resp = squeeze(nanmean(nanmean(Traces(:,tAna,:),1),2));
        mu_base = squeeze(nanmean(nanmean(Traces(:,tBase,:),1),2));
        sigma_resp = squeeze(std(nanmean(Traces(:,tAna,:),2),[],1,'omitnan'));
        sigma_base = squeeze(std(nanmean(Traces(:,tBase,:),2),[],1,'omitnan'));
        TrialsAVG(:,i) = (mu_resp - mu_base)./sqrt(sigma_resp.^2/size(Traces,1) + sigma_base.^2/size(Traces,1));
        TrialsBASE(:,i) = (mu_base - mu_base)./sqrt(sigma_base.^2/size(Traces,1) + sigma_base.^2/size(Traces,1));
        %  z-score = (mu,resp - mu,baseline) / (sqrt [ (sigma,resp)^{2}/n + (sigma,baseline)^{2}/n ]).
    end
end

%% select odor
if selec_odor
    h_odor = false(2,size(Data.(field{i}),3));
    neg = false(2,size(Data.(field{1}),3));
    figure;
    subplot(1,3,2);hold on
    for i = 1:2
        disp(field{i+1})
        if metric == 1 || metric == 2
            for ii = 1:size(Data.(field{i}),3)
                Resp(:,i,ii) = squeeze(nanmean(Data.(field{i+1})(:,tAna,ii),2));
                Base = squeeze(nanmean(Data.(field{i+1})(:,tBase,ii),2));
                [h,p]=ttest(Resp(:,i,ii),Base,'Alpha',alpha_odor);
                h_odor(i,ii) = h;
                p_odor(i,ii) = p;
                if h_odor(i,ii)
                    if nanmean(Resp(:,i,ii),1)<nanmean(Base)
                        neg(i,ii) = true;
                        plot(1,nanmean(Resp(:,i,ii),1),'bo')
                    else
                        plot(2,nanmean(Resp(:,i,ii),1),'ro')
                    end
                else
                    plot(3,nanmean(Resp(:,i,ii),1),'o','Color',[0.5 0.5 0.5])
                end
            end
        elseif metric == 3
            if alpha_odor == 0.05
                zscore_th = 1.96;
            elseif alpha_odor == 0.01
                zscore_th = 2.58;
            end
            h_odor(i,TrialsAVG(:,i+1) < -zscore_th|TrialsAVG(:,i+1) > zscore_th)=true;
            for ii = 1:size(Data.(field{i}),3)
                if h_odor(i,ii)
                    if TrialsAVG(ii,i+1)<0
                        neg(i,ii) = true;
                        plot(1,TrialsAVG(ii,i+1),'bo')
                    else
                        plot(2,TrialsAVG(ii,i+1),'ro')
                    end
                else
                    plot(3,TrialsAVG(ii,i+1),'o','Color',[0.5 0.5 0.5])
                end
            end
        end
    end
    
    h_odor = logical(h_odor);    neg = logical(neg);

    % - plot the odor responses
    xlim([0 4])
    plot([0.5 3.5],[0 0],'k--')
    xticks(1:3)
    if metric == 1 || metric == 2
        meanResp = squeeze(nanmean(Resp,1));
        meanResp = cat(2,meanResp(1,:),meanResp(2,:));
    elseif metric == 3
        meanResp = cat(2,TrialsAVG(:,2)',TrialsAVG(:,3)');
    end
    h_odor2 = [h_odor(1,:),h_odor(2,:)];
    neg2 = [neg(1,:),neg(2,:)];
    sem_neg = std(meanResp(neg2),[],'omitnan')/sqrt(size(meanResp(neg2),2));
    sem_pos = std(meanResp(h_odor2&~neg2),[],'omitnan')/sqrt(size(meanResp(h_odor2&~neg2),2));
    sem_no = std(meanResp(~h_odor2),[],'omitnan')/sqrt(size(meanResp(~h_odor2),2));
    
    plot([0.8 1.2],[nanmean(meanResp(neg2)) nanmean(meanResp(neg2))],'k-','LineWidth',5)
    plot([1.8 2.2],[nanmean(meanResp(h_odor2&~neg2)) nanmean(meanResp(h_odor2&~neg2))],'k-','LineWidth',5)
    plot([2.8 3.2],[nanmean(meanResp(~h_odor2)) nanmean(meanResp(~h_odor2))],'k-','LineWidth',5)
    bar(1.5,nanmean(meanResp(neg2)),'FaceColor','none','EdgeColor',[0.5 0.5 0.5])
    bar(2.5,nanmean(meanResp(h_odor2&~neg2)),'FaceColor','none','EdgeColor',[0.5 0.5 0.5])
    bar(3.5,nanmean(meanResp(~h_odor2)),'FaceColor','none','EdgeColor',[0.5 0.5 0.5])
    plot([1.5 1.5],[nanmean(meanResp(neg2))+sem_neg nanmean(meanResp(neg2))-sem_neg],'k-')
    plot([2.5 2.5],[nanmean(meanResp(h_odor2&~neg2))+sem_pos nanmean(meanResp(h_odor2&~neg2))-sem_pos],'k-')
    plot([3.5 3.5],[nanmean(meanResp(~h_odor2))+sem_no nanmean(meanResp(~h_odor2))-sem_no],'k-')
        ylim([-20 30])
        yl=ylim;
    text(1,yl(2)-0.05,[num2str(sum(neg2)),' | ',num2str(sum(neg2)/sum(h_odor2),3)],'HorizontalAlignment','Center','VerticalAlignment','top')
    text(2,yl(2),[num2str(sum(h_odor2&~neg2)),' | ',num2str(sum(h_odor2&~neg2)/sum(h_odor2),3)],'HorizontalAlignment','Center','VerticalAlignment','top')
    text(3,yl(2)-0.05,[num2str(sum(~h_odor2)),' | ',num2str(sum(~h_odor2)/length(h_odor2),3)],'HorizontalAlignment','Center','VerticalAlignment','top')
    xticklabels({'Neg','Pos','No Change'})
    xtickangle(45)
    if metric == 1
        ylabel('AUC')
    elseif metric == 2
        ylabel('dFF')
    elseif metric == 3
        ylabel('z score')
    end
    title({['Fraction cell-odor pair passing: ',num2str(mean(sum(h_odor,2)./length(h_odor)),2)],['(Cell | Fraction)']})
    
    subplot(1,3,3); hold on
    xl = xlim; yl=ylim;
    text(xl(1)+0.1,yl(2),['Inhibited: ' num2str(nanmean(meanResp(neg2)),3) ...
        '+-' num2str(std(meanResp(neg2),[],'omitnan')./sqrt(size(meanResp(neg2),2)),3)],...
        'VerticalAlignment','top')
    text(xl(1)+0.1,(yl(2)-yl(1))/2,['Excited: ' num2str(nanmean(meanResp(h_odor2&~neg2)),3)...
        '+-' num2str(std(meanResp(h_odor2&~neg2),[],'omitnan')./sqrt(size(meanResp(h_odor2&~neg2),2)),3)])
    text(xl(1)+0.1,yl(1),['No Change: ' num2str(nanmean(meanResp(~h_odor2)),3)...
        '+-' num2str(std(meanResp(~h_odor2),[],'omitnan')./sqrt(size(meanResp(~h_odor2),2)),3)],...
        'VerticalAlignment','bottom')
    title('Mean +- sem')
    
    subplot(1,3,1);hold on
    violin_data{:,1} = meanResp(neg2);
    violin_data{:,2} = meanResp(h_odor2&~neg2);
    violin_data{:,3} = meanResp(~h_odor2);
    violin_plot(violin_data);
    clear violin_data
    xticks(1:3)
    xticklabels({'Neg','Pos','No Change'})
    xtickangle(45)
    plot([0.5 3.5],[0 0],'k--')
    ylim([-20 30])
    title({'Odor responses',['Total cell-odor pair, n = ',num2str(length(h_odor2))]})
    
    set(gcf,'Units','Normalized','OuterPosition',[0.05 0.05 0.7 0.7])
    if saveFig
        %         saveas(gcf,[saveName_full,'OdorResp.m']);
        screen2png([saveName_full,'OdorResp']);
        set(gcf,'Renderer','Painter')
        hgexport(gcf,[saveName_full,'OdorResp.eps']);
    end
else
    h_odor = true(2,size(Data.(field{i}),3));
    h_odor2 = [h_odor(1,:),h_odor(2,:)];
    disp('no selection on odor responses');
end


selec2 = logical(repmat(selec,1,2));
%% 1) x = light, y = delta("odor+light"-"odor")
% plot the scatter plot as for Supplementary Fig 8c,f
deltaAUC_odor(:,1) = TrialsAVG(:,4)-TrialsAVG(:,2);
deltaAUC_odor(:,2) = TrialsAVG(:,5)-TrialsAVG(:,3);

figure; hold on
s1 = scatter(TrialsAVG(selec&h_odor(1,:),1),deltaAUC_odor(selec&h_odor(1,:),1),'b','filled','markerfacealpha',0.2);
s2 = scatter(TrialsAVG(selec&h_odor(2,:),1),deltaAUC_odor(selec&h_odor(2,:),2),'r','filled','markerfacealpha',0.2);
xl = xlim; yl = ylim;
axis_max = max(xl(2),yl(2));
axis_min = min(xl(1),yl(1));
axis equal
xlim([axis_min axis_max]);
ylim([axis_min axis_max]);
xlabel('Light, spont');
ylabel('Light, odor-evoked')
plot([xl(1) xl(2)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
plot([xl(1) xl(2)],[xl(1) xl(2)],'k');

% correlation 
a=[TrialsAVG(selec&h_odor(1,:),1);TrialsAVG(selec&h_odor(2,:),1)];
b=[deltaAUC_odor(selec&h_odor(1,:),1);deltaAUC_odor(selec&h_odor(2,:),2)];
rh0=corr(a,b,'type','Spearman');

% mean and error
mean_CL = nanmean(a);
sem_CL = std(a);%./sqrt(size(deltaAUC(selec2&h_odor2,1),1));
mean_delta = nanmean(b);
sem_delta = std(b);%/sqrt(size(deltaAUC(selec2&h_odor2,2),1));
plot([mean_CL+sem_CL mean_CL-sem_CL],[mean_delta mean_delta],'-','color',[1 1 1])
plot([mean_CL mean_CL],[mean_delta+sem_delta mean_delta-sem_delta],'-','color',[1 1 1])

[~,p]=ttest(a,b,'alpha',0.05);

% finalize plot
axis equal;
xlim([-30 10]); ylim([-30 10])
title({ExpType, ['spearman R=' num2str(rh0,3)],...
    ['p=' num2str(p,3) ' | n=',num2str(length(a)),' cell-odor pairs']},...
    'Interpreter','none')
if saveFig
    %     saveas(gcf,[saveName_full,'LightVsdeltaAUC.m']);
    screen2png([saveName_full,'LightVsdeltaAUC']);
    set(gcf,'Renderer','Painter')
    hgexport(gcf,[saveName_full,'LightVsdeltaAUC.eps']);
end

%% 2) x = odor, y = odor+light
% plot the scatter plot as for Fig. 6f (MC, TC) or Supplementary Fig 8f (GL)
figure; hold on
scatter(TrialsAVG(selec&h_odor(1,:),2),TrialsAVG(selec&h_odor(1,:),4),'b','filled','markerfacealpha',0.2);
scatter(TrialsAVG(selec&h_odor(2,:),3),TrialsAVG(selec&h_odor(2,:),5),'r','filled','markerfacealpha',0.2);
xl = xlim; yl = ylim;
axis_max = max(xl(2),yl(2));
axis_min = min(xl(1),yl(1));
axis equal
xlim([axis_min axis_max]);
ylim([axis_min axis_max]);
plot([xl(1) xl(2)],[0 0],'k--')
plot([0 0],[yl(1) yl(2)],'k--')
plot([xl(1) xl(2)],[xl(1) xl(2)],'k');
xlabel('odor');ylabel('odor+light');

% - linear fit
odor = cat(1,TrialsAVG(:,2),TrialsAVG(:,3));
odor_CL = cat(1,TrialsAVG(:,4),TrialsAVG(:,5));

k = find(odor>0);
pos = false(1,length(odor));
pos(k) = true;

MCmdl1= fitlm(odor(selec2&h_odor2&~pos),odor_CL(selec2&h_odor2&~pos),'linear');
[MCpred1, MCci1] = predict(MCmdl1,odor(selec2&h_odor2&~pos));
plot(odor(selec2&h_odor2&~pos),MCpred1,'k-');
plot(odor(selec2&h_odor2&~pos),MCci1, 'color',[0.5 0.5 0.5],'LineStyle',':');

MCmdl2= fitlm(odor(selec2&h_odor2&pos),odor_CL(selec2&h_odor2&pos),'linear');
[MCpred2, MCci2] = predict(MCmdl2,odor(selec2&h_odor2&pos));
plot(odor(selec2&h_odor2&pos),MCpred2,'k-');
plot(odor(selec2&h_odor2&pos),MCci2, 'color',[0.5 0.5 0.5],'LineStyle',':');

% - plot the mean and error
mean_odorLight_neg = nanmean(odor_CL(selec2&h_odor2&~pos));
mean_odorLight_pos = nanmean(odor_CL(selec2&h_odor2&pos));
sem_odorLight_neg = std(odor_CL(selec2&h_odor2&~pos));%./sqrt(length(odor_CL(selec2&h_odor2&~pos)));
sem_odorLight_pos = std(odor_CL(selec2&h_odor2&pos));%./sqrt(length(odor_CL(selec2&h_odor2&pos)));
mean_odor_neg = nanmean(odor(selec2&h_odor2&~pos));
mean_odor_pos = nanmean(odor(selec2&h_odor2&pos));
sem_odor_pos = std(odor(selec2&h_odor2&pos));%./sqrt(length(odor(selec2&h_odor2&pos)));
sem_odor_neg = std(odor(selec2&h_odor2&~pos));%./sqrt(length(odor(selec2&h_odor2&~pos)));

mean_odor = nanmean(odor(selec2&h_odor2));
sem_odor = std(odor(selec2&h_odor2));%./sqrt(length(odor(selec2&h_odor2)));

plot([mean_odor_neg+sem_odor_neg mean_odor_neg-sem_odor_neg],...
    [mean_odorLight_neg mean_odorLight_neg],'-','color',[1 1 1])
plot([mean_odor_neg mean_odor_neg],...
    [mean_odorLight_neg+sem_odorLight_neg mean_odorLight_neg-sem_odorLight_neg],'-','color',[1 1 1])
plot([mean_odor_pos+sem_odor_pos mean_odor_pos-sem_odor_pos],...
    [mean_odorLight_pos mean_odorLight_pos],'-','color',[1 1 1])
plot([mean_odor_pos mean_odor_pos],...
    [mean_odorLight_pos+sem_odorLight_pos mean_odorLight_pos-sem_odorLight_pos],'-','color',[1 1 1])

[~,p]=ttest(odor(selec2&h_odor2),odor_CL(selec2&h_odor2),'alpha',0.05);

title({ExpType,...
    ['slope: ' num2str(table2array(MCmdl1.Coefficients(2,1)),3) ' adj r2 = ' num2str(MCmdl1.Rsquared.Adjusted,3),'; n = ',num2str(sum(selec2&h_odor2&~pos)),' odor-inhibited cells'],...
    ['slope: ' num2str(table2array(MCmdl2.Coefficients(2,1)),3) ' adj r2 = ' num2str(MCmdl2.Rsquared.Adjusted,3),'; n = ',num2str(sum(selec2&h_odor2&pos)),' odor-excited cells'],...
    ['ttest, p = ',num2str(p,3),'; n = ',num2str(sum(selec2&h_odor2)),' cell-odor pair']},...
    'Interpreter','none')
xlim([-20 30])
ylim([-40 30])

if saveFig
    %     saveas(gcf,[saveName_full,'OdorVsOdor+Light'],'epsc');
    %     saveas(gcf,[saveName_full,'OdorVsOdor+Light.m']);
    screen2png([saveName_full,'OdorVsOdor+Light']);
    set(gcf,'Renderer','Painter')
    hgexport(gcf,[saveName_full,'OdorVsOdor+Light.eps']);
end

end
