%% Script to categorize the cell responses to optogentic stimulation of GABAergic feedback on spontaneous activity
% Please refere to Mazo et al., Nat Comm 2022
% srcipt used to analyze and plot data presented in Supplementary Fig. 8b (MC and TC)
% written by Camille Mazo

% Inputs: 
% data is the data to analyze: MC or TC. For light control, one has to
%   concatenate the data first. For instance for MC CL ctrl:
%   Light_control = cat(3,MClightControl{2,:});
% tAna and tBase are the frames to calculate the response and baseline
% alpha is the criterion for light-responsive

% for the analyzes in the paper, we used tAna=[135:150], tBase = [90:105],
% alpha = 0.01
% a typical call would be: 
% stats = ResponsiveCells(data.CL.light,[135:150],[90:105],0.01,'saveName')
%%

function stats = ResponsiveCells(data,tAna,tBase,alpha,saveName)
%% paired ttest to determine whether a cell is responsive or not
    
nCells = size(data,3);
resp = squeeze(nanmean(nanmean(data(:,tAna,:),1),2));
base = squeeze(nanmean(nanmean(data(:,tBase,:),1),2));
diff = resp - base;

neg = false(1,nCells); h_all = false(1,nCells);  p_all = h_all;
for i = 1:nCells
    [h,p] = ttest(nanmean(data(:,tAna,i),2),...
        nanmean(data(:,tBase,i),2),'alpha',alpha);
    h_all(i) = h; p_all(i) = p;
    if diff(i)<0
        neg = true;
    end
end

disp(['non-responsive cells: ', num2str(sum(~h_all)), ' = ',num2str(sum(~h_all)/nCells)])
disp(['excited: ', num2str(sum(h_all&~neg)), ' = ',num2str(sum(h_all&~neg)/nCells)])
disp(['inhibited: ',num2str(sum(h_all&neg)), ' = ', num2str(sum(h_all&neg)/nCells)])

stats.Info.tBase = tBase;
stats.Info.tAna = tAna;
stats.Info.alpha = alpha;
stats.Info.h_all = h_all;
stats.Info.p_all = p_all;
stats.nCells = nCells;
stats.nonResponsive = sum(~h_all);
stats.excited = sum(h_all&~neg);
stats.inhibited = sum(h_all&neg);

% save(['YourFolder' filesep,saveName],'stats')
end