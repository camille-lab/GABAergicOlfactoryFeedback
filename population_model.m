%% build a population rate model of mitral cell-granule cell interactions
% Please refere to Mazo et al., Nat Comm 2022
% srcipt used to build and plot the model presented in Fig. 5

% population dynamics:
% dMC = -MC + I - x*GC - y*FB (I is FF odor input)
% dGC = -GC + w*MC - z*FB
% We don't need to simulate this because we're just interested
% in the steady states, i.e. where the firing rates of both
% populations go as they stabilize over time.
% If we set dMC and dGC to zero (i.e., both populations don't
% change their rates, we can find the solution to this, and then
% we just need to rearrange terms to get equations of two lines,
% y = ax + b, with y being the GC population rate, and x being
% the MC population rate. The equations when dMC and dGC are set
% to zero are called nullclines.

% null-clines (steady-states of population equations)
% dMC = -MC + I - x*GC - y*FB
% --> MC = (1/x)*(I-MC-y*FB)
% dGC = -GC + w*MC - z*FB
% --> GC = w*MC - z*FB

%% parameters
v = 1;    % GC to dSAC strength
w = 1;    % MC to GC strength
x = .5;   % GC to MC strength
y = 2.0;  % FB to MC strength
z = 4.0;  % FB to GC strength

% FF odor input can be fixed, results won't depend on it.
I = 15;

% For plotting, we can consider MC population rates of 10-20 Hz
MCvals = linspace(10,20,101);

w_vector = linspace(0,2,41); w_vector = w_vector(2:end);
FBvector = linspace(0,2,41); FBvector = FBvector(2:end);
for i = 1:length(w_vector)
    for ii = 1:length(FBvector)
        w=w_vector(i); x=w_vector(i);
        y=FBvector(ii); z=2*y;
        
        % case 1: no feedback stimulation (just MC, GC dynamics)
        FB = 0; SAC = 1;
        MC_NC = (1/x)*(I - MCvals - y*FB);
        GC_NC = w*MCvals - z*FB - v*SAC;
        fixpt_MC = 1/(w*x+1)*(I+(z*x-y)*FB + v*SAC*x);
        fixpt_GC = w*fixpt_MC - z*FB - v*SAC;
        
        % case 2: perturbation of inhibitory FB, so FB=1
        FB = 1; SAC = .1/z;
        MC_NC2 = (1./x)*(I - MCvals - y*FB);
        GC_NC2 = w*MCvals - z*FB - v*SAC;
        fixpt_MC2 = 1.0/(w*x+1)*(I+(z*x-y)*FB + v*SAC*x);
        fixpt_GC2 = w*fixpt_MC2 - z*FB - v*SAC;
                
        dMC(i,ii) = fixpt_MC2-fixpt_MC;
        dGC(i,ii) = fixpt_GC2-fixpt_GC;
        if w_vector(i) == .3 && (FBvector(ii) == .1 || FBvector(ii) == .9 || FBvector(ii) == 1.8)
            if FBvector(ii) == .1
                figure;
            subplot(1,2,1);
            hold on;
            plot(MCvals,MC_NC,'Color','red');
            plot(MCvals,GC_NC,'Color','blue');
            end
            subplot(1,2,1);
            plot(MCvals,MC_NC2,'LineStyle','--','Color','red');
            plot(MCvals,GC_NC2,'LineStyle','--','Color','blue');
            plot(fixpt_MC,fixpt_GC,'.','MarkerSize',20,'Color','black');
            plot(fixpt_MC2,fixpt_GC2,'.','MarkerSize',20,'Color','magenta');
            xlabel('Mitral cell population rate (Hz)');
            ylabel('Granule cell population rate (Hz)');
        end
    end
    MCfr(i) = fixpt_MC;
    GCfr(i) = fixpt_GC;
end

wgc_value = find(w_vector == .3);
fb_values(1) = find(w_vector == .1); fb_values(2) = find(w_vector == .9);
fb_values(3) = find(w_vector == 1.8); %fb_values(4) = find(w_vector == .8);
subplot(1,2,2);
bar([1:2],[dMC(wgc_value,fb_values);dGC(wgc_value,fb_values)]);
xticklabels({'MC','GC'});
ylabel('Change in rate with FB');

%% heatmap
maxMC = max(dMC,[],'all'); minMC = min(dMC,[],'all'); climMC = max(abs(minMC),maxMC);
maxGC = max(dGC,[],'all'); minGC = min(dGC,[],'all'); climGC = max(abs(minGC),maxGC);
figure;
subplot(2,2,1); hold on
imagesc(dMC,[-climMC climMC]);
xlabel('FB strength'); ylabel('GC-MC interaction strength')
xticks([1:length(FBvector)]);xticklabels(FBvector);xlim([1 length(FBvector)]);
yticks([1:length(w_vector)]);yticklabels(w_vector);ylim([1 length(FBvector)]);
scatter(fb_values,[wgc_value wgc_value wgc_value],'k')
title('MC'); axis square; set(gca, 'YDir','normal'); colorbar ;

subplot(2,2,2); hold on
imagesc(dGC,[-climGC climGC]); title('GC'); axis square; set(gca, 'YDir','normal'); colorbar
xlabel('FB strength')
xticks([1:length(FBvector)]);xticklabels(FBvector);xlim([1 length(FBvector)]);
yticks([1:length(w_vector)]);yticklabels(w_vector);ylim([1 length(FBvector)]);
scatter(fb_values,[wgc_value wgc_value wgc_value],'k')

colormap diverging_bwr_40_95_c42_n256

subplot(2,2,3); hold on;
p1 = plot(w_vector,MCfr,'r'); p2 = plot(w_vector,GCfr,'b');
plot([.3 .3],[0 15],'k--')

axis square;
text(15,2,['red:MC' newline 'blue: GC'],...
    'horizontalalignment','left','verticalalignment','top')
xlabel('basal FR'); ylabel('GC-MC interaction strength')
