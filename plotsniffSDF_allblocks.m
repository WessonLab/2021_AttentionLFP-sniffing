clear;
load('RAT066-032419-1allblocks_finalized2_sniffs.mat');
plot_psthwin = [-2000 2000];
SDFxvals = (-2000:50:1999);
variance = 1; %1 for SEM, 2 for 2stds, 3 for CI
varianceoptions = {'_meanpmSEM','_meanpm2std','_meanpm95CI'};          
OOcolor = [0.4660,0.6740,0.1880];
TAcolor = [0,0.4470,0.7410];
OAcolor = [0.8500,0.3250,0.0980];

switchblocks = [11 12 13];
OAblocks = [14 15 16 17 18 19 20];

switchindex_CONG =   [min(find(alldata.toneoff.odorattn.correct.congruent.block==min(switchblocks)))...
                      max(find(alldata.toneoff.odorattn.correct.congruent.block==max(switchblocks)))];
switchindex_INCONG = [min(find(alldata.toneoff.odorattn.correct.incongruent.block==min(switchblocks)))...
                      max(find(alldata.toneoff.odorattn.correct.incongruent.block==max(switchblocks)))];
OAindex_CONG = [min(find(alldata.toneoff.odorattn.correct.congruent.block==min(OAblocks)))...
                max(find(alldata.toneoff.odorattn.correct.congruent.block==max(OAblocks)))];
OAindex_INCONG = [min(find(alldata.toneoff.odorattn.correct.incongruent.block==min(OAblocks)))...
                  max(find(alldata.toneoff.odorattn.correct.incongruent.block==max(OAblocks)))];

switchCONGsniffs = mean(alldata.toneoff.odorattn.correct.congruent.SDFsmooth_allTrials(switchindex_CONG(1):switchindex_CONG(2),:));
switchCONGsniffs_SEM = std(alldata.toneoff.odorattn.correct.congruent.SDFsmooth_allTrials(switchindex_CONG(1):switchindex_CONG(2),:))/sqrt(size(alldata.toneoff.odorattn.correct.congruent.SDFsmooth_allTrials(switchindex_CONG(1):switchindex_CONG(2),:),1));
switchCONGsniffs_error(1,:) = switchCONGsniffs - switchCONGsniffs_SEM;
switchCONGsniffs_error(2,:) = switchCONGsniffs + switchCONGsniffs_SEM;

switchINCONGsniffs = mean(alldata.toneoff.odorattn.correct.incongruent.SDFsmooth_allTrials(switchindex_INCONG(1):switchindex_INCONG(2),:));
switchINCONGsniffs_SEM = std(alldata.toneoff.odorattn.correct.incongruent.SDFsmooth_allTrials(switchindex_INCONG(1):switchindex_INCONG(2),:))/sqrt(size(alldata.toneoff.odorattn.correct.incongruent.SDFsmooth_allTrials(switchindex_INCONG(1):switchindex_INCONG(2),:),1));
switchINCONGsniffs_error(1,:) = switchINCONGsniffs - switchINCONGsniffs_SEM;
switchINCONGsniffs_error(2,:) = switchINCONGsniffs + switchINCONGsniffs_SEM;

OACONGsniffs = mean(alldata.toneoff.odorattn.correct.congruent.SDFsmooth_allTrials(OAindex_CONG(1):OAindex_CONG(2),:));
OACONGsniffs_SEM = std(alldata.toneoff.odorattn.correct.congruent.SDFsmooth_allTrials(OAindex_CONG(1):OAindex_CONG(2),:))/sqrt(size(alldata.toneoff.odorattn.correct.congruent.SDFsmooth_allTrials(OAindex_CONG(1):OAindex_CONG(2),:),1));
OACONGsniffs_error(1,:) = OACONGsniffs - OACONGsniffs_SEM;
OACONGsniffs_error(2,:) = OACONGsniffs + OACONGsniffs_SEM;

OAINCONGsniffs = mean(alldata.toneoff.odorattn.correct.incongruent.SDFsmooth_allTrials(OAindex_INCONG(1):OAindex_INCONG(2),:));
OAINCONGsniffs_SEM = std(alldata.toneoff.odorattn.correct.incongruent.SDFsmooth_allTrials(OAindex_INCONG(1):OAindex_INCONG(2),:))/sqrt(size(alldata.toneoff.odorattn.correct.incongruent.SDFsmooth_allTrials(OAindex_INCONG(1):OAindex_INCONG(2),:),1));
OAINCONGsniffs_error(1,:) = OAINCONGsniffs - OAINCONGsniffs_SEM;
OAINCONGsniffs_error(2,:) = OAINCONGsniffs + OAINCONGsniffs_SEM;
          
figure; %this are plotted in a weird order so the legend comes out how I want...
set(gcf,'renderer','Painters'); %changes the matlab renderer so that when I save files as .eps, they will actually be in vector format
%I think the patch thing makes them use a different renderer? Without
%this line, saving as an eps results in a tiled image rather than a vector
title(strcat('RAT066 crit/noncrit congruent/incongruent'));

plot(SDFxvals,switchCONGsniffs,'Color','k','LineWidth',2);
hold on;
plot(SDFxvals,switchINCONGsniffs,'Color','k','LineStyle','--','LineWidth',2);
plot(SDFxvals,OACONGsniffs,'Color',OAcolor,'LineWidth',2);
plot(SDFxvals,OAINCONGsniffs,'Color',OAcolor,'LineStyle','--','LineWidth',2);

%plot the error first
patch([SDFxvals fliplr(SDFxvals)],[switchCONGsniffs_error(2,:) fliplr(switchCONGsniffs_error(1,:))],'k','EdgeColor','none','FaceAlpha',0.2);
patch([SDFxvals fliplr(SDFxvals)],[switchINCONGsniffs_error(2,:) fliplr(switchINCONGsniffs_error(1,:))],'k','EdgeColor','none','FaceAlpha',0.2);
patch([SDFxvals fliplr(SDFxvals)],[OACONGsniffs_error(2,:) fliplr(OACONGsniffs_error(1,:))],OAcolor,'EdgeColor','none','FaceAlpha',0.2);
patch([SDFxvals fliplr(SDFxvals)],[OAINCONGsniffs_error(2,:) fliplr(OAINCONGsniffs_error(1,:))],OAcolor,'EdgeColor','none','FaceAlpha',0.2);


xlabel('Time from odor onset');
ylabel('Frequency(Hz)');
legend('congruent, <80%','incongruent, <80%','congruent, >80%','incongruent, >80%');
ylim([0 10]);


