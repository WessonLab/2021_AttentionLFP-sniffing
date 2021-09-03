%Make plots for PAC figure

%% setup stuff from computePAC script
clear;
addpath(genpath('C:\TDT\TDTMatlabSDK\TDTSDK'));
addpath(genpath('C:\Users\smelluser\Documents\MATLAB\chronux_2_12'));
addpath(genpath('H:\Dropbox (Wesson Lab)\Wesson Lab general\code\circstat-matlab-master'));

%user input
spreadsheet = 'allfilenames_allblocks.xlsx';
[~,files, ~] = xlsread(spreadsheet, 'A2:A24');
loadfiletag = '_6andswitchblocks.mat';
baseline = [-2 -1.7]; %relative to odor onset
specgramwindow = [0.6 0.01]; %first value is the window size, second value is window step
tasktypes = {'baselinepower','odoronly','toneattn','odorattn','switch'};
outcomes = {'correct','incorrect','somecorrect','all','prevcorrect','previncorrect','someprevcorrect','incorrect_all','somecorrect_all'};
structures = {'OBLFP','PFCLFP','OTLFP','OBPFC_cohgram','OBOT_cohgram','PFCOT_cohgram'};
structureslong = {'OBLFP','PFCLFP','OTLFP','OBLFP2','PFCLFP2','OTLFP2','Cpoke','Rpoke','Lpoke','odo_to_cout','cout_to_choice','block','trial','odor','tone','trialonset','blockpcor','OBLFP_theta'};%everything for combining correct/incorrect structures
cohgramstructures = {'noshuf','shuf1','shuf2','shufboth'};
agstructures = {'specgrams','cohgrams','coherence','power'};
trialepochs = {'cohgram','pre','hold','odor'};
trialtypes = {'congruent','incongruent','all'};
bands = {'all','theta','beta','lowgamma','highgamma'};
makeplots = 0; %set to zero if you don't want plots
computecomodulograms = 1; %set to zero if you don't want to do the whole comodulogram, just the specific phase/amp pair
thetarange = [2 12];
betarange = [15 35];
lowgammarange = [40 60];
highgammarange = [60 80];
epochindex = [0,0;2035,3052;3053,4071;4072,4477;];

%% PAC stuff
% Define the amplitude- and phase-frequencies for PAC
PhaseFreqVector=2:2:20;
AmpFreqVector=10:5:120;

PhaseFreq_BandWidth=4;
AmpFreq_BandWidth=20;

% Define phase bins for PAC
nbin = 51; % number of phase bins was 18 making it 51
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end

%for MI/MeanAmp calculation for specific phase/frequency pair
Pf1 = 2;
Pf2 = 10;
Af1 = 65;
Af2 = 100;

%% set up struct of all files divided by rat
allrats = {'RAT064','RAT066','RAT071','RAT073','RAT074'};
allfiles.RAT064 = files(1:4);
allfiles.RAT066 = files(5:9);
allfiles.RAT071 = files(10:15);
allfiles.RAT073 = files(16:20);
allfiles.RAT074 = files(21:23);

%% 

for r = 1:length(allrats)
    clear PeakPhaseAngle_allfiles
    for f=1:length(allfiles.(allrats{r})) %for each file per rat
        load(strcat(allfiles.(allrats{r}){f},loadfiletag)); %load the file
        %For this file, plot a trial by trial of the MeanAmp
        clear MeanAmp_alltrls Outcome_alltrls Rstat_alltrls AllStructTrls sorted lfp wsbtrials;

        %randomly delete half of the correct odor only trials, so there are a
        %comparable number to the number of TA and OA trials
        OOtrls = size(alldata.toneoff.odoronly.correct.all.OBLFP,1); %how many OO trls are there
        structures2del = {'OBLFP','PFCLFP','OTLFP','OBLFP2','PFCLFP2','OTLFP2','Cpoke','Rpoke','Lpoke',...
                           'odo_to_cout','cout_to_choice','block','trial','odor','tone','trialonset',...
                           'blockpcor','OBLFP_theta','OBmorlet','OBmorletpeaks','OBmorletpeaklocs','OBmorletfreqs',...
                           'OBmorletcycles','sniffs','sniffPSTH','SDF_allTrials','SDFsmooth_allTrials'}; % delete not only the neural data but all dat aassociated with the trials. Note that the sniff SDF means will still contain all trials
        for t=1:(round(OOtrls/2)) %for half the number of trials
            randtrl = randperm(size(alldata.toneoff.odoronly.correct.all.OBLFP,1)); %how many trials remain (will be same first round, and decrease by 1 with each row deleted
            for s=1:length(structures2del)
                alldata.toneoff.odoronly.correct.all.(structures2del{s})(randtrl(1),:)=[];
            end
        end
        
        % for switch data, combine corrects and incorrects
        for x=1:length(structureslong)%for switch data, combine corrects and incorrects
            alldata.toneoff.switch.all.all.(structureslong{x}) = [alldata.toneoff.switch.correct.all.(structureslong{x});alldata.toneoff.switch.incorrect.all.(structureslong{x})];
        end
        
        %add switch corrects and incorrects into the all structure
        for x=1:length(structureslong)
            alldata.toneoff.all.(structureslong{x}) = [alldata.toneoff.all.(structureslong{x});alldata.toneoff.switch.all.all.(structureslong{x})];
        end
        
        %reorder the all structure
        AllStructTrls(:,1) = alldata.toneoff.all.block;
        AllStructTrls(:,2) = alldata.toneoff.all.trial;
        AllStructTrls(:,3:6108) = alldata.toneoff.all.OBLFP;
        sorted = sortrows(AllStructTrls,[1,2]); 
        
        lfp = sorted(:,612:5190);
        data_length = length(lfp);
        srate = alldata.info.LFPfs;
        dt = 1/srate;
        t = (1:data_length)*dt;
        t=t-3; %making sure 0 = odor onset
        phaseangle = 10:7:360; %pretty sure this is right
        for ttt = 1:length(alldata.toneoff.all.trial)
            [MI,MeanAmp] = ModIndex_v1(lfp(ttt,:),srate,Pf1,Pf2,Af1,Af2,position);
            MI_alltrls(ttt) = MI;
            MeanAmp_alltrls(ttt,:) = smoothdata(MeanAmp/sum(MeanAmp));
            
            %find the outcome of the current trial and save it for later
            %plotting. 2=correct, 1=incorrect 0=omit
            currtrl = [sorted(ttt,1) sorted(ttt,2)];
            wsbdata = wholeshebang(2:length(wholeshebang),:); %delete header row of wsb
            wsbtrials(:,1) = [wsbdata{:,2}]'; %just block/trial info
            wsbtrials(:,2) = [wsbdata{:,3}]';
            wsbindex = intersect(find(wsbtrials(:,1)==currtrl(1)),find(wsbtrials(:,2)==currtrl(2))) +1;
            if strcmp(wholeshebang(wsbindex,7),'correct')
                Outcome_alltrls(ttt,1) = 1.5;
            elseif strcmp(wholeshebang(wsbindex,7),'incorrect')
                Outcome_alltrls(ttt,1) = 6;
            else
                Outcome_alltrls(ttt,1) = 0;
            end
            
            %put the trial type in the same outcome matrix
            %3=odoronly 4=toneattn 5=switch 6=odorattn
            if strcmp(wholeshebang{wsbindex,1},'odoronly')
                Outcome_alltrls(ttt,2) = 3.75;
            elseif strcmp(wholeshebang{wsbindex,1},'toneattn')
                Outcome_alltrls(ttt,2) = 2.75;
            elseif strcmp(wholeshebang{wsbindex,1},'switch')
                Outcome_alltrls(ttt,2) = 5.3;
            elseif strcmp(wholeshebang{wsbindex,1},'odorattn')
                Outcome_alltrls(ttt,2) = 4.8;
            end
        end
        
        %find peak phase angle for each trial type
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                outcomes = {'correct','incorrect','all'};
            else
                outcomes = {'correct','incorrect'};
            end
            for o=1:length(outcomes)
                clear PeakPhaseAngle_alltrls;
                if ~isempty(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.trial)
                    lfp = alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.OBLFP(:,610:5188);
                    for ttt = 1:size(lfp,1)
                        %Peak phase angle
                        [MI,MeanAmp] = ModIndex_v1(lfp(ttt,:),srate,Pf1,Pf2,Af1,Af2,position);
                        phaseangle = 10:7:360;
                        [peak,loc] = max(MeanAmp);%find highest value in MeanAmp
                        PeakPhaseAngle_alltrls(ttt) = phaseangle(loc);%use loc to find the angle for the peakvalue from the phaseangle variable and save it
                    end
                    PeakPhaseAngle.(tasktypes{a}).(outcomes{o}) = PeakPhaseAngle_alltrls;
                    Rstat.(tasktypes{a}).(outcomes{o}) = circ_rtest(PeakPhaseAngle_alltrls*(pi/180));
                
                    PeakPhaseAngle_allfiles.(tasktypes{a}).(outcomes{o})(f,1:100) = NaN;
                    PeakPhaseAngle_allfiles.(tasktypes{a}).(outcomes{o})(f,1:length(PeakPhaseAngle_alltrls)) = PeakPhaseAngle_alltrls;
                else
                    PeakPhaseAngle.(tasktypes{a}).(outcomes{o}) = NaN;
                    Rstat.(tasktypes{a}).(outcomes{o}) = NaN;
                    PeakPhaseAngle_allfiles.(tasktypes{a}).(outcomes{o})(f,1:100) = NaN;
                end
            end
        end
        
        figure('Name',strcat('Phase Amplitude relationship for: ',allrats{r},' file ',num2str(f)));
        set(gcf,'renderer','Painters');
        subplot(2,2,1);
        imagesc(Outcome_alltrls);
        ylabel('Trial #');
        xlabel ('  Outcome               Task     ');
        subplot(2,2,2);
        imagesc(phaseangle,1:length(MeanAmp_alltrls),MeanAmp_alltrls);
        colormap jet;
        colorbar;
        caxis([.015 .035]);
        xlabel('Phase angle');
        
        subplot(2,2,4)
        colors_tasktypes = [0.4660, 0.6740, 0.1880; 0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;1 .2 .8]; %green: 0.4660, 0.6740, 0.1880; blue: 0, 0.4470, 0.7410; orange: 0.8500, 0.3250, 0.0980; pink: 1 .2 .8 - match these up with the tasktypes variable.
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                PeakAngle = PeakPhaseAngle.(tasktypes{a}).all;
            else
                PeakAngle = PeakPhaseAngle.(tasktypes{a}).correct;
            end
            if ~isnan(PeakAngle)
                PeakAngle = PeakAngle*(pi/180);
                polarhistogram(PeakAngle,24,'FaceColor',colors_tasktypes(a-1,:),'EdgeColor',colors_tasktypes(a-1,:),'FaceAlpha',.3);
            end
            hold on;
        end
        legend(strcat(tasktypes{2},', p=',num2str(Rstat.(tasktypes{2}).correct)),strcat(tasktypes{3},', p=',num2str(Rstat.(tasktypes{3}).correct)),strcat(tasktypes{4},', p=',num2str(Rstat.(tasktypes{4}).correct)), strcat(tasktypes{5},', p=',num2str(Rstat.(tasktypes{5}).all)),'Location','southoutside');
        
        subplot(2,2,3)
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                PeakAngle = [PeakPhaseAngle.(tasktypes{a}).all,PeakPhaseAngle.(tasktypes{a}).all+360];
            else
                PeakAngle = [PeakPhaseAngle.(tasktypes{a}).correct,PeakPhaseAngle.(tasktypes{a}).correct+360];
            end
            if ~isnan(PeakAngle)
                histogram(PeakAngle,48,'FaceColor',colors_tasktypes(a-1,:),'EdgeColor',colors_tasktypes(a-1,:),'FaceAlpha',.3);
            end
            hold on;
        end
        %set x axis limits to peak phase angle +-60 degrees
        midpoint = median(PeakPhaseAngle.odorattn.correct(PeakPhaseAngle.odorattn.correct<360));
        xlim([midpoint-100,midpoint+100]);
        
        %ks test
        %OO vs. TA, OO vs. SW, OO vs. OA
        x1=PeakPhaseAngle.odoronly.correct;
        x2=PeakPhaseAngle.toneattn.correct;
        
        [h,p]=kstest2(x1,x2);
        text(mean(PeakPhaseAngle.odorattn.correct)-100,14,strcat('OO vs. TA ks: h=',num2str(h),' p=',num2str(p)));
        
        x2=PeakPhaseAngle.switch.all;
        if ~isnan(x2)
        [h,p]=kstest2(x1,x2);
        text(mean(PeakPhaseAngle.odorattn.correct)-100,13,strcat('OO vs. SW ks: h=',num2str(h),' p=',num2str(p)));
        end
        
        x2=PeakPhaseAngle.odorattn.correct;
        [h,p]=kstest2(x1,x2);
        text(mean(PeakPhaseAngle.odorattn.correct)-100,12,strcat('OO vs. OA ks: h=',num2str(h),' p=',num2str(p)));
        
        %TA vs. SW, TA vs. OA
        x1=PeakPhaseAngle.toneattn.correct;
        x2=PeakPhaseAngle.switch.all;
        
        if ~isnan(x2)
        [h,p]=kstest2(x1,x2);
        text(mean(PeakPhaseAngle.odorattn.correct)-100,11,strcat('TA vs. SW ks: h=',num2str(h),' p=',num2str(p)));
        end
        
        x2=PeakPhaseAngle.odorattn.correct;
        [h,p]=kstest2(x1,x2);
        text(mean(PeakPhaseAngle.odorattn.correct)-100,10,strcat('TA vs. OA ks: h=',num2str(h),' p=',num2str(p)));
        
        %SW vs. OA
        x1=PeakPhaseAngle.switch.all;
        x2=PeakPhaseAngle.odorattn.correct;
        
        if ~isnan(x1)
        [h,p]=kstest2(x1,x2);
        text(mean(PeakPhaseAngle.odorattn.correct)-100,9,strcat('SW vs. OA ks: h=',num2str(h),' p=',num2str(p)));
        end
        
    end
    
    %compute Rstat and plot polar histogram for all incorrects across files
    %for this rat
    figure('Name',strcat(allrats{r},' all files, correct vs. incorrect'));
    set(gcf,'renderer','Painters');
    for a=2:length(tasktypes)
        corrects = PeakPhaseAngle_allfiles.(tasktypes{a}).correct;
        incorrects = PeakPhaseAngle_allfiles.(tasktypes{a}).incorrect;
        PeakAngle_correct = (corrects(~isnan(corrects)));
        PeakAngle_incorrect = (incorrects(~isnan(incorrects)));
        
        %randomly delete corrects until it's same size as incorrects
        for ttt=1:(length(PeakAngle_correct)-length(PeakAngle_incorrect)) %for half the number of trials
            randtrl = randperm(length(PeakAngle_correct)); %how many trials remain (will be same first round, and decrease by 1 with each row deleted
            PeakAngle_correct(randtrl(1)) = [];
        end
        
        subplot(2,2,a-1);
        PeakAngle_correct = PeakAngle_correct*(pi/180);
        polarhistogram(PeakAngle_correct,24,'FaceColor',colors_tasktypes(a-1,:),'EdgeColor',colors_tasktypes(a-1,:),'FaceAlpha',.3);
        hold on;
        PeakAngle_incorrect = PeakAngle_incorrect*(pi/180);
        polarhistogram(PeakAngle_incorrect,24,'FaceColor',colors_tasktypes(a-1,:),'EdgeColor','k','FaceAlpha',.3);
        Rstat_correct = circ_rtest(PeakAngle_correct*(pi/180));
        Rstat_incorrect = circ_rtest(PeakAngle_incorrect*(pi/180));
        %legend(strcat(tasktypes{a},' correct, p=',num2str(Rstat_correct)),strcat(tasktypes{a},' incorrect, p=',num2str(Rstat_incorrect)));
        
        %ks test for correct/incorrect trials across all sessions
        x1=PeakPhaseAngle_allfiles.(tasktypes{a}).correct;
        x1=(x1(~isnan(x1)));
        x2=PeakPhaseAngle_allfiles.(tasktypes{a}).incorrect;
        x2=(x2(~isnan(x2)));
        [h,p]=kstest2(x1,x2);
        %plot a dummy so that results of ks test can go in legend
        dummy = 0; %polar histogram wont plot a dummy of Nan, so doing 0
        polarhistogram(dummy,24);
        legend(strcat(tasktypes{a},' correct, p=',num2str(Rstat_correct)),strcat(tasktypes{a},' incorrect, p=',num2str(Rstat_incorrect)),strcat((tasktypes{a}),', correct vs. incorrect ks: h=',num2str(h),' p=',num2str(p)));
        %legend(strcat((tasktypes{a}),', correct vs. incorrect ks: h=',num2str(h),' p=',num2str(p)));
        
    end

    figure('Name',strcat('All trials across sessions for:  ',allrats{r}));
    set(gcf,'renderer','Painters');
    subplot(1,3,1);
     for a=2:length(tasktypes)
        if strcmp(tasktypes{a},'switch')
            alltrls = PeakPhaseAngle_allfiles.(tasktypes{a}).all;
            PeakAngle = (alltrls(~isnan(alltrls)));
            Rstat_allfiles.(tasktypes{a}).all = circ_rtest(PeakAngle*(pi/180));
        else
            alltrls = PeakPhaseAngle_allfiles.(tasktypes{a}).correct;
            PeakAngle = (alltrls(~isnan(alltrls)));
            Rstat_allfiles.(tasktypes{a}).correct = circ_rtest(PeakAngle*(pi/180));
        end
        if ~isnan(PeakAngle)
            PeakAngle = PeakAngle*(pi/180);
            polarhistogram(PeakAngle,24,'FaceColor',colors_tasktypes(a-1,:),'EdgeColor',colors_tasktypes(a-1,:),'FaceAlpha',.3);

        end
        hold on;
    end
        legend(strcat(tasktypes{2},', p=',num2str(Rstat_allfiles.(tasktypes{2}).correct)),strcat(tasktypes{3},', p=',num2str(Rstat.(tasktypes{3}).correct)),strcat(tasktypes{4},', p=',num2str(Rstat.(tasktypes{4}).correct)), strcat(tasktypes{5},', p=',num2str(Rstat.(tasktypes{5}).all)),'Location','southoutside');   

    subplot(1,3,2);
    for a=2:length(tasktypes)
        if strcmp(tasktypes{a},'switch')
            alltrls = PeakPhaseAngle_allfiles.(tasktypes{a}).all;
            PeakAngle = (alltrls(~isnan(alltrls)));
            Rstat_allfiles.(tasktypes{a}).all = circ_rtest(PeakAngle*(pi/180));
        else
            alltrls = PeakPhaseAngle_allfiles.(tasktypes{a}).correct;
            PeakAngle = (alltrls(~isnan(alltrls)));
            Rstat_allfiles.(tasktypes{a}).correct = circ_rtest(PeakAngle*(pi/180)); %calculate rstat based on peak angle without duplicating it + 360
        end
        if ~isnan(PeakAngle)
            PeakAngle = [PeakAngle;PeakAngle+360]; %for plotting purposes
            histogram(PeakAngle,48,'FaceColor',colors_tasktypes(a-1,:),'EdgeColor',colors_tasktypes(a-1,:),'FaceAlpha',.3);
        end
        hold on;
    end
    %set x axis limits to peak phase angle +-60 degrees
    
    alltrls = [PeakPhaseAngle_allfiles.odorattn.correct,PeakPhaseAngle_allfiles.odorattn.correct+360];
    PeakAngle = (alltrls(~isnan(alltrls)));
    midpoint = median(PeakAngle(PeakAngle<360));
    xlim([midpoint-100,midpoint+100]);
    subplot(1,3,3)
    xlim([1 2]);
    %ks test
    %OO vs. TA, OO vs. SW, OO vs. OA
    x1=PeakPhaseAngle_allfiles.odoronly.correct;
    x1=(x1(~isnan(x1)));
    x2=PeakPhaseAngle_allfiles.toneattn.correct;
    x2=(x2(~isnan(x2)));

    [h,p]=kstest2(x1,x2);
    text(1,1,strcat('OO vs. TA ks: h=',num2str(h),' p=',num2str(p)));

    x2=PeakPhaseAngle_allfiles.switch.all;
    x2=(x2(~isnan(x2)));
    if ~isnan(x2)
    [h,p]=kstest2(x1,x2);
    text(1,0.9,strcat('OO vs. SW ks: h=',num2str(h),' p=',num2str(p)));
    end

    x2=PeakPhaseAngle_allfiles.odorattn.correct;
    x2=(x2(~isnan(x2)));
    [h,p]=kstest2(x1,x2);
    text(1,0.8,strcat('OO vs. OA ks: h=',num2str(h),' p=',num2str(p)));

    %TA vs. SW, TA vs. OA
    x1=PeakPhaseAngle_allfiles.toneattn.correct;
    x1=(x1(~isnan(x1)));
    x2=PeakPhaseAngle_allfiles.switch.all;
    x2=(x2(~isnan(x2)));

    if ~isnan(x2)
    [h,p]=kstest2(x1,x2);
    text(1,0.7,strcat('TA vs. SW ks: h=',num2str(h),' p=',num2str(p)));
    end

    x2=PeakPhaseAngle_allfiles.odorattn.correct;
    x2=(x2(~isnan(x2)));
    [h,p]=kstest2(x1,x2);
    text(1,0.6,strcat('TA vs. OA ks: h=',num2str(h),' p=',num2str(p)));

    %SW vs. OA
    x1=PeakPhaseAngle_allfiles.switch.all;
    x1=(x1(~isnan(x1)));
    x2=PeakPhaseAngle_allfiles.odorattn.correct;
    x2=(x2(~isnan(x2)));

    if ~isnan(x1)
    [h,p]=kstest2(x1,x2);
    text(1,0.5,strcat('SW vs. OA ks: h=',num2str(h),' p=',num2str(p)));
    end
end

%save all figures
FigList = findobj(allchild(0),'flat','Type','figure');
savefig(FigList,'ThetaGammaPACplots_withkstest');
