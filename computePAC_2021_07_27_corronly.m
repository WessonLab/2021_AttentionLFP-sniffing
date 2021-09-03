%% calculates and saves (1) spectrograms (2) spectral power based on spectrograms
%(3) coherograms (4) coherence (5) spectral power based on mtspectrumc function 
%7.27.21 modifying this to only work for correct trials for OO TA OA and
%all trials for SW so I can pull things out faster. also only do the OB, no
%PFC or OT
clear;
addpath(genpath('C:\TDT\TDTMatlabSDK\TDTSDK'));
addpath(genpath('C:\Users\smelluser\Documents\MATLAB\chronux_2_12'));

%user input
spreadsheet = 'allfilenames_allblocks.xlsx';
[~,files, ~] = xlsread(spreadsheet, 'A2:A24');
loadfiletag = '_6andswitchblocks.mat';
%baseline = [-2 -1.7]; %relative to odor onset
%specgramwindow = [0.6 0.01]; %first value is the window size, second value is window step
tasktypes = {'baselinepower','odoronly','toneattn','odorattn','switch'};
outcomes = {'correct','all'}; %outcomes = {'correct','incorrect','somecorrect','all','prevcorrect','previncorrect','someprevcorrect','incorrect_all','somecorrect_all'};
structures = {'OBLFP'}; %structures = {'OBLFP','PFCLFP','OTLFP','OBPFC_cohgram','OBOT_cohgram','PFCOT_cohgram'};
structureslong = {'OBLFP','PFCLFP','OTLFP','OBLFP2','PFCLFP2','OTLFP2','Cpoke','Rpoke','Lpoke','odo_to_cout','cout_to_choice','block','trial','odor','tone','trialonset','blockpcor','OBLFP_theta'};%everything for combining correct/incorrect structures
%cohgramstructures = {'noshuf','shuf1','shuf2','shufboth'};
%agstructures = {'specgrams','cohgrams','coherence','power'};
trialepochs = {'cohgram','pre','hold','odor'};
trialtypes = {'congruent','incongruent','all'};
%bands = {'all','theta','beta','lowgamma','highgamma'};
makeplots = 1; %set to zero if you don't want plots
computecomodulograms = 1; %set to zero if you don't want to do the whole comodulogram, just the specific phase/amp pair
%thetarange = [2 12];
%betarange = [15 35];
%lowgammarange = [40 60];
%highgammarange = [60 80];
%epochindex = [0,0;2035,3052;3053,4071;4072,4477;];

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
Pf1_MI = 2;
Pf2_MI = 10;
Af1_MI_beta = 15;
Af2_MI_beta = 35;
Af1_MI_gamma = 65;
Af2_MI_gamma = 100;

%% set up struct of all files divided by rat
allrats = {'RAT064','RAT066','RAT071','RAT073','RAT074'};
allfiles.RAT064 = files(1:4);
allfiles.RAT066 = files(5:9);
allfiles.RAT071 = files(10:15);
allfiles.RAT073 = files(16:20);
allfiles.RAT074 = files(21:23);


for r = 1:length(allrats) %for each rat
    %in between rats, clear variables because some rats will have more files than others
    clear Comodulograms shuffled subtracted incorrects;
    tasktypes = {'baselinepower','odoronly','toneattn','odorattn','switch'}; %reset order of this, which was changed to plot individual files
    for f=1:length(allfiles.(allrats{r})) %for each file per rat
        load(strcat(allfiles.(allrats{r}){f},loadfiletag)); %load the file
        params=struct('tapers',[3 5],'pad',0,'Fs',alldata.info.LFPfs,'fpass',[0.5 100],'trialave',1);
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
        for s = 1:3
            for x=1:length(structureslong)%for switch data, combine corrects and incorrects
                alldata.toneoff.switch.all.(structureslong{x}) = [alldata.toneoff.switch.correct.all.(structureslong{x});alldata.toneoff.switch.incorrect.all.(structureslong{x})];
            end
        end
        
        %PAC
        for a=2:length(tasktypes)
            clear Comodulogram_alltrls MI_alltrls MeanAmp_alltrls PeakPhaseAngle_alltrls
            if strcmp(tasktypes{a},'switch') %if switch, set outcome to all, otherwise correct
                outcomes = {'all'};
            else
                outcomes = {'correct'};
            end
            for s=1 %just OB
                for o=1:length(outcomes)
                    clear Comodulogram_alltrls MI_alltrls MeanAmp_alltrls PeakPhaseAngle_alltrls
                    if a<5
                        nblocks = length(unique(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.block));
                    elseif a==5
                        nblocks = length(unique(alldata.toneoff.(tasktypes{a}).(outcomes{o}).block));
                    end    
                    if nblocks>1 
                        if a<5
                            lfp = alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s})(:,610:5188); %need 4572 samples to use 2 hz as lowest phase frequency. This goes from -3 to + 1.5 sec from odor
                        elseif a==5
                            lfp = alldata.toneoff.(tasktypes{a}).(outcomes{o}).(structures{s})(:,610:5188); %need 4572 samples to use 2 hz as lowest phase frequency. This goes from -3 to + 1.5 sec from odor
                        end
                        data_length = length(lfp);
                        srate = alldata.info.LFPfs;
                        dt = 1/srate;
                        t = (1:data_length)*dt;
                        t=t-3; %making sure 0 = odor onset
                        for ttt = 1:size(lfp,1) %for each trial, compute the PAC
                            if computecomodulograms >0
                            %'CPU filtering'
                            Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
                            AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
                            PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

                            for ii=1:length(AmpFreqVector)
                                Af1 = AmpFreqVector(ii);
                                Af2=Af1+AmpFreq_BandWidth;
                                AmpFreq=eegfilt(lfp(ttt,:),srate,Af1,Af2); % filtering
                                AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
                            end

                            for jj=1:length(PhaseFreqVector)
                                Pf1 = PhaseFreqVector(jj);
                                Pf2 = Pf1 + PhaseFreq_BandWidth;
                                PhaseFreq=eegfilt(lfp(ttt,:),srate,Pf1,Pf2); % filtering 
                                PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % getting the phase time series
                            end

                            %'Comodulation loop'
                            counter1=0;
                            for ii=1:length(PhaseFreqVector)
                            counter1=counter1+1;

                                Pf1 = PhaseFreqVector(ii);
                                Pf2 = Pf1+PhaseFreq_BandWidth;

                                counter2=0;
                                for jj=1:length(AmpFreqVector)
                                counter2=counter2+1;

                                    Af1 = AmpFreqVector(jj);
                                    Af2 = Af1+AmpFreq_BandWidth;
                                    [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
                                    Comodulogram(counter1,counter2)=MI;
                                end
                            end
                            %save comodulogram for the trial
                            Comodulogram_alltrls(:,:,ttt) = Comodulogram;
                            end
                            
                            %THETA-BETA
                            %Comodulogram done: now compute MI/MeanAmp for theta/beta specifically and save
                            [MI_beta,MeanAmp_beta] = ModIndex_v1(lfp(ttt,:),srate,Pf1_MI,Pf2_MI,Af1_MI_beta,Af2_MI_beta,position);
                            MI_beta_alltrls(ttt) = MI_beta;
                            MeanAmp_beta_alltrls(ttt,:) = MeanAmp_beta;

                            %Peak phase angle
                            phaseangle = 10:7:360; %pretty sure this is right
                            %find highest value in MeanAmp
                            [peak_beta,loc_beta] = max(MeanAmp_beta);
                            %use loc to find the angle for the peakvalue from the phaseangle variable and save it
                            PeakPhaseAngle_beta_alltrls(ttt) = phaseangle(loc_beta);
                            
                            %THETA-GAMMA
                            %Comodulogram done: now compute MI/MeanAmp for theta/gamma specifically and save
                            [MI_gamma,MeanAmp_gamma] = ModIndex_v1(lfp(ttt,:),srate,Pf1_MI,Pf2_MI,Af1_MI_gamma,Af2_MI_gamma,position);
                            MI_gamma_alltrls(ttt) = MI_gamma;
                            MeanAmp_gamma_alltrls(ttt,:) = MeanAmp_gamma;

                            %Peak phase angle
                            phaseangle = 10:7:360; %pretty sure this is right
                            %find highest value in MeanAmp
                            [peak_gamma,loc_gamma] = max(MeanAmp_gamma);
                            %use loc to find the angle for the peakvalue from the phaseangle variable and save it
                            PeakPhaseAngle_gamma_alltrls(ttt) = phaseangle(loc_gamma);
                        end

                        %compute and save mean comodulogram
                        if computecomodulograms > 0
                        Comodulogram_mean = mean(Comodulogram,3);
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram(:,:,f) = Comodulogram_mean;
                        end

                        %compute and save mean MI and MeanAmp for BETA
                        MI_beta_mean = mean(MI_beta_alltrls);
                        MeanAmp_beta_mean = mean(MeanAmp_beta_alltrls);
                        PPAvar_beta = var(PeakPhaseAngle_beta_alltrls);
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta(f) = MI_beta_mean;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta(f,:) = MeanAmp_beta_mean;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_beta{f} = {PeakPhaseAngle_beta_alltrls};
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_beta(f) = PPAvar_beta;
                        
                        %compute and save mean MI and MeanAmp for GAMMA
                        MI_gamma_mean = mean(MI_gamma_alltrls);
                        MeanAmp_gamma_mean = mean(MeanAmp_gamma_alltrls);
                        PPAvar_gamma = var(PeakPhaseAngle_gamma_alltrls);
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma(f) = MI_gamma_mean;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma(f,:) = MeanAmp_gamma_mean;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_gamma{f} = {PeakPhaseAngle_gamma_alltrls};
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_gamma(f) = PPAvar_gamma;

                    elseif nblocks==1 | nblocks==0 | nblocks==NaN
                        if computecomodulograms > 0
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram(:,:,f) = NaN;
                        end
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta(f) = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta(f,:) = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_beta{f} = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_beta(f) = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma(f) = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma(f,:) = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_gamma{f} = NaN;
                        Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_gamma(f) = NaN;
                    end
                end
            end
        end
    
        %Count up the number of trials contributing to each type
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                outcomes = {'all'};%if switch, cycle through outcomes. otherwise, just correct
                for o=1
                    Comodulograms.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = length(alldata.toneoff.(tasktypes{a}).(outcomes{o}).trial);
                end
            else
                outcomes = {'correct'};
                for o=1
                    Comodulograms.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = length(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.trial);
                end
            end
        end

    end %end each file per rat

    %-----------------------------------------------------------------------
    %ALL FILES FOR CURRENT RAT DONE
    %COMPUTE MEAN FOR THIS RAT AND SAVE IN Comodulograms_ALLRATS
    %-----------------------------------------------------------------------
    
    %plot each file
    for f=1:length(allfiles.(allrats{r}))
        %plot the mean comodulograms for each file
        figure('Name',strcat(allrats{r},' file# ',num2str(f)));
        colormap('jet');
        set(gcf,'renderer','Painters');
        sp=1;
        tasktypes = {'baselinepower','odoronly','toneattn','switch','odorattn'}; %change the order of task types so switch is plotted where I want
        spranges = [0,0.018;0,0.005;0,0.005];

        for a=2:length(tasktypes)%start with a=2 to skip baselinepower
            for s=1 %:3
                subplot(4,1,sp);
                if strcmp(tasktypes{a},'switch')
                    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms.(tasktypes{a}).all.all.(structures{s}).comodulogram(:,:,f)',30,'lines','none'); %switch, plot all trials
                else
                    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms.(tasktypes{a}).correct.all.(structures{s}).comodulogram(:,:,f)',30,'lines','none'); %other tasktypes, plot correct only
                end
                ylabel('Amplitude Frequency (Hz)')
                xlabel('Phase Frequency (Hz)')
                colorbar;
                title (strcat(tasktypes{a},structures{s}));
                caxis(spranges(s,:));
                xlim([4 12]);
                sp=sp+1;
            end
        end
    end
    
    %save data in Comodulograms_allsessions
    for a=2:length(tasktypes)
        for s=1
            if strcmp(tasktypes{a},'switch')
                outcomes = {'all'};
            else
                outcomes = {'correct'};
            end
            for o=1
                Comodulograms_allsessions.(allrats{r}).(tasktypes{a}).trlnums = Comodulograms.(tasktypes{a}).trlnums;
                Comodulograms_allsessions.(allrats{r}).(tasktypes{a}).(outcomes{o}).all.(structures{s}) = Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s});
            end
        end
    end
    
    for a=2:length(tasktypes)
        for s=1
            if strcmp(tasktypes{a},'switch')
                outcomes = {'all'};
                for o=1
                    if computecomodulograms > 0
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram(:,:,r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram,3);
                    end
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta(r,:) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_beta{r} = {Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_beta}; %save all peak angles don't average
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_beta(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_beta);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma(r,:) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_gamma{r} = {Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_gamma}; %save all peak angles don't average
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_gamma(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_gamma);
                end
            else
                outcomes = {'correct'};
                for o=1
                    if computecomodulograms > 0
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram(:,:,r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram,3);
                    end
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta(r,:) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_beta{r} = Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_beta; %save all peak angles don't average
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_beta(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_beta);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma(r,:) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma);
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_gamma{r} = Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngle_gamma; %save all peak angles don't average
                    Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_gamma(r) = nanmean(Comodulograms.(tasktypes{a}).(outcomes{o}).all.(structures{s}).PeakPhaseAngleVariance_gamma);
                end
            end
        end
        
        %average the trl nums
        if strcmp(tasktypes{a},'switch')
            outcomes = {'all'};
            for o=1
                for tt=3 %trial types - just all for now
                    Comodulograms_allrats.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt})(r) = nanmean(Comodulograms.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt}));
                end
            end
        else
            outcomes = {'correct'};
            for o=1
                for tt=3
                    Comodulograms_allrats.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt})(r) = nanmean(Comodulograms.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt}));
                end
            end
        end

    end

end %end all rats

Comodulograms_allrats.info.PAC = 'Computed using routines from Tort lab J Neurophys 2010 paper.';
if computecomodulograms > 0
    Comodulograms_allrats.info.AmpFreq = AmpFreq;
    Comodulograms_allrats.info.AmpFreqTransformed = AmpFreqTransformed;
    Comodulograms_allrats.info.PhaseFreq = PhaseFreq;
    Comodulograms_allrats.info.PhaseFreqTransformed = PhaseFreqTransformed;
end
Comodulograms_allrats.info.AmpFreq_BandWidth = AmpFreq_BandWidth;
Comodulograms_allrats.info.AmpFreqVector = AmpFreqVector;
Comodulograms_allrats.info.PhaseFreq_BandWidth = PhaseFreq_BandWidth;
Comodulograms_allrats.info.PhaseFreqVector = PhaseFreqVector;
Comodulograms_allrats.info.nbin = nbin;
Comodulograms_allrats.info.position = position;
Comodulograms_allrats.info.MI_Pf1 = Pf1_MI;
Comodulograms_allrats.info.MI_Pf2 = Pf2_MI;
Comodulograms_allrats.info.MI_Af1_beta = Af1_MI_beta;
Comodulograms_allrats.info.MI_Af2_beta = Af2_MI_beta;
Comodulograms_allrats.info.MI_Af1_gamma = Af1_MI_gamma;
Comodulograms_allrats.info.MI_Af2_gamma = Af2_MI_gamma;
Comodulograms_allrats.info.PhaseAngle = phaseangle;

%---------------------------------------------------------------------------
%AVERAGE ACROSS ALL RATS
%---------------------------------------------------------------------------
for a=2:length(tasktypes)
        
        for s=1 %average the comodulograms across rats
            if strcmp(tasktypes{a},'switch')
                outcomes = {'all'};
                for o=1
                    if computecomodulograms > 0
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram,3);
                    end
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta);
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta);
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma);
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma);
                end
            else
                outcomes = {'correct'};
                for o=1
                    if computecomodulograms > 0
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).comodulogram,3);
                    end
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_beta);
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_beta);
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MI_gamma);
                    Comodulograms_allrats_mean.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma = nanmean(Comodulograms_allrats.(tasktypes{a}).(outcomes{o}).all.(structures{s}).MeanAmp_gamma);
                end
            end
        end
end


save('Comodulograms_6andswitchblocks_2021_07_28b.mat','Comodulograms_allrats_mean','Comodulograms_allrats');    

%% 

%----------------
% MAKE PLOTS
%----------------

if makeplots > 0
   
    %plot the mean comodulograms ACROSS all rats
    figure('Name','Mean Comodulograms across rats');
    colormap('jet');
    set(gcf,'renderer','Painters');
    sp=1;
    tasktypes = {'baselinepower','odoronly','toneattn','switch','odorattn'}; %change the order of task types so switch is plotted where I want
    colors_tasktypes = [0.4660, 0.6740, 0.1880; 0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;1 .2 .8]; %green: 0.4660, 0.6740, 0.1880; blue: 0, 0.4470, 0.7410; orange: 0.8500, 0.3250, 0.0980; pink: 1 .2 .8 - match these up with the tasktypes variable.
    spranges = [0,0.015;0,0.005;0,0.005];
    
    for a=2:length(tasktypes)%start with a=2 to skip baselinepower
        for s=1 %:3
            subplot(4,1,sp);
            if strcmp(tasktypes{a},'switch')
                contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms_allrats_mean.(tasktypes{a}).all.all.(structures{s}).comodulogram',30,'lines','none'); %switch, plot all trials
            else
                contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms_allrats_mean.(tasktypes{a}).correct.all.(structures{s}).comodulogram',30,'lines','none'); %other tasktypes, plot correct only
            end
            ylabel('Amplitude Frequency (Hz)')
            xlabel('Phase Frequency (Hz)')
            colorbar;
            title (strcat(tasktypes{a},structures{s}));
            caxis(spranges(s,:));
            xlim([4 12]);
            sp=sp+1;
        end
    end
    
    %plot the mean comodulograms for EACH rat
    for r=1:length(allrats)
        figure('Name',strcat('Mean comodulograms for: ',allrats{r}));
        set(gcf,'renderer','Painters');
        colormap('jet');
        sp=1;
        
        %plot all the mean data
        for a=2:length(tasktypes)%start with a=2 to skip baselinepower
            for s=1 
                subplot(4,1,sp);
                if strcmp(tasktypes{a},'switch')
                    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms_allrats.(tasktypes{a}).all.all.(structures{s}).comodulogram(:,:,r)',30,'lines','none'); %switch, plot all trials
                else
                    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms_allrats.(tasktypes{a}).correct.all.(structures{s}).comodulogram(:,:,r)',30,'lines','none'); %other tasktypes, plot correct only
                end
                ylabel('Amplitude Frequency (Hz)')
                xlabel('Phase Frequency (Hz)')
                colorbar;
                title (strcat(tasktypes{a},structures{s}));
                caxis(spranges(s,:));
                xlim([4 12]);
                sp=sp+1;
            end
        end
    end
    
    %for each rat, plot mean MI for each tasktype for correct (OO/TA/OA) or all(switch) trials
    for r = 1:length(allrats) %for each rat
        figure('Name',strcat('Theta gamma PAC for: ',allrats{r}));
        for a = 2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                MeanAmp = Comodulograms_allrats.(tasktypes{a}).all.all.OBLFP.MeanAmp_gamma(r,:);
            else
                MeanAmp = Comodulograms_allrats.(tasktypes{a}).correct.all.OBLFP.MeanAmp_gamma(r,:);
            end
            plot(10:7:720,[MeanAmp,MeanAmp]/sum(MeanAmp),'Color',colors_tasktypes(a-1,:),'LineWidth',2) % was 10:20:720. needs to be twice the length of MeanAmp, which is the same length as position, which is determined by nbins
            xlim([0 720])
            set(gca,'xtick',0:360:720)
            xlabel('Phase (Deg)')
            ylabel('Amplitude')
            ylim([0.01 0.035]);
            hold on;
        end
        legend(tasktypes{2}, tasktypes{3}, tasktypes{4}, tasktypes{5})
    end
    
    %for each rat, plot mean MI for each tasktype for correct (OO/TA/OA) or all(switch) trials
    for r = 1:length(allrats) %for each rat
        figure('Name',strcat('Theta beta PAC for: ',allrats{r}));
        for a = 2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                MeanAmp = Comodulograms_allrats.(tasktypes{a}).all.all.OBLFP.MeanAmp_beta(r,:);
            else
                MeanAmp = Comodulograms_allrats.(tasktypes{a}).correct.all.OBLFP.MeanAmp_beta(r,:);
            end
            plot(10:7:720,[MeanAmp,MeanAmp]/sum(MeanAmp),'Color',colors_tasktypes(a-1,:),'LineWidth',2) % was 10:20:720. needs to be twice the length of MeanAmp, which is the same length as position, which is determined by nbins
            xlim([0 720])
            set(gca,'xtick',0:360:720)
            xlabel('Phase (Deg)')
            ylabel('Amplitude')
            ylim([0.01 0.035]);
            hold on;
        end
        legend(tasktypes{2}, tasktypes{3}, tasktypes{4}, tasktypes{5})
    end
end

%plot comodulograms for all sessions
for r=1:length(allrats)
    for f=1:length(allfiles.(allrats{r}))
        figure('Name',strcat('Session comodulogram for: ',allrats{r},', File#',num2str(f)));
        set(gcf,'renderer','Painters');
        colormap('jet');
        sp=1;
        
        %plot all the mean data
        for a=2:length(tasktypes)%start with a=2 to skip baselinepower
            for s=1 
                subplot(4,1,sp);
                if strcmp(tasktypes{a},'switch')
                    if isnan(Comodulograms_allsessions.(allrats{r}).(tasktypes{a}).all.all.(structures{s}).comodulogram(:,:,f))
                        Comodulograms_allsessions.(allrats{r}).(tasktypes{a}).all.all.(structures{s}).comodulogram(:,:,f)=zeros;
                    end
                    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms_allsessions.(allrats{r}).(tasktypes{a}).all.all.(structures{s}).comodulogram(:,:,f)',30,'lines','none'); %switch, plot all trials
                else
                    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulograms_allsessions.(allrats{r}).(tasktypes{a}).correct.all.(structures{s}).comodulogram(:,:,f)',30,'lines','none'); %other tasktypes, plot correct only
                end
                ylabel('Amplitude Frequency (Hz)')
                xlabel('Phase Frequency (Hz)')
                colorbar;
                title (strcat(tasktypes{a},structures{s}));
                caxis(spranges(s,:));
                xlim([4 12]);
                sp=sp+1;
            end
        end
    end
end

%save all figures
FigList = findobj(allchild(0),'flat','Type','figure');
savefig(FigList,'OBComodulogramFigs_07_28_21c');