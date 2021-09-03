%script to convert detected sniff peak data into raster format

clear;

spreadsheet = 'allfilenames_sniffing.xlsx';
[~,files, ~] = xlsread(spreadsheet, 'A2:A19');
loadfiletag = 'RAT137-210327-095904_toneoff_morlet_rejected.mat';
savefiletag = 'toneon_allblocks_finalized2_sniffs';

structs = {'OBLFP','PFCLFP','OTLFP','OBLFP2','PFCLFP2','OTLFP2','Cpoke','Rpoke','Lpoke',...
    'odo_to_cout','cout_to_choice','block','trial','odor','tone','trialonset','blockpcor',...
    'OBLFP_theta','OBmorlet','OBmorletpeaks','OBmorletpeaklocs','OBmorletfreqs','OBmorletcycles'};
morletstructs = {'OBmorlet','OBmorletpeaks','OBmorletpeaklocs','OBmorletfreqs','OBmorletcycles'};
outcomes = {'correct','incorrect'};
trialtypes = {'congruent','incongruent','all'};

if contains(loadfiletag,'toneon')
    tonetype = 'toneon';
    tasktypes = {'toneattn','odorattn'};
else
    tonetype = 'toneoff';
    tasktypes = {'odoronly','toneattn','odorattn'};
end

%setup from Marie's script
%--------------------------------------------------
% initialize important variables
%--------------------------------------------------
 
%---set window and bin sizes---
psthwin=[-2000 2000];       %ms, window of spike times to look at, relative to alignment choice; used to calc regular bins, bkgdsize, sdf. Needs to encompass time ranges of psth_bkgd and window FR comparisons
psthwin_small=[-2000 2000];   %ms, allows option to have diff psth sizes based on specific alignment choice
psthwin_bkgd=[-2000 -1800];  %ms, defines bkgd window; needs to be within range of psthwin; can have option to calc relative to a diff bkgd alignment option for better measure; 
%psthwin_bkgd used for BKGD_SDF(to calc better "bkgd" that might need to be adjusted/aligned differently), some ROC options, aveFR and t-tests that compare to bkgd
plot_psthwin=[-2000 2000];  %ms, size of PSTH to plot in figure, relative to alignment choice
binsize=50;         %ms, bin size for psth and plotting
nbins = ((psthwin(2)-psthwin(1))/binsize)+1;
binsize_SDF=10;     %ms, optional different bin size to use for spike density function (if needed)
%binsize_licks=200;  %ms, bin size for lick SDF (if needed)
%BaselinePeriod=600; %ms, for restricting significance testing in ROC (if needed)

%---setup for SDM analysis---
SaveAllSpkTrains=0;   %will save spike trains for later use in SDM analysis
UseLicksForSpikes=0;  %will save timestamps of licks instead of spikes
RunSpecificFiles=0;   %need to load xls sheet of which text files to run
nConditions=7;        %# of stimuli separated in xls sheet
AllSpkTrains={};      %initialize cell array
AllSpkTrains_index=0; %initialize index to zero
SDMlower=0; %sets lower limit for spiketimes to save, in sec
SDMupper=1; %sets upper limit for spiketimes to save, in sec
maxTrials=300; %used to initializeAllSpkTrains to NaN
AllSpkPrefix='AllSpkTrains_'; %for output filename
%read in selected units for SpikeTrain array analysis (needs to be ordered
%by UnitNumber and then by StimulusType
%[NUMorig,TXTorig,ORIG]=xlsread('All3at2mag_ResponsiveLickFiles.xlsx',1,'A2:AF78');

minTrials = 2;

for f=12
    load(strcat(files{f},loadfiletag));
    
    %% for each trial, convert times to desired format
    %ALL
    for t=1:length(alldata.(tonetype).all.trial)
        sniffs = alldata.(tonetype).all.OBmorletpeaklocs(t,:);
        snifftimes = alldata.(tonetype).all.OBmorletpeaklocs(t,1:(find(isnan(sniffs),1)-1));
        snifftimes = snifftimes.*1000; %convert to ms for compatability with fxns below
        alldata.(tonetype).all.sniffs(t,1) = {[snifftimes]};
    end
    
    %ALL4BASELINE
    for t=1:length(alldata.(tonetype).all4baseline.trial)
        sniffs = alldata.(tonetype).all4baseline.OBmorletpeaklocs(t,:);
        snifftimes = alldata.(tonetype).all4baseline.OBmorletpeaklocs(t,1:(find(isnan(sniffs),1)-1));
        snifftimes = snifftimes.*1000; %convert to ms for compatability with fxns below
        alldata.(tonetype).all4baseline.sniffs(t,1) = {[snifftimes]};
    end
    
    %TRIAL TYPES
    for a=1:length(tasktypes)
        for b=1:length(outcomes)
            for c=1:length(trialtypes)
                for t=1:length(alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).trial)
                    sniffs = alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).OBmorletpeaklocs(t,:);
                    snifftimes = alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).OBmorletpeaklocs(t,1:(find(isnan(sniffs),1)-1));
                    snifftimes = snifftimes.*1000; %convert to ms for compatability with fxns below
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).sniffs(t,1)={[snifftimes]};
                end
            end
        end
    end
 
    
%---------------------------------------------------------
                % COMPUTE PSTH & SPIKE DENSITY FUNCTION 
                % finds spike array to use based on set decision logic 
%---------------------------------------------------------
%----------------------------------------------    
    %FOR ALL TRIALS
%----------------------------------------------
    nSweeps = length(alldata.(tonetype).all.trial);
    %only run analysis on TrialTypes that have at least X sweeps
    if nSweeps >=minTrials
        %need to clear these local variables. (If previous array was larger than next, will not overwrite with empty values)
        clear spikearray lickarray PSTH tbins PSTH_licks tbins_licks spikearray_reordered lickarray_reordered spk spk_licks...
            sdf sdf_licks tmpSDF tmpSDF_licks tmpSDF_10ms tmpSDF_licks_10ms SDF_allTrials SDF_allTrials_10ms SDFsmooth_allTrials...
            SDF_licks_allTrials SDFsmooth_licks_allTrials PSTH PSTH_licks sdf_mean SDF_meanminus95CI SDF_meanplus95CI sdf_mean_licks...
            SDF_licks_meanminus95CI SDF_licks_meanplus95CI  tmpSDF_licks_200ms SDF_licks_allTrials_200ms sdf_mean_licks_200ms...
            SDF_meanplus2std SDF_meanminus2std SDF_meanplusSEM SDF_meanminusSEM

        spikearray = alldata.(tonetype).all.sniffs;
        %array for sdf                    
        spk=alldata.(tonetype).all.OBmorletpeaklocs;
        spk = spk.*1000; %convert to ms for use below

        %loop through sweeps to get PSTH and SDF trial-by-trial
        for t=1:nSweeps
            %----------------------------------------------
            % calculate PSTH/SDF for spikes
            %----------------------------------------------  
            % compute psth
            [PSTH(t,1:nbins), tbins] = psthMG(spikearray{t}, binsize, psthwin);

            %run SDF, convolving with function that resembles a postsynaptic potential on each trial
            %function reads in (spikearray, align value, time bins,kernal type)
            %kernal type: 1=PSP, 2=Gaussian
            %output is in 1ms resolution
            sdf(t,:)=SDFConv(spk(t,:),0,[psthwin(1) (psthwin(2)-1)],1); %PSP

            %average SDF to Binsize resolution
            %takes average of X bins and then shifts to next X
            %(e.g. 50 ms)
            Tstart=1;
            Tend=0;
            nAveBins=size(sdf,2)/binsize;
            for n=1:nAveBins
                Tend=Tend+binsize;
                tmpSDF(n)=mean(sdf(t, Tstart:Tend));
                Tstart=Tstart+binsize;
            end

            %saves SDF average result for each sweep (spikes)
            SDF_allTrials(t,:)=tmpSDF;
            SDFsmooth_allTrials(t,:)=slidingavg(SDF_allTrials(t,:),5);

        end %end all trials 
        %----------------------------------------------
        % Average across trials for spikes
        %----------------------------------------------
        %compute means
        sdf_mean = mean(SDF_allTrials);
        sdfsmooth_mean = mean(SDFsmooth_allTrials);

        %compute 95% CI, STD and SEM computed below upon
        %placing in matrix
        [mu,sigma,muci,sigmaci] = normfit(SDF_allTrials); %based on 50-ms bins
        [~,~,muci_smooth,~] = normfit(SDFsmooth_allTrials); %based on smooth data

        %save sniff data and PSTH
        alldata.(tonetype).all.sniffPSTH = PSTH;     

        %save all info to a new variable in the alldata struct
        %50ms bin data
        alldata.(tonetype).all.SDF_allTrials = SDF_allTrials;
        alldata.(tonetype).all.SDF_mean = mean(SDF_allTrials);
        alldata.(tonetype).all.SDF_meanpm95CI = muci;
        alldata.(tonetype).all.SDF_meanpm2std(1,:) = sdf_mean-(2*std(SDF_allTrials));
        alldata.(tonetype).all.SDF_meanpm2std(2,:) = sdf_mean+(2*std(SDF_allTrials));
        alldata.(tonetype).all.SDF_meanpmSEM(1,:) = sdf_mean-(std(SDF_allTrials)/sqrt(size(SDF_allTrials,1)));
        alldata.(tonetype).all.SDF_meanpmSEM(2,:) = sdf_mean+(std(SDF_allTrials)/sqrt(size(SDF_allTrials,1)));

        %smooth data
        alldata.(tonetype).all.SDFsmooth_allTrials = SDFsmooth_allTrials;
        alldata.(tonetype).all.SDFsmooth_mean = mean(SDFsmooth_allTrials);
        alldata.(tonetype).all.SDFsmooth_meanpm95CI = muci_smooth;
        alldata.(tonetype).all.SDFsmooth_meanpm2std(1,:) = sdfsmooth_mean-(2*std(SDFsmooth_allTrials));
        alldata.(tonetype).all.SDFsmooth_meanpm2std(2,:) = sdfsmooth_mean+(2*std(SDFsmooth_allTrials));
        alldata.(tonetype).all.SDFsmooth_meanpmSEM(1,:) = sdfsmooth_mean-(std(SDFsmooth_allTrials)/sqrt(size(SDFsmooth_allTrials,1)));
        alldata.(tonetype).all.SDFsmooth_meanpmSEM(2,:) = sdfsmooth_mean+(std(SDFsmooth_allTrials)/sqrt(size(SDFsmooth_allTrials,1)));
        
    end %end if statement about min trials
    
%----------------------------------------------    
    %FOR ALL4BASELINE TRIALS
%----------------------------------------------    
    nSweeps = length(alldata.(tonetype).all4baseline.trial);
    %only run analysis on TrialTypes that have at least X sweeps
    if nSweeps >=minTrials
        %need to clear these local variables. (If previous array was larger than next, will not overwrite with empty values)
        clear spikearray lickarray PSTH tbins PSTH_licks tbins_licks spikearray_reordered lickarray_reordered spk spk_licks...
            sdf sdf_licks tmpSDF tmpSDF_licks tmpSDF_10ms tmpSDF_licks_10ms SDF_allTrials SDF_allTrials_10ms SDFsmooth_allTrials...
            SDF_licks_allTrials SDFsmooth_licks_allTrials PSTH PSTH_licks sdf_mean SDF_meanminus95CI SDF_meanplus95CI sdf_mean_licks...
            SDF_licks_meanminus95CI SDF_licks_meanplus95CI  tmpSDF_licks_200ms SDF_licks_allTrials_200ms sdf_mean_licks_200ms...
            SDF_meanplus2std SDF_meanminus2std SDF_meanplusSEM SDF_meanminusSEM

        spikearray = alldata.(tonetype).all4baseline.sniffs;
        %array for sdf                    
        spk=alldata.(tonetype).all4baseline.OBmorletpeaklocs;
        spk = spk.*1000; %convert to ms for use below

        %loop through sweeps to get PSTH and SDF trial-by-trial
        for t=1:nSweeps
            %----------------------------------------------
            % calculate PSTH/SDF for spikes
            %----------------------------------------------  
            % compute psth
            [PSTH(t,1:nbins), tbins] = psthMG(spikearray{t}, binsize, psthwin);

            %run SDF, convolving with function that resembles a postsynaptic potential on each trial
            %function reads in (spikearray, align value, time bins,kernal type)
            %kernal type: 1=PSP, 2=Gaussian
            %output is in 1ms resolution
            sdf(t,:)=SDFConv(spk(t,:),0,[psthwin(1) (psthwin(2)-1)],1); %PSP

            %average SDF to Binsize resolution
            %takes average of X bins and then shifts to next X
            %(e.g. 50 ms)
            Tstart=1;
            Tend=0;
            nAveBins=size(sdf,2)/binsize;
            for n=1:nAveBins
                Tend=Tend+binsize;
                tmpSDF(n)=mean(sdf(t, Tstart:Tend));
                Tstart=Tstart+binsize;
            end

            %saves SDF average result for each sweep (spikes)
            SDF_allTrials(t,:)=tmpSDF;
            SDFsmooth_allTrials(t,:)=slidingavg(SDF_allTrials(t,:),5);

        end %end all trials 
        %----------------------------------------------
        % Average across trials for spikes
        %----------------------------------------------
        %compute means
        sdf_mean = mean(SDF_allTrials);
        sdfsmooth_mean = mean(SDFsmooth_allTrials);

        %compute 95% CI, STD and SEM computed below upon
        %placing in matrix
        [mu,sigma,muci,sigmaci] = normfit(SDF_allTrials); %based on 50-ms bins
        [~,~,muci_smooth,~] = normfit(SDFsmooth_allTrials); %based on smooth data

        %save sniff data and PSTH
        alldata.(tonetype).all4baseline.sniffPSTH = PSTH;     

        %save all info to a new variable in the alldata struct
        %50ms bin data
        alldata.(tonetype).all4baseline.SDF_allTrials = SDF_allTrials;
        alldata.(tonetype).all4baseline.SDF_mean = mean(SDF_allTrials);
        alldata.(tonetype).all4baseline.SDF_meanpm95CI = muci;
        alldata.(tonetype).all4baseline.SDF_meanpm2std(1,:) = sdf_mean-(2*std(SDF_allTrials));
        alldata.(tonetype).all4baseline.SDF_meanpm2std(2,:) = sdf_mean+(2*std(SDF_allTrials));
        alldata.(tonetype).all4baseline.SDF_meanpmSEM(1,:) = sdf_mean-(std(SDF_allTrials)/sqrt(size(SDF_allTrials,1)));
        alldata.(tonetype).all4baseline.SDF_meanpmSEM(2,:) = sdf_mean+(std(SDF_allTrials)/sqrt(size(SDF_allTrials,1)));

        %smooth data
        alldata.(tonetype).all4baseline.SDFsmooth_allTrials = SDFsmooth_allTrials;
        alldata.(tonetype).all4baseline.SDFsmooth_mean = mean(SDFsmooth_allTrials);
        alldata.(tonetype).all4baseline.SDFsmooth_meanpm95CI = muci_smooth;
        alldata.(tonetype).all4baseline.SDFsmooth_meanpm2std(1,:) = sdfsmooth_mean-(2*std(SDFsmooth_allTrials));
        alldata.(tonetype).all4baseline.SDFsmooth_meanpm2std(2,:) = sdfsmooth_mean+(2*std(SDFsmooth_allTrials));
        alldata.(tonetype).all4baseline.SDFsmooth_meanpmSEM(1,:) = sdfsmooth_mean-(std(SDFsmooth_allTrials)/sqrt(size(SDFsmooth_allTrials,1)));
        alldata.(tonetype).all4baseline.SDFsmooth_meanpmSEM(2,:) = sdfsmooth_mean+(std(SDFsmooth_allTrials)/sqrt(size(SDFsmooth_allTrials,1)));

    end %end if statement about min trials
%----------------------------------------------
    %FOR EACH TASK/OUTCOME/TRIAL TYPES
%----------------------------------------------    
    for a=1:length(tasktypes)
        for b=1:length(outcomes)
            for c=1:length(trialtypes)
                nSweeps = length(alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).trial);
                %only run analysis on TrialTypes that have at least X sweeps
                if nSweeps >=minTrials
                    %need to clear these local variables. (If previous array was larger than next, will not overwrite with empty values)
                    clear spikearray lickarray PSTH tbins PSTH_licks tbins_licks spikearray_reordered lickarray_reordered spk spk_licks...
                        sdf sdf_licks tmpSDF tmpSDF_licks tmpSDF_10ms tmpSDF_licks_10ms SDF_allTrials SDF_allTrials_10ms SDFsmooth_allTrials...
                        SDF_licks_allTrials SDFsmooth_licks_allTrials PSTH PSTH_licks sdf_mean SDF_meanminus95CI SDF_meanplus95CI sdf_mean_licks...
                        SDF_licks_meanminus95CI SDF_licks_meanplus95CI  tmpSDF_licks_200ms SDF_licks_allTrials_200ms sdf_mean_licks_200ms...
                        SDF_meanplus2std SDF_meanminus2std SDF_meanplusSEM SDF_meanminusSEM
                    
                    spikearray = alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).sniffs;
                    %array for sdf                    
                    spk=alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).OBmorletpeaklocs;
                    spk = spk.*1000; %convert to ms for use below
                    
                    %loop through sweeps to get PSTH and SDF trial-by-trial
                    for t=1:nSweeps
                        %----------------------------------------------
                        % calculate PSTH/SDF for spikes
                        %----------------------------------------------  
                        % compute psth
                        [PSTH(t,1:nbins), tbins] = psthMG(spikearray{t}, binsize, psthwin);

                        %run SDF, convolving with function that resembles a postsynaptic potential on each trial
                        %function reads in (spikearray, align value, time bins,kernal type)
                        %kernal type: 1=PSP, 2=Gaussian
                        %output is in 1ms resolution
                        sdf(t,:)=SDFConv(spk(t,:),0,[psthwin(1) (psthwin(2)-1)],1); %PSP
                        
                        %average SDF to Binsize resolution
                        %takes average of X bins and then shifts to next X
                        %(e.g. 50 ms)
                        Tstart=1;
                        Tend=0;
                        nAveBins=size(sdf,2)/binsize;
                        for n=1:nAveBins
                            Tend=Tend+binsize;
                            tmpSDF(n)=mean(sdf(t, Tstart:Tend));
                            Tstart=Tstart+binsize;
                        end

                        %saves SDF average result for each sweep (spikes)
                        SDF_allTrials(t,:)=tmpSDF;
                        SDFsmooth_allTrials(t,:)=slidingavg(SDF_allTrials(t,:),5);
                        
                    end %end all trials 
                    %----------------------------------------------
                    % Average across trials for spikes
                    %----------------------------------------------
                    %compute means
                    sdf_mean = mean(SDF_allTrials);
                    sdfsmooth_mean = mean(SDFsmooth_allTrials);
                    
                    %compute 95% CI, STD and SEM computed below upon
                    %placing in matrix
                    [mu,sigma,muci,sigmaci] = normfit(SDF_allTrials); %based on 50-ms bins
                    [~,~,muci_smooth,~] = normfit(SDFsmooth_allTrials); %based on smooth data

                    %save sniff data and PSTH
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).sniffPSTH = PSTH;     
                    
                    %save all info to a new variable in the alldata struct
                    %50ms bin data
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_allTrials = SDF_allTrials;
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_mean = mean(SDF_allTrials);
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_meanpm95CI = muci;
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_meanpm2std(1,:) = sdf_mean-(2*std(SDF_allTrials));
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_meanpm2std(2,:) = sdf_mean+(2*std(SDF_allTrials));
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_meanpmSEM(1,:) = sdf_mean-(std(SDF_allTrials)/sqrt(size(SDF_allTrials,1)));
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDF_meanpmSEM(2,:) = sdf_mean+(std(SDF_allTrials)/sqrt(size(SDF_allTrials,1)));
                    
                    %smooth data
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_allTrials = SDFsmooth_allTrials;
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_mean = mean(SDFsmooth_allTrials);
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_meanpm95CI = muci_smooth;
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_meanpm2std(1,:) = sdfsmooth_mean-(2*std(SDFsmooth_allTrials));
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_meanpm2std(2,:) = sdfsmooth_mean+(2*std(SDFsmooth_allTrials));
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_meanpmSEM(1,:) = sdfsmooth_mean-(std(SDFsmooth_allTrials)/sqrt(size(SDFsmooth_allTrials,1)));
                    alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).SDFsmooth_meanpmSEM(2,:) = sdfsmooth_mean+(std(SDFsmooth_allTrials)/sqrt(size(SDFsmooth_allTrials,1)));

                end %end if statement about min trials
                
            end %end trial types
        end %end outcomes
    end %end task types

                
    %save new info
    save(strcat(files{f},savefiletag),'alldata','wholeshebang');

end %end all files