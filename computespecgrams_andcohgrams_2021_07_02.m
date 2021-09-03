%% calculates and saves (1) spectrograms (2) spectral power based on spectrograms
%(3) coherograms (4) coherence (5) spectral power based on mtspectrumc function 
%7.2.21 adding in code to calculate/organize things based on the outcome of
%the PREVIOUS trial, in addition to the outcome of the future trial

clear;
addpath(genpath('C:\TDT\TDTMatlabSDK\TDTSDK'));
addpath(genpath('C:\Users\smelluser\Documents\MATLAB\chronux_2_12'));

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
plotmeans = 0; %set to zero if you don't want to create the plots of the band power means
thetarange = [2 12];
betarange = [15 35];
lowgammarange = [40 60];
highgammarange = [60 80];
epochindex = [0,0;2035,3052;3053,4071;4072,4477;];

%set up struct of all files divided by rat
allrats = {'RAT064','RAT066','RAT071','RAT073','RAT074'};
allfiles.RAT064 = files(1:4);
allfiles.RAT066 = files(5:9);
allfiles.RAT071 = files(10:15);
allfiles.RAT073 = files(16:20);
allfiles.RAT074 = files(21:23);


for r = 1:length(allrats) %for each rat
    %in between rats, clear variables because some rats will have more files than others
    clear allgrams shuffled subtracted SWITCHspecgram SWITCHdbconverted incorrects;
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
        
        %delete correct trials until they equal the number of incorrect
        %trials and save in a new struture "somecorrect"
        for a = 2:length(tasktypes)
            for tt = 1:2 %for congruent and incongruent trials
                incorrtrls = size(alldata.toneoff.(tasktypes{a}).incorrect.(trialtypes{tt}).OBLFP,1); %how many incorrect trls are there
                corrtrls = size(alldata.toneoff.(tasktypes{a}).correct.(trialtypes{tt}).OBLFP,1);
                randtrls = randperm(size(alldata.toneoff.(tasktypes{a}).correct.(trialtypes{tt}).OBLFP,1)); %random numbers for total num correct/tt trls
                if incorrtrls>0 && corrtrls>=incorrtrls %if there are incorrects, and there are more corrects than inocorrects
                    for t=1:incorrtrls %using random numbers above, place data for a randomly selected correct trial of the same type in the somecorrect structure
                        if strcmp(tasktypes{a},'switch')
                            for s=1:18 %not all the sniffing related data are in the switch structure
                                alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt}).(structures2del{s})(t,:)=alldata.toneoff.(tasktypes{a}).correct.(trialtypes{tt}).(structures2del{s})(randtrls(t),:);
                            end
                        else
                            for s=1:length(structures2del)
                                alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt}).(structures2del{s})(t,:)=alldata.toneoff.(tasktypes{a}).correct.(trialtypes{tt}).(structures2del{s})(randtrls(t),:);
                            end
                        end
                    end
                elseif incorrtrls>0 && incorrtrls>corrtrls %if there are more incorrects (should only happen for switch blocks), add incorr trials until they equal corr
                    alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt})=alldata.toneoff.(tasktypes{a}).correct.(trialtypes{tt});
                else %if there are no incorrect trials of this trial type
                    if strcmp(tasktypes{a},'switch')
                        for s=1:18
                            alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt}).(structures2del{s})(1,:)=[NaN];
                            alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt}).(structures2del{s})(1,:)=[];
                        end
                    else
                        for s=1:length(structures2del)
                            alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt}).(structures2del{s})(1,:)=[NaN];
                            alldata.toneoff.(tasktypes{a}).somecorrect.(trialtypes{tt}).(structures2del{s})(1,:)=[];
                        end
                    end
                end
            end
        end
        
        %combine congruent and incongruent into 'all' in the somecorrect structure
        for a = 2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                for s=1:18
                    alldata.toneoff.(tasktypes{a}).somecorrect.all.(structures2del{s}) = [alldata.toneoff.(tasktypes{a}).somecorrect.congruent.(structures2del{s});alldata.toneoff.(tasktypes{a}).somecorrect.incongruent.(structures2del{s})];
                end
            else
                for s=1:length(structures2del)
                    alldata.toneoff.(tasktypes{a}).somecorrect.all.(structures2del{s}) = [alldata.toneoff.(tasktypes{a}).somecorrect.congruent.(structures2del{s});alldata.toneoff.(tasktypes{a}).somecorrect.incongruent.(structures2del{s})];
                end
            end
        end
        
        % for switch data, combine corrects and incorrects
        for s = 1:3
            for x=1:length(structureslong)%for switch data, combine corrects and incorrects
                alldata.toneoff.switch.all.(structureslong{x}) = [alldata.toneoff.switch.correct.all.(structureslong{x});alldata.toneoff.switch.incorrect.all.(structureslong{x})];
            end
        end
        
        %For just switch blocks, create "prevcorrect" and "previncorrect"
        %structures that organize trials based on the outcome of the
        %previous trial
        for a=5 %just switch
            switchblks = unique(alldata.toneoff.switch.all.block);
            wholeshebangswitch = ismember(cell2mat(wholeshebang(2:length(wholeshebang),2)),switchblks);
            switchtrials = find(wholeshebangswitch);
            switchtrials = switchtrials +1;
            prevcorrcount = 1;
            previncorrcount = 1;
            for t = 2:length(switchtrials)
                currtrl = [cell2mat(wholeshebang(switchtrials(t),2)),cell2mat(wholeshebang(switchtrials(t),3))]; %current trial
                curroutcome = wholeshebang{switchtrials(t),7}; %where to find the data for the current trial
                prevoutcome = wholeshebang{switchtrials(t-1),7};
                currtone = cell2mat(wholeshebang(switchtrials(t),5));
                %find the row for the current trial in the switch.all structure
                block = find(alldata.toneoff.switch.all.block==currtrl(1));
                trial = find(alldata.toneoff.switch.all.trial==currtrl(2));
                trlindex = intersect(block,trial);
                if currtone == 0 && strcmp(prevoutcome,'correct') && ~strcmp(curroutcome,'omit') %if previous trial was correct and current trial is tone off (and therefore data is avail)
                    %put all data in the switch.prevcorrect structure
                    for s=1:length(structureslong)
                        alldata.toneoff.switch.prevcorrect.all.(structureslong{s})(prevcorrcount,:) = alldata.toneoff.switch.all.(structureslong{s})(trlindex,:);
                    end
                    prevcorrcount = prevcorrcount+1;
                elseif currtone == 0 && strcmp(prevoutcome,'incorrect') && ~strcmp(curroutcome,'omit') %if previous trial was incorrect
                    %put all data in the switch.previncorrect.structure
                    for s=1:length(structureslong)
                        alldata.toneoff.switch.previncorrect.all.(structureslong{s})(previncorrcount,:) = alldata.toneoff.switch.all.(structureslong{s})(trlindex,:);
                    end
                    previncorrcount = previncorrcount+1;
                end
            end
        end
        %randomly delete trials so that n prevcorrect trials == n
        %previncorrect trials and save this in someprevcorrect structure
        for a = 5 %just switch
            if ~isempty(switchblks)
                previncorrtrls = size(alldata.toneoff.(tasktypes{a}).previncorrect.all.OBLFP,1); %how many incorrect trls are there
                prevcorrtrls = size(alldata.toneoff.(tasktypes{a}).prevcorrect.all.OBLFP,1);
                randtrls = randperm(size(alldata.toneoff.(tasktypes{a}).prevcorrect.all.OBLFP,1)); %random numbers for total num correct/tt trls
                if previncorrtrls>0 && prevcorrtrls>=previncorrtrls %if there are incorrects, and there are more corrects than inocorrects
                    for t=1:previncorrtrls %using random numbers above, place data for a randomly selected correct trial of the same type in the somecorrect structure
                        for s=1:18 %not all the sniffing related data are in the switch structure
                            alldata.toneoff.(tasktypes{a}).someprevcorrect.all.(structures2del{s})(t,:)=alldata.toneoff.(tasktypes{a}).prevcorrect.all.(structures2del{s})(randtrls(t),:);
                        end
                    end
                elseif previncorrtrls>0 && previncorrtrls>prevcorrtrls %if there are more incorrects (should only happen for switch blocks), add incorr trials until they equal corr
                    alldata.toneoff.(tasktypes{a}).someprevcorrect.all=alldata.toneoff.(tasktypes{a}).prevcorrect.all;
                else %if there are no incorrect trials of this trial type
                    for s=1:18
                        alldata.toneoff.(tasktypes{a}).someprevcorrect.all.(structures2del{s})(1,:)=[NaN];
                        alldata.toneoff.(tasktypes{a}).someprevcorrect.all.(structures2del{s})(1,:)=[];
                    end
                end
            else
            end
        end

        %compute the OB specgram for OO, TA, SWITCH, and OA, and all trials
        for s = 1:3 %for the first 3 elements in the "structures" variable: OB, PFC, OT
        
            [OOspecgram,t_spec,freq_spec] = mtspecgramc(alldata.toneoff.odoronly.correct.all.(structures{s})(:,1:6105)', specgramwindow, params);
            [TAspecgram,~,~] = mtspecgramc(alldata.toneoff.toneattn.correct.all.(structures{s})(:,1:6105)', specgramwindow, params);
            [OAspecgram,~,~] = mtspecgramc(alldata.toneoff.odorattn.correct.all.(structures{s})(:,1:6105)', specgramwindow, params);
            [all4bspecgram,~,~] = mtspecgramc(alldata.toneoff.all4baseline.(structures{s})(:,1:6105)', specgramwindow, params);
            if length(unique(alldata.toneoff.switch.all.block))>1 %some files will not have switch blocks, some will have only 1 switch block - exclude these
                [SWITCHspecgram,~,~] = mtspecgramc(alldata.toneoff.switch.all.(structures{s})(:,1:6105)', specgramwindow, params);
            elseif length(unique(alldata.toneoff.switch.all.block))==1
                SWITCHspecgram = NaN(550,100);
            elseif isempty(alldata.toneoff.switch.all.(structures{s}))
                SWITCHspecgram = NaN(550,100);
            end
            
            %subtracting 3.6 b/c time window is -4 to +2 from trial onset (which is 0.4 after odor onset)
            %now t is relative to odor onset.
            t_spec=t_spec+(alldata.info.times(2,1))+0.4;

            %convert the baseline window into an index
            [~,baseline_index(1)]=min(abs(t_spec-baseline(1)));% and convert it to an index
            [~,baseline_index(2)]=min(abs(t_spec-baseline(2)));

            %compute the baseline power based on these indices
            baseline_power = mean(all4bspecgram(baseline_index(1):baseline_index(2),:),1);
            OOdbconverted = 10*log10( bsxfun(@rdivide,OOspecgram,baseline_power));
            TAdbconverted = 10*log10( bsxfun(@rdivide,TAspecgram,baseline_power));
            OAdbconverted = 10*log10( bsxfun(@rdivide,OAspecgram,baseline_power));
            if ~isnan(SWITCHspecgram(1,1))
                SWITCHdbconverted = 10*log10( bsxfun(@rdivide,SWITCHspecgram,baseline_power));
            elseif isnan(SWITCHspecgram(1,1))
                SWITCHdbconverted = NaN(550,100);
            end
            
            allgrams.baselinepower.(strcat(structures{s},'_specgram'))(f,:) = baseline_power;
            allgrams.odoronly.specgrams.correct.(strcat(structures{s},'_specgram'))(:,:,f) = OOdbconverted;
            allgrams.toneattn.specgrams.correct.(strcat(structures{s},'_specgram'))(:,:,f) = TAdbconverted;
            allgrams.odorattn.specgrams.correct.(strcat(structures{s},'_specgram'))(:,:,f) = OAdbconverted;
            if ~isnan(SWITCHdbconverted(1,1))
                allgrams.switch.specgrams.all.(strcat(structures{s},'_specgram'))(:,:,f) = SWITCHdbconverted;
            elseif isnan(SWITCHdbconverted(1,1))
                allgrams.switch.specgrams.all.(strcat(structures{s},'specgram'))(:,:,f) = NaN(500,100);
            end
            allgrams.info = {'oo/ta/oa specgrams contain correct trials only. switch contains all tone-off trials'};
        end
        
        %average within each trial and spectral band from the specgram just created
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')
                if ~isnan(SWITCHdbconverted(1))
                for s=1:3
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).pre.theta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(133:234,thetarange(1):thetarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).hold.theta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(235:336,thetarange(1):thetarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).odor.theta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(337:378,thetarange(1):thetarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).pre.beta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(133:234,betarange(1):betarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).hold.beta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(235:336,betarange(1):betarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).odor.beta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(337:378,betarange(1):betarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).pre.lowgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(133:234,lowgammarange(1):lowgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).hold.lowgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(235:336,lowgammarange(1):lowgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).odor.lowgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(337:378,lowgammarange(1):lowgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).pre.highgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(133:234,highgammarange(1):highgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).hold.highgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(235:336,highgammarange(1):highgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).odor.highgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(337:378,highgammarange(1):highgammarange(2),f)));
                end
                end
            else
                for s=1:3
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).pre.theta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(133:234,thetarange(1):thetarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).hold.theta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(235:336,thetarange(1):thetarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).odor.theta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(337:378,thetarange(1):thetarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).pre.beta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(133:234,betarange(1):betarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).hold.beta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(235:336,betarange(1):betarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).odor.beta(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(337:378,betarange(1):betarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).pre.lowgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(133:234,lowgammarange(1):lowgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).hold.lowgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(235:336,lowgammarange(1):lowgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).odor.lowgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(337:378,lowgammarange(1):lowgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).pre.highgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(133:234,highgammarange(1):highgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).hold.highgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(235:336,highgammarange(1):highgammarange(2),f)));
                    allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).odor.highgamma(f,:) = mean(mean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(337:378,highgammarange(1):highgammarange(2),f)));
                end
            end
        end
        
        %Compute the coherence-grams for each combo of regions
        %This requires subtracted LFPs, so compute subtracted LFPs for each
        %region and save in structure "subtracted"
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch') %for switch cycle through the outcomes. for other task types use correct only
                if ~isnan(SWITCHdbconverted(1))
                for s=1:3 %three regions
                    for o=1:7 %correct, incorrect, somecorrect, all, prevcorrect, previncorrect
                        if o<4 | o>4 %correct, incorrect, all, prevcorrect, previncorrect
                            subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}) = (alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s}))-(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(strcat((structures{s}),'2')));
                        elseif o==4 %all
                            subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}) = (alldata.toneoff.(tasktypes{a}).(outcomes{o}).(structures{s}))-(alldata.toneoff.(tasktypes{a}).(outcomes{o}).(strcat((structures{s}),'2')));
                        end
                    end
                end
                end
            else
                for s=1:3 %three regions
                    for o=1:3 %for non switch, do just correct/incorrect/somecorrect
                        if ~isempty(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s}));%there are trials
                            subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}) = (alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s}))-(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(strcat((structures{s}),'2')));
                        else
                            subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}) = NaN;
                        end
                    end
                end
            end
        end
        
        %create shuffled versions of each data set and save in "shuffled"
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch') %for switch, cycle through the outcomes. for other task types use correct only
                if ~isnan(SWITCHdbconverted(1))
                for s=1:3 %three regions
                    for o=1:7
                        [m,n]=size(subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}));
                        shuffled.(tasktypes{a}).(outcomes{o}).(structures{s}) = subtracted.(tasktypes{a}).(outcomes{o}).(structures{s});
                        for i=1:m
                            idx=randperm(n);
                            shuffled.(tasktypes{a}).(outcomes{o}).(structures{s})(i,idx) = subtracted.(tasktypes{a}).(outcomes{o}).(structures{s})(1,:);
                        end
                    end
                end
                end
            else
                for s=1:3 %three regions
                    for o=1:3 %correct/incorrect/somecorrect
                        if ~isnan(subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}))
                            [m,n]=size(subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}));
                            shuffled.(tasktypes{a}).(outcomes{o}).(structures{s}) = subtracted.(tasktypes{a}).(outcomes{o}).(structures{s});
                            for i=1:m
                                idx=randperm(n);
                                shuffled.(tasktypes{a}).(outcomes{o}).(structures{s})(i,idx) = subtracted.(tasktypes{a}).(outcomes{o}).(structures{s})(1,:);
                            end
                        else
                            shuffled.(tasktypes{a}).(outcomes{o}).(structures{s}) = NaN;
                        end
                    end
                end
            end
        end
        
        %compute all cohgrams and shuffled cohgrams and save in allgrams
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch')%if switch, cycle through all outcomes
                if ~isnan(SWITCHdbconverted(1))
                    for o=1:7
                        for s=4:6 % 3 region combos in 'structures'
                            if contains(structures{s},'OBPFC')
                                [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                                [cohgramshuf1,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                                [cohgramshuf2,~,~,~,~,~,~] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                                [cohgramshufboth,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                            elseif contains(structures{s},'OBOT')
                                [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf1,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf2,~,~,~,~,~,~] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshufboth,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                            elseif contains(structures{s},'PFCOT')
                                [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf1,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf2,~,~,~,~,~,~] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshufboth,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                            end
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).noshuf.cohgram(:,:,f) = cohgram;
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf1.cohgram(:,:,f) = cohgramshuf1;
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf2.cohgram(:,:,f) = cohgramshuf2;
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shufboth.cohgram(:,:,f) = cohgramshufboth;
                        end
                    end
                elseif isnan(SWITCHdbconverted(1))
                    for o=1:7
                        for s = 4:6 % 3 region combos in 'structures'
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).noshuf.cohgram(:,:,f) = NaN(550,100);
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf1.cohgram(:,:,f) = NaN(550,100);
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf2.cohgram(:,:,f) = NaN(550,100);
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shufboth.cohgram(:,:,f) = NaN(550,100);
                        end
                    end
                end
            else
                for s=4:6 % 3 region combos in 'structures'
                    for o=1:3
                        if ~isnan(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP)
                            if contains(structures{s},'OBPFC')
                                [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                                [cohgramshuf1,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                                [cohgramshuf2,~,~,~,~,~,~] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                                [cohgramshufboth,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                            elseif contains(structures{s},'OBOT')
                                [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf1,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf2,~,~,~,~,~,~] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshufboth,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                            elseif contains(structures{s},'PFCOT')
                                [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf1,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshuf2,~,~,~,~,~,~] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                                [cohgramshufboth,~,~,~,~,~,~] = cohgramc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                            end
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).noshuf.cohgram(:,:,f) = cohgram;
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf1.cohgram(:,:,f) = cohgramshuf1;
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf2.cohgram(:,:,f) = cohgramshuf2;
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shufboth.cohgram(:,:,f) = cohgramshufboth;
                        else
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).noshuf.cohgram(:,:,f) = NaN(550,100);
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf1.cohgram(:,:,f) = NaN(550,100);
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shuf2.cohgram(:,:,f) = NaN(550,100);
                            allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).shufboth.cohgram(:,:,f) = NaN(550,100);
                        end
                    end
                end
            end
        end
        
        %Now that you have computed specgrams and cohgrams, compute the 
        %coherence within each spectral band and trial epoch for each task
        %type/region/regioncombo
        
        %compute coherence within each spectral band and trial epoch using function coherencyc
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch') %if switch, cycle through outcomes. otherwise, just correct
                if ~isnan(SWITCHdbconverted(1))
                    for s=4:6 %3 region combos: 'OBPFCcohgram','OBOTcohgram','PFCOTcohgram'
                        for o=1:7
                            %for each epoch
                            for e=2:length(trialepochs)
                                %compute the coherence
                                if contains(structures{s},'OBPFC')
                                    [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf1,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf2,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshufboth,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                elseif contains(structures{s},'OBOT')
                                    [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf1,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf2,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshufboth,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                elseif contains(structures{s},'PFCOT')
                                    [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf1,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf2,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshufboth,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                end
                                %now compute the means for each band and save to allgrams
                                if e<4
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).all(f,:) = coherence;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).theta(f) = mean(coherence(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).beta(f) = mean(coherence(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).lowgamma(f) = mean(coherence(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).highgamma(f) = mean(coherence(60:80));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).all(f,:) = coherenceshuf1;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).theta(f) = mean(coherenceshuf1(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).beta(f) = mean(coherenceshuf1(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf1(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).highgamma(f) = mean(coherenceshuf1(60:80));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).all(f,:) = coherenceshuf2;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).theta(f) = mean(coherenceshuf2(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).beta(f) = mean(coherenceshuf2(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf2(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).highgamma(f) = mean(coherenceshuf2(60:80));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).all(f,:) = coherenceshufboth;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).theta(f) = mean(coherenceshufboth(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).beta(f) = mean(coherenceshufboth(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).lowgamma(f) = mean(coherenceshufboth(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).highgamma(f) = mean(coherenceshufboth(60:80));
                                elseif e==4
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).all(f,:) = coherence;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).theta(f) = mean(coherence(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).beta(f) = mean(coherence(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).lowgamma(f) = mean(coherence(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).highgamma(f) = mean(coherence(31:40));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).all(f,:) = coherenceshuf1;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).theta(f) = mean(coherenceshuf1(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).beta(f) = mean(coherenceshuf1(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf1(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).highgamma(f) = mean(coherenceshuf1(31:40));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).all(f,:) = coherenceshuf2;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).theta(f) = mean(coherenceshuf2(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).beta(f) = mean(coherenceshuf2(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf2(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).highgamma(f) = mean(coherenceshuf2(31:40));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).all(f,:) = coherenceshufboth;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).theta(f) = mean(coherenceshufboth(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).beta(f) = mean(coherenceshufboth(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).lowgamma(f) = mean(coherenceshufboth(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).highgamma(f) = mean(coherenceshufboth(31:40));
                                end
                            end
                        end
                    end
                elseif isnan(SWITCHdbconverted(1)) %if no switch blocks, just save NaNs
                    for s=4:6 %3 region combos: 'OBPFCcohgram','OBOTcohgram','PFCOTcohgram'
                        for o=1:7
                            %for each epoch
                            for e=2:length(trialepochs)
                                %now compute the means for each band and save to allgrams
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).highgamma(f) = NaN;

                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).highgamma(f) = NaN;

                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).highgamma(f) = NaN;

                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).highgamma(f) = NaN;
                            end
                        end
                    end
                end
            else %if not switch blocks
                for o=1:3
                    if ~isnan(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP)
                        for s=4:6%for each structure/combo
                            for e=2:length(trialepochs)%for each epoch
                                %compute the coherence
                                if contains(structures{s},'OBPFC')
                                    [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf1,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf2,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshufboth,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                elseif contains(structures{s},'OBOT')
                                    [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf1,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf2,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshufboth,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                elseif contains(structures{s},'PFCOT')
                                    [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf1,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshuf2,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                    [coherenceshufboth,phi,S12,S1,S2,freq] = coherencyc(shuffled.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',shuffled.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                                end
                                %now compute the means for each band and save to allgrams
                                if e<4
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).all(f,:) = coherence;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).theta(f) = mean(coherence(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).beta(f) = mean(coherence(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).lowgamma(f) = mean(coherence(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).highgamma(f) = mean(coherence(60:80));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).all(f,:) = coherenceshuf1;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).theta(f) = mean(coherenceshuf1(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).beta(f) = mean(coherenceshuf1(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf1(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).highgamma(f) = mean(coherenceshuf1(60:80));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).all(f,:) = coherenceshuf2;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).theta(f) = mean(coherenceshuf2(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).beta(f) = mean(coherenceshuf2(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf2(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).highgamma(f) = mean(coherenceshuf2(60:80));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).all(f,:) = coherenceshufboth;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).theta(f) = mean(coherenceshufboth(2:12));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).beta(f) = mean(coherenceshufboth(15:35));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).lowgamma(f) = mean(coherenceshufboth(40:60));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).highgamma(f) = mean(coherenceshufboth(60:80));
                                elseif e==4
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).all(f,:) = coherence;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).theta(f) = mean(coherence(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).beta(f) = mean(coherence(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).lowgamma(f) = mean(coherence(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).highgamma(f) = mean(coherence(31:40));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).all(f,:) = coherenceshuf1;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).theta(f) = mean(coherenceshuf1(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).beta(f) = mean(coherenceshuf1(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf1(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).highgamma(f) = mean(coherenceshuf1(31:40));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).all(f,:) = coherenceshuf2;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).theta(f) = mean(coherenceshuf2(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).beta(f) = mean(coherenceshuf2(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).lowgamma(f) = mean(coherenceshuf2(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).highgamma(f) = mean(coherenceshuf2(31:40));

                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).all(f,:) = coherenceshufboth;
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).theta(f) = mean(coherenceshufboth(1:6));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).beta(f) = mean(coherenceshufboth(8:18));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).lowgamma(f) = mean(coherenceshufboth(21:30));
                                    allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).highgamma(f) = mean(coherenceshufboth(31:40));
                                end
                            end
                        end
                    elseif isnan(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP)
                        for s=4:6 %for each structure/combo
                            for e=2:length(trialepochs) %for each epoch
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).noshuf.(trialepochs{e}).highgamma(f) = NaN;

                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf1.(trialepochs{e}).highgamma(f) = NaN;

                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shuf2.(trialepochs{e}).highgamma(f) = NaN;

                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).shufboth.(trialepochs{e}).highgamma(f) = NaN;
                            end
                        end
                    end
                end
            end
        end
        
        %compute bandpowers using mtspectrumc in addition to the specgram
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch') %if switch, cycle through outcomes. otherwise, just correct
                if ~isnan(SWITCHdbconverted(1))
                for s=1:3 %3 region combos: 'OBPFCcohgram','OBOTcohgram','PFCOTcohgram'
                    for o=1:7
                        %for each epoch
                        for e=2:length(trialepochs)
                            %compute the power spectrum
                            if o<4 | o>4
                                power = mtspectrumc(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s})(:,epochindex(e,1):epochindex(e,2))',params);
                            elseif o==4
                               power = mtspectrumc(alldata.toneoff.(tasktypes{a}).(outcomes{o}).(structures{s})(:,epochindex(e,1):epochindex(e,2))',params);
                            end
                            %now compute the means for each band and save to allgrams
                            if e<4
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).all(f,:) = power;
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).theta(f) = mean(power(2:12));
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).beta(f) = mean(power(15:35));
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).lowgamma(f) = mean(power(40:60));
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).highgamma(f) = mean(power(60:80));
                            elseif e==4
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).all(f,:) = power;
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).theta(f) = mean(power(1:6));
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).beta(f) = mean(power(8:18));
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).lowgamma(f) = mean(power(21:30));
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).highgamma(f) = mean(power(31:40));
                            end
                        end
                    end
                end
                elseif isnan(SWITCHdbconverted(1))
                    for s=1:3 %3 region combos: 'OBPFCcohgram','OBOTcohgram','PFCOTcohgram'
                        for o=1:7
                            %for each epoch
                            for e=2:length(trialepochs)
                                %no switch trials, save nothing
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).all(f,:) = NaN;
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).theta(f) = NaN;
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).beta(f) = NaN;
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).lowgamma(f) = NaN;
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).highgamma(f) = NaN;
                             end
                        end
                    end
                end
            else
                for s=1:3 %for each structure/combo
                    for o=1:3
                        if ~isempty(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s}))
                            for e=2:length(trialepochs)%for each epoch
                                %compute the power spectrum
                                power = mtspectrumc(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.(structures{s})(:,epochindex(e,1):epochindex(e,2))',params);
                                %now compute the means for each band and save to allgrams
                                if e<4
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).all(f,:) = power;
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).theta(f) = mean(power(2:12));
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).beta(f) = mean(power(15:35));
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).lowgamma(f) = mean(power(40:60));
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).highgamma(f) = mean(power(60:80));
                                elseif e==4
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).all(f,:) = power;
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).theta(f) = mean(power(1:6));
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).beta(f) = mean(power(8:18));
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).lowgamma(f) = mean(power(21:30));
                                    allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).highgamma(f) = mean(power(31:40));
                                end
                            end
                        else
                            for e=2:length(trialepochs)
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).all(f,:) = NaN;
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).theta(f) = NaN;
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).beta(f) = NaN;
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).lowgamma(f) = NaN;
                                allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).highgamma(f) = NaN;
                            end
                        end
                    end
                end
            end
        end
    
        %Count up the number of trials contributing to each type
        for a=2:length(tasktypes)
            if strcmp(tasktypes{a},'switch') %if switch, cycle through outcomes. otherwise, just correct
                if ~isnan(SWITCHdbconverted(1))
                    for o=1:7
                        if o<4 | o>4
                            allgrams.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = length(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.trial);
                        elseif o==4
                            allgrams.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = length(alldata.toneoff.(tasktypes{a}).(outcomes{o}).trial);
                        end
                    end
                elseif isnan(SWITCHdbconverted(1))
                    for o=1:7
                        if o<4 | o>4
                            allgrams.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = NaN;
                        elseif o==4
                            allgrams.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = NaN;
                        end
                    end
                end
            else
            %for each structure/combo
            for o=1:3 %correct trials only.
                if ~isempty(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.trial)
                    allgrams.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = length(alldata.toneoff.(tasktypes{a}).(outcomes{o}).all.trial);
                else
                    allgrams.(tasktypes{a}).trlnums.(outcomes{o}).all(f) = NaN;
                end
            end
            end
        end

        %For this file, set aside the incorrects and corresponding somecorrects to be averaged across files
        for a=2:length(tasktypes)
            incorrects.(strcat('file',num2str(f))).(tasktypes{a}).incorrect = alldata.toneoff.(tasktypes{a}).incorrect.all;
            incorrects.(strcat('file',num2str(f))).(tasktypes{a}).somecorrect = alldata.toneoff.(tasktypes{a}).somecorrect.all;
        end
        

        
        
    end %end each file per rat

    %-----------------------------------------------------------------------
    %ALL FILES FOR CURRENT RAT DONE
    %COMPUTE MEAN FOR THIS RAT AND SAVE IN ALLGRAMS_ALLRATS
    %-----------------------------------------------------------------------
    
    for a=2:length(tasktypes)
        
        for s=1:3 %average the specgrams
            if strcmp(tasktypes{a},'switch')
                allgrams_allrats.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(:,:,r)=nanmean(allgrams.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram')),3);
            else
                allgrams_allrats.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(:,:,r)=nanmean(allgrams.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram')),3);
            end
        end
        
        for s=1:3 %average the band means computed from the specgram
            if strcmp(tasktypes{a},'switch')
                for e=2:length(trialepochs)
                    for b=2:length(bands)
                        allgrams_allrats.(tasktypes{a}).power_specgrams.all.(structures{s}).(trialepochs{e}).(bands{b})(r,:) = nanmean(allgrams.(tasktypes{a}).power_specgrams.all.(structures{s}).(trialepochs{e}).(bands{b}));
                    end
                end
            else
                for e=2:length(trialepochs)
                    for b=2:length(bands)
                        allgrams_allrats.(tasktypes{a}).power_specgrams.correct.(structures{s}).(trialepochs{e}).(bands{b})(r,:) = nanmean(allgrams.(tasktypes{a}).power_specgrams.correct.(structures{s}).(trialepochs{e}).(bands{b}));
                    end
                end
            end
        end

        for s=4:6 %average the cohgrams
            if strcmp(tasktypes{a},'switch') && exist('SWITCHdbconverted','var')
                for o=1:7
                    for c=1:length(cohgramstructures)
                        allgrams_allrats.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram(:,:,r) = nanmean(allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram,3);
                    end
                end
            else
                for o=1:3
                    for c=1:length(cohgramstructures)
                        allgrams_allrats.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram(:,:,r) = nanmean(allgrams.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram,3);
                    end
                end
            end
        end
        
        for s=4:6 %average the coherence values
            if strcmp(tasktypes{a},'switch')
                for o=1:7
                    for c=1:length(cohgramstructures)
                        for e=2:length(trialepochs)
                            for b=1:length(bands)
                                allgrams_allrats.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b})(r,:) = nanmean(allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}));
                            end
                        end
                    end
                end
            else
                for o=1:3
                    for c=1:length(cohgramstructures)
                        for e=2:length(trialepochs)
                            for b=1:length(bands)
                                allgrams_allrats.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b})(r,:) = nanmean(allgrams.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}));
                            end
                        end
                    end
                end
            end
        end
        
        for s=1:3 %average the band powers
            if strcmp(tasktypes{a},'switch')
                for o=1:7
                    for e=2:length(trialepochs)
                        for b=1:length(bands)
                            allgrams_allrats.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b})(r,:) = nanmean(allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}));
                        end
                    end
                end
            else
                for o=1:3
                    for e=2:length(trialepochs)
                        for b=1:length(bands)
                            allgrams_allrats.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b})(r,:) = nanmean(allgrams.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}));
                        end
                    end
                end
            end
        end
        
        %average the trl nums
        if strcmp(tasktypes{a},'switch') && exist('SWITCHdbconverted','var')
            for o=1:7
                for tt=3 %trial types - just all for now
                    allgrams_allrats.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt})(r) = nanmean(allgrams.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt}));
                end
            end
        else
            for o=1:3
                for tt=3
                    allgrams_allrats.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt})(r) = nanmean(allgrams.(tasktypes{a}).trlnums.(outcomes{o}).(trialtypes{tt}));
                end
            end
        end

    end
    
    %-----------------------------------------------------------------------
    %ALL INCORRECTS/SOMECORRECTS
    %COMBINE ACROSS FILES AND COMPUTE MEANS, THEN SAVE IN ALLGRAMS_ALLRATS
    %-----------------------------------------------------------------------
    
    %first create empty all structures to be filled with each file
    incorrects.all = incorrects.file1;
    for a=2:length(tasktypes)
        for o=2:3
            for s=1:18
                incorrects.all.(tasktypes{a}).(outcomes{o}).(structures2del{s}) = [];
            end
        end
    end
    
    %go through each file combine incorrects/somecorrects
    for a=2:length(tasktypes)
        for o=2:3 %incorrect, some correct
            for i=1:length(allfiles.(allrats{r}))
                for s=1:18
                    incorrects.all.(tasktypes{a}).(outcomes{o}).(structures2del{s}) = [incorrects.all.(tasktypes{a}).(outcomes{o}).(structures2del{s});incorrects.(strcat('file',num2str(i))).(tasktypes{a}).(outcomes{o}).(structures2del{s})];
                end
            end
        end
    end
    
    %for these corrects/incorrects, compute coherence and spectra
    %save the results under
    %allgrams_allrats.tasktype.coherence.incorrect_all/somecorrect_all
    
    for a=2:length(tasktypes) %specgrams
        for o=2:3 %incorrect/somecorrect
            for s=1:3
                [specgram,t_spec,freq_spec] = mtspecgramc(incorrects.all.(tasktypes{a}).(outcomes{o}).(structures{s})(:,1:6105)', specgramwindow, params);
                baselinepower = mean(allgrams.baselinepower.(strcat(structures{s},'_specgram')));
                specgramdbconverted = 10*log10( bsxfun(@rdivide,specgram,baseline_power));
                allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(:,:,r) = specgram;
            end
        end
    end
    
    for a=2:length(tasktypes) %band power from specgrams
        for o=2:3 %incorrect/somecorrect
            for s=1:3
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).pre.theta(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(133:234,thetarange(1):thetarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).hold.theta(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(235:336,thetarange(1):thetarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).odor.theta(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(337:378,thetarange(1):thetarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).pre.beta(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(133:234,betarange(1):betarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).hold.beta(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(235:336,betarange(1):betarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).odor.beta(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(337:378,betarange(1):betarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).pre.lowgamma(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(133:234,lowgammarange(1):lowgammarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).hold.lowgamma(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(235:336,lowgammarange(1):lowgammarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).odor.lowgamma(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(337:378,lowgammarange(1):lowgammarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).pre.highgamma(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(133:234,highgammarange(1):highgammarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).hold.highgamma(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(235:336,highgammarange(1):highgammarange(2))));
                allgrams_allrats.(tasktypes{a}).power_specgrams.(strcat(outcomes{o},'_all')).(structures{s}).odor.highgamma(r,:) = mean(mean(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(outcomes{o},'_all')).(strcat(structures{s},'_specgram'))(337:378,highgammarange(1):highgammarange(2))));
            end
        end
    end
    
    %cohgrams - first create subtracted
    clear subtracted; 
    for a=2:length(tasktypes)
        for o=2:3 %incorrect/somecorrect
            for s=1:3
                subtracted.(tasktypes{a}).(outcomes{o}).(structures{s}) = (incorrects.all.(tasktypes{a}).(outcomes{o}).(structures{s})-incorrects.all.(tasktypes{a}).(outcomes{o}).(strcat(structures{s},'2')));
            end
        end
    end
    
    for a=2:length(tasktypes)%now make the cohgrams and save
        for o=2:3 %incorrect/somecorrect
            for s=4:6
                if contains(structures{s},'OBPFC')
                    [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',specgramwindow,params);
                elseif contains(structures{s},'OBOT')
                    [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                elseif contains(structures{s},'PFCOT')
                    [cohgram,phi,S12,S1,S2,t_coh,f_coh] = cohgramc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP',specgramwindow,params);
                end
                allgrams_allrats.(tasktypes{a}).cohgrams.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.cohgram(:,:,r) = cohgram;
            end
        end
    end
    
    for a=2:length(tasktypes) %compute coherence with coherencyc and save
        for o=2:3 %incorrect, some correct
            for s=4:6
                for e=2:4 %trial epochs
                    if contains(structures{s},'OBPFC')
                        [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',params);
                    elseif contains(structures{s},'OBOT')
                        [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).OBLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                    elseif contains(structures{s},'PFCOT')
                        [coherence,phi,S12,S1,S2,freq] = coherencyc(subtracted.(tasktypes{a}).(outcomes{o}).PFCLFP(:,epochindex(e,1):epochindex(e,2))',subtracted.(tasktypes{a}).(outcomes{o}).OTLFP(:,epochindex(e,1):epochindex(e,2))',params);
                    end
                    %save
                    if e<4
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).all(r,:) = coherence;
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).theta(r) = mean(coherence(2:12));
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).beta(r) = mean(coherence(15:35));
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).lowgamma(r) = mean(coherence(40:60));
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).highgamma(r) = mean(coherence(60:80));

                    elseif e==4
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).all(r,:) = coherence;
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).theta(r) = mean(coherence(1:6));
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).beta(r) = mean(coherence(8:18));
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).lowgamma(r) = mean(coherence(21:30));
                        allgrams_allrats.(tasktypes{a}).coherence.(strcat(outcomes{o},'_all')).(structures{s}).noshuf.(trialepochs{e}).highgamma(r) = mean(coherence(31:40));
                    end
                end
            end
        end
    end
    
    for a=2:length(tasktypes) %compute bandpowers with mtspectrumc and save
        for o=2:3 %incorrect/somecorrect
            for s=1:3 %lfp structures
                %for each epoch
                for e=2:length(trialepochs)
                    %compute the power spectrum
                    power = mtspectrumc(incorrects.all.(tasktypes{a}).(outcomes{o}).(structures{s})(:,epochindex(e,1):epochindex(e,2))',params);
                    %now compute the means for each band and save to allgrams
                    if e<4
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).all(r,:) = power;
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).theta(r) = mean(power(2:12));
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).beta(r) = mean(power(15:35));
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).lowgamma(r) = mean(power(40:60));
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).highgamma(r) = mean(power(60:80));
                    elseif e==4
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).all(r,:) = power;
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).theta(r) = mean(power(1:6));
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).beta(r) = mean(power(8:18));
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).lowgamma(r) = mean(power(21:30));
                        allgrams_allrats.(tasktypes{a}).power.(strcat(outcomes{o},'_all')).(structures{s}).(trialepochs{e}).highgamma(r) = mean(power(31:40));
                    end
                end
            end
        end
    end

    for a=2:length(tasktypes) %count up trial nums
        for o=2:3
            allgrams_allrats.(tasktypes{a}).trlnums.(strcat(outcomes{o},'_all')).all(r,:) = length(incorrects.all.(tasktypes{a}).(outcomes{o}).trial);
        end
    end
    
    %Save the allgrams for the specific rat
    allgrams_allsessions.(allrats{r})=allgrams;
    %save the incorrects for this specific rat
    incorrects_allrats.(allrats{r})=incorrects;
    
end %end all rats

allgrams_allrats.info.params = params;
allgrams_allrats.info.specgrams = 'made with mtspecgramc function, baseline -2-1.7, condition averaged baseline including all correct oo/ta/oa trls';
allgrams_allrats.info.cohgrams = 'first subtracted traces, then made cohgrams with fxn cohgramc';
allgrams_allrats.info.coherence = 'computed using function coherencyc';
allgrams_allrats.info.power = 'computed using function mtspectrumc';
allgrams_allrats.info.specgrampower = 'computed by averaging from specgram';

%---------------------------------------------------------------------------
%AVERAGE ACROSS ALL RATS
%---------------------------------------------------------------------------
for a=2:length(tasktypes)
        
        for s=1:3 %average the specgrams
            if strcmp(tasktypes{a},'switch')
                for o=[4 8 9] %all, incorrect_all, somecorrect_all
                    allgrams_allrats_mean.(tasktypes{a}).specgrams.(outcomes{o}).(strcat(structures{s},'_specgram')) = nanmean(allgrams_allrats.(tasktypes{a}).specgrams.(outcomes{o}).(strcat(structures{s},'_specgram')),3);
                end
            else
                for o=[1 8 9] %correct, incorrect_all, somecorrect_all
                    allgrams_allrats_mean.(tasktypes{a}).specgrams.(outcomes{o}).(strcat(structures{s},'_specgram')) = nanmean(allgrams_allrats.(tasktypes{a}).specgrams.(outcomes{o}).(strcat(structures{s},'_specgram')),3);
                end
            end
        end
        
        for s=1:3 %average the bandpowers computed from specgrams
            if strcmp(tasktypes{a},'switch')
                for o=[4 8 9] %all, incorrect_all, somecorrect_all
                    for e=2:length(trialepochs)
                        for b=2:length(bands)
                            allgrams_allrats_mean.(tasktypes{a}).power_specgrams.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}) = mean(allgrams_allrats.(tasktypes{a}).power_specgrams.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}));
                        end
                    end
                end
            else
                for o= [1 8 9] %correct, incorrect_all, somecorrect_all
                    for e=2:length(trialepochs)
                        for b=2:length(bands)
                            allgrams_allrats_mean.(tasktypes{a}).power_specgrams.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}) = mean(allgrams_allrats.(tasktypes{a}).power_specgrams.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}));
                        end
                    end
                end
            end
        end
        
        for s=4:6 %average the cohgrams
            if strcmp(tasktypes{a},'switch')
                for o=1:9
                    if o<8
                        for c=1:length(cohgramstructures) %incorrect_all and somecorrect_all only contain "noshuf" structure
                            allgrams_allrats_mean.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram = nanmean(allgrams_allrats.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram,3);
                        end
                    else
                        for c=1
                            allgrams_allrats_mean.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram = nanmean(allgrams_allrats.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram,3);
                        end
                    end
                end
            else
                for o=[1 2 3 8 9]
                    if o<8 %incorrect_all and somecorrect_all only contain "noshuf" structure
                        for c=1:length(cohgramstructures)
                            allgrams_allrats_mean.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram = nanmean(allgrams_allrats.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram,3);
                        end
                    else
                        for c=1
                            allgrams_allrats_mean.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram = nanmean(allgrams_allrats.(tasktypes{a}).cohgrams.(outcomes{o}).(structures{s}).(cohgramstructures{c}).cohgram,3);
                        end
                    end
                end
            end
        end
        
        for s=4:6 %average the coherence values
            if strcmp(tasktypes{a},'switch')
                for o=1:9
                    if o<8 %incorrect_all and somecorrect_all only contain "noshuf" structure
                        for c=1:length(cohgramstructures)
                            for e=2:length(trialepochs)
                                for b=1:length(bands)
                                    allgrams_allrats_mean.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}) = nanmean(allgrams_allrats.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}));
                                end
                            end
                        end
                    else
                        for c=1
                            for e=2:length(trialepochs)
                                for b=1:length(bands)
                                    allgrams_allrats_mean.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}) = nanmean(allgrams_allrats.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}));
                                end
                            end
                        end

                    end
                end
            else
                for o=[1 2 3 8 9]
                    if o<8 %incorrect_all and somecorrect_all only contain noshuf structure
                        for c=1:length(cohgramstructures)
                            for e=2:length(trialepochs)
                                for b=1:length(bands)
                                    allgrams_allrats_mean.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}) = nanmean(allgrams_allrats.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}));
                                end
                            end
                        end
                    else
                        for c=1
                            for e=2:length(trialepochs)
                                for b=1:length(bands)
                                    allgrams_allrats_mean.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}) = nanmean(allgrams_allrats.(tasktypes{a}).coherence.(outcomes{o}).(structures{s}).(cohgramstructures{c}).(trialepochs{e}).(bands{b}));
                                end
                            end
                        end
                    end
                end
            end
        end
        
        for s=1:3 %average the band powers
            if strcmp(tasktypes{a},'switch')
                for o=1:9
                    for e=2:length(trialepochs)
                        for b=1:length(bands)
                            allgrams_allrats_mean.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}) = nanmean(allgrams_allrats.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}));
                        end
                    end
                end
            else
                for o=[1 2 3 8 9]
                    for e=2:length(trialepochs)
                        for b=1:length(bands)
                            allgrams_allrats_mean.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}) = nanmean(allgrams_allrats.(tasktypes{a}).power.(outcomes{o}).(structures{s}).(trialepochs{e}).(bands{b}));
                        end
                    end
                end
            end
        end
end


save('specgramsandcohgrams_6andswitchblocks_2021_07_09.mat','allgrams_allrats','allgrams_allrats_mean','allgrams_allsessions','incorrects_allrats','t_coh','t_spec');    

%----------------
% MAKE PLOTS
%----------------

if makeplots > 0
    
    %subtracting 3.6 b/c time window is -4 to +2 from trial onset (which is 0.4 after odor onset)
    %now t is relative to odor onset.
    t_coh=t_coh+(alldata.info.times(2,1))+0.4;

    %plot the mean specgram and cohgram across all rats
    figure('Name','Mean specgrams across rats');
    colormap('jet');
    set(gcf,'renderer','Painters');
    sp=1;
    spranges = [-5,3;-2,4;-3,5;0,0.25;0,0.8;0,0.25];
    plotindex = [133 439]; %plot only data from -2 to +2 from odor onset
    tasktypes = {'baselinepower','odoronly','toneattn','switch','odorattn'}; %change the order of task types so switch is plotted where I want
    colors_tasktypes = [0.4660, 0.6740, 0.1880; 0, 0.4470, 0.7410;1 .2 .8;0.8500, 0.3250, 0.0980]; %green: 0.4660, 0.6740, 0.1880; blue: 0, 0.4470, 0.7410; orange: 0.8500, 0.3250, 0.0980; pink: 1 .2 .8 - match these up with the tasktypes variable.


    %plot all the mean data
    for a=2:length(tasktypes)%start with a=2 to skip baselinepower
        for s=1:3
            subplot(4,3,sp);
            if strcmp(tasktypes{a},'switch')
                contourf(t_spec(plotindex(1):plotindex(2)),freq_spec,(allgrams_allrats_mean.(tasktypes{a}).specgrams.all.(strcat(structures{s},'_specgram'))(plotindex(1):plotindex(2),:))',20,'linecolor','none');
            else
                contourf(t_spec(plotindex(1):plotindex(2)),freq_spec,(allgrams_allrats_mean.(tasktypes{a}).specgrams.correct.(strcat(structures{s},'_specgram'))(plotindex(1):plotindex(2),:))',20,'linecolor','none');
            end
            title (strcat(tasktypes{a},structures{s}));
            colorbar;
            caxis(spranges(s,:));
            sp=sp+1;
        end
    end
    figure('Name',strcat('Mean cohgrams across rats'));
    colormap('jet');
    set(gcf,'renderer','Painters');
    sp=1;
    for a=2:length(tasktypes)
        for s=4:6 %for the cohgrams
            subplot(4,3,sp);
            if strcmp(tasktypes{a},'switch')
                contourf(t_coh(plotindex(1):plotindex(2)),f_coh(2:12),(allgrams_allrats_mean.(tasktypes{a}).cohgrams.all.(structures{s}).noshuf.cohgram(plotindex(1):plotindex(2),2:12))',20,'linecolor','none');
            else
                contourf(t_coh(plotindex(1):plotindex(2)),f_coh(2:12),(allgrams_allrats_mean.(tasktypes{a}).cohgrams.correct.(structures{s}).noshuf.cohgram(plotindex(1):plotindex(2),2:12))',20,'linecolor','none');
            end
            title (strcat(tasktypes{a},structures{s}));
            colorbar;
            caxis(spranges(s,:));
            sp=sp+1;
        end
    end

    %plot mean specgram and cohgram for each rat
    %and plot a spectrogram for the rat
    for r=1:length(allrats)
        figure('Name',strcat(allrats{r},'mean specgrams'));
        colormap('jet');
        sp=1;
        %plot all the mean data
        for a=2:length(tasktypes)%start with a=2 to skip baselinepower
            for s=1:3 %for the specgrams
                subplot(4,3,sp);
                contourf(t_spec,freq_spec,(allgrams_allrats.(tasktypes{a}).specgrams.(strcat(structures{s},'_specgram'))(:,:,r))',20,'linecolor','none');
                title (strcat(tasktypes{a},structures{s}));
                colorbar;
                caxis(spranges(s,:));
                sp=sp+1;
            end
        end
        figure('Name',strcat(allrats{r},'mean cohgrams'));
        colormap('jet');
        sp=1;
        for a=2:length(tasktypes);
            for s=4:6 %for the cohgrams
                subplot(4,3,sp);
                contourf(t_coh,f_coh,(allgrams_allrats.(tasktypes{a}).cohgrams.correct.(structures{s}).noshuf.cohgram(:,:,r))',20,'linecolor','none');
                title (strcat(tasktypes{a},structures{s}));
                colorbar;
                caxis(spranges(s,:));
                sp=sp+1;
            end
        end
    end

    %plot spectra for each rat
    figure('Name','Spectra - all rats');
    sp=1;

    for e=2:length(trialepochs)
        for s=1:3
            for a=2:length(tasktypes)
                subplot(3,3,sp);
                hold on;
                plot(allgrams_allrats_mean.(tasktypes{a}).power.correct.(structures{s}).(trialepochs{e}).all,'Color',colors_tasktypes((a-1),:));
            end
            title(strcat(trialepochs{e},structures{s}));
            sp=sp+1;
        end
    end

    %plot coherence spectra for each rat
    figure('Name','Coherence spectra - all rats');
    sp=1;

    for e=2:length(trialepochs)
        for s=4:6
            for a=2:length(tasktypes)
                subplot(3,3,sp);
                hold on;
                plot(allgrams_allrats_mean.(tasktypes{a}).coherence.correct.(structures{s}).noshuf.(trialepochs{e}).all,'Color',colors_tasktypes((a-1),:));
            end
            title(strcat(trialepochs{e},structures{s}));
            sp=sp+1;
        end
    end
    
end
