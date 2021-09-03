%script to do the following, after completing trial inspection/rejection using
%the trialGUI_wavelets script:
% 1) update morlet info in all relevant structs (initially it is only
% edited in the alldata struct, but not the all4baseline or individual
% trial type structs
% 2) delete rejected trials using the rejected struct

clear;

spreadsheet = 'allfilenames_sniffing.xlsx';
[~,files, ~] = xlsread(spreadsheet, 'A2:A19');
loadfiletag = '_toneoff_morlet_rejected.mat';
savefiletag = '_toneoff_finalized.mat';

%list of structs to cycle through when deleting
structs = {'RESP','Cpoke','Rpoke','Lpoke',...
    'odo_to_cout','cout_to_choice','trialonset','blockpcor',...
    'RESP_filt','RESPmorlet','RESPmorletpeaks','RESPmorletpeaklocs','RESPmorletfreqs','RESPmorletcycles'};
morletstructs = {'RESPmorlet','RESPmorletpeaks','RESPmorletpeaklocs','RESPmorletfreqs','RESPmorletcycles'};
outcomes = {'correct','incorrect'};
trialtypes = {'congruent','incongruent','all'};
if contains(loadfiletag,'toneon')
    tonetype = 'toneon';
    tasktypes = {'toneattn','odorattn','switch'};
else
    tonetype = 'toneoff';
    tasktypes = {'odoronly','toneattn','odorattn','switch'};
end

for f=1:16 
    load(strcat(files{f},loadfiletag));
    clear wsbdata wsbtrials;
    wsbdata = wholeshebang(2:length(wholeshebang),:); %delete header row of wsb
    wsbtrials(:,1) = [wsbdata{:,2}]'; %just block/trial info
    wsbtrials(:,2) = [wsbdata{:,3}]';
    
    %each trial type struct needs RESPmorletfreqs and RESPmorlet cycles
    for i=1:length(tasktypes)
        for j=1:length(outcomes)
            for k=1:length(trialtypes)
                ntrials=length(alldata.(tonetype).(tasktypes{i}).(outcomes{j}).(trialtypes{k}).trial);
                alldata.(tonetype).(tasktypes{i}).(outcomes{j}).(trialtypes{k}).RESPmorletfreqs=NaN(ntrials,1);
                alldata.(tonetype).(tasktypes{i}).(outcomes{j}).(trialtypes{k}).RESPmorletcycles=NaN(ntrials,1);
            end
        end
    end
    %all4baseline struct needs all morlet structs
    ntrials=length(alldata.(tonetype).all4baseline.trial);
    alldata.(tonetype).all4baseline.RESPmorlet = NaN(ntrials,3663);
    alldata.(tonetype).all4baseline.RESPmorletpeaks = NaN(ntrials,100);
    alldata.(tonetype).all4baseline.RESPmorletpeaklocs = NaN(ntrials,100);
    alldata.(tonetype).all4baseline.RESPmorletfreqs = NaN(ntrials,1);
    alldata.(tonetype).all4baseline.RESPmorletcycles = NaN(ntrials,1);
    
    %% FIRST UPDATE ALL MORLET INFO FOR ALL TRIALS
    for a=1:length(alldata.(tonetype).all.trial)
        %get the current block/trial info
        currtrial = [alldata.(tonetype).all.block(a), alldata.(tonetype).all.trial(a)];
        wsbindex = intersect(find(wsbtrials(:,1)==currtrial(1)),find(wsbtrials(:,2)==currtrial(2)));
        tasktype = char(wsbdata{wsbindex,1});
        correct = char(wsbdata{wsbindex,7});
        congruent = char(wsbdata{wsbindex,6});
        
        %update info in all4baseline struct, if necessary
        all4baselineindex = intersect(find(alldata.(tonetype).all4baseline.block==currtrial(1)),...
            find(alldata.(tonetype).all4baseline.trial==currtrial(2)));
        if ~isempty(all4baselineindex)
            for m=1:length(morletstructs)
                alldata.(tonetype).all4baseline.(morletstructs{m})(all4baselineindex,:)=alldata.(tonetype).all.(morletstructs{m})(a,:);
            end
        else
        end

        %update in trial type struct
        trialtypeindex = intersect(find(alldata.(tonetype).(tasktype).(correct).(congruent).block==currtrial(1)),...
            find(alldata.(tonetype).(tasktype).(correct).(congruent).trial==currtrial(2)));
        for m=1:length(morletstructs)
            alldata.(tonetype).(tasktype).(correct).(congruent).(morletstructs{m})(trialtypeindex,:)=alldata.(tonetype).all.(morletstructs{m})(a,:);
        end
        
        %and in trial type all struct
        trialtypeallindex = intersect(find(alldata.(tonetype).(tasktype).(correct).all.block==currtrial(1)),...
            find(alldata.(tonetype).(tasktype).(correct).all.trial==currtrial(2)));
        for m=1:length(morletstructs)
            alldata.(tonetype).(tasktype).(correct).all.(morletstructs{m})(trialtypeallindex,:)=alldata.(tonetype).all.(morletstructs{m})(a,:);
        end
    end
    
    %% NOW DELETE REJECTED TRIALS
    %find all indices for trials that were rejected
    clear rejectedindex;
    rejectedindex = find(rejected.trials);
    
    %DETERMINE BLOCK/TRIAL AND WHOLESHEBANG INDEX FOR EACH REJECTED TRIAL
    nrejectedtrials = length(rejectedindex);
    trialinfo = NaN(nrejectedtrials,2);
    wsbindex = NaN(nrejectedtrials,1);

    for r=1:nrejectedtrials
        %determine the block/trial
        trialinfo(r,1)=alldata.(tonetype).all.block(rejectedindex(r));
        trialinfo(r,2)=alldata.(tonetype).all.trial(rejectedindex(r));
        %determine the corresponding wholeshebang index (for use with
        %header-less matrix)
        wsbindex(r) = intersect(find(wsbtrials(:,1)==trialinfo(r,1)),find(wsbtrials(:,2)==trialinfo(r,2)));
    end
    
    %DELETE TRIALS FROM ALLDATA STRUCT
    for r=1:nrejectedtrials
        %find index for alldata struct
        alldataindex = intersect(find(alldata.(tonetype).all.block==trialinfo(r,1)),...
            find(alldata.(tonetype).all.trial==trialinfo(r,2)));
        for s=1:length(structs)
            alldata.(tonetype).all.(structs{s})(alldataindex,:)=NaN;
        end
    end
    
    %DELETE TRIALS FROM ALL4BASELINE STRUCT(if necessary)
    for r=1:nrejectedtrials
        %find index for all4baseline struct
        all4baselineindex = intersect(find(alldata.(tonetype).all4baseline.block==trialinfo(r,1)),...
            find(alldata.(tonetype).all4baseline.trial==trialinfo(r,2)));
        %if the index is empty, do nothing, if it contains a value, delete
        %those trials
        if ~isempty(all4baselineindex)
            for s=1:length(structs)
            alldata.(tonetype).all4baseline.(structs{s})(all4baselineindex,:)=NaN;
            end
        else
        end
    end
    
    %DELETE TRIALS FROM SPECIFIC TRIAL TYPE STRUCT
    for r=1:nrejectedtrials
        tasktype = char(wsbdata{wsbindex(r),1});
        correct = char(wsbdata{wsbindex(r),7});
        congruent = char(wsbdata{wsbindex(r),6});
        %find index for trial type struct and delete
        trialtypeindex = intersect(find(alldata.(tonetype).(tasktype).(correct).(congruent).block==trialinfo(r,1)),...
            find(alldata.(tonetype).(tasktype).(correct).(congruent).trial==trialinfo(r,2)));
        for s=1:length(structs)
            alldata.(tonetype).(tasktype).(correct).(congruent).(structs{s})(trialtypeindex,:)=NaN;
        end
        
        %now do the same for the 'all' struct (which includes both
        %congruent and incongruent trials
        trialtypeallindex = intersect(find(alldata.(tonetype).(tasktype).(correct).all.block==trialinfo(r,1)),...
            find(alldata.(tonetype).(tasktype).(correct).all.trial==trialinfo(r,2)));
        for s=1:length(structs)
            alldata.(tonetype).(tasktype).(correct).all.(structs{s})(trialtypeallindex,:)=NaN;
        end
    end
    
    %add rejected struct to info struct
    alldata.info.rejected = rejected;
    alldata.info.rejected.trianinfoheader = {'block','trial'};
    alldata.info.rejected.trialinfo = trialinfo;
    
    %save adjusted file with new name
    save(strcat(files{f},savefiletag),'alldata','wholeshebang');
end