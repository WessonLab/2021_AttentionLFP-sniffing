%previous script 'deletetrials' replaced all rejected trials with NaNs.
%This was done so the 4blocks files could be updated based on the info in
%all blocks
%this fixes that so they are fully deleted now.

clear;

spreadsheet = 'allfilenames_sniffing.xlsx';
[~,files, ~] = xlsread(spreadsheet, 'A2:A19');
loadfiletag = '_toneoff_finalized.mat';
savefiletag = '_toneoff_finalized2.mat';

structs = {'RESP','Cpoke','Rpoke','Lpoke','odo_to_cout','cout_to_choice',...
    'block','trial','odor','tone','trialonset','blockpcor',...
    'RESP_filt','RESPmorlet','RESPmorletpeaks','RESPmorletpeaklocs','RESPmorletfreqs','RESPmorletcycles'};
outcomes = {'correct','incorrect'};
trialtypes = {'congruent','incongruent','all'};
if contains(loadfiletag,'toneon')
    tonetype = 'toneon';
    tasktypes = {'toneattn','odorattn','switch'};
else
    tonetype = 'toneoff';
    tasktypes = {'odoronly','toneattn','odorattn','switch'};
end

for f=1:11
    load(strcat(files{f},loadfiletag));

    %FOR ALL TRIALS
    %find index of deleted trial by finding NaNs in the trial onset
    %struct(could be any)
    deleteindex = find(isnan(alldata.(tonetype).all.trialonset));
    %subtract 1 to the second item in deleteindex, 2 to the 3rd, and so
    %on, so that the indices remain accurate as previous trials are
    %deleted
    for n=1:length(deleteindex)
        deleteindex(n)=deleteindex(n)-(n-1);
    end
    %delete the trials using the modified index info
    for n=1:length(deleteindex)
        for s=1:length(structs)
            alldata.(tonetype).all.(structs{s})(deleteindex(n),:) = [];
        end
    end

    %FOR ALL4BASELINE
    %find index of deleted trial by finding NaNs in the trial onset
    %struct(could be any)
    deleteindex = find(isnan(alldata.(tonetype).all4baseline.trialonset));
    %subtract 1 to the second item in deleteindex, 2 to the 3rd, and so
    %on, so that the indices remain accurate as previous trials are
    %deleted
    for n=1:length(deleteindex)
        deleteindex(n)=deleteindex(n)-(n-1);
    end
    %delete the trials using the modified index info
    for n=1:length(deleteindex)
        for s=1:length(structs)
            alldata.(tonetype).all4baseline.(structs{s})(deleteindex(n),:) = [];
        end
    end

    %FOR TASK TYPE STRUCTS
    for a=1:length(tasktypes)
        for b=1:length(outcomes)
            for c=1:length(trialtypes)

                %find index of deleted trial by finding NaNs in the trial onset
                %struct(could be any)
                deleteindex = find(isnan(alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).trialonset));
                %subtract 1 to the second item in deleteindex, 2 to the 3rd, and so
                %on, so that the indices remain accurate as previous trials are
                %deleted
                for n=1:length(deleteindex)
                    deleteindex(n)=deleteindex(n)-(n-1);
                end
                %delete the trials using the modified index info
                for n=1:length(deleteindex)
                    for s=1:length(structs)
                        alldata.(tonetype).(tasktypes{a}).(outcomes{b}).(trialtypes{c}).(structs{s})(deleteindex(n),:) = [];
                    end
                end
            end
        end
    end
    
    %save things
    save(strcat(files{f},savefiletag),'alldata', 'wholeshebang');
    
end %end all files