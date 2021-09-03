%% user input, import file
%Script for pulling out LFP data for individual trials and sort them by trial type,
%as well as the accompanying behavior data
%Hillary Cansler 2019

clear;
addpath(genpath('C:\TDT\TDTMatlabSDK\TDTSDK'));
addpath(genpath('C:\Users\smelluser\Documents\MATLAB\chronux_2_12'));
spreadsheet='allfilenames_sniffing.xlsx';
savefiletag='_alltrls';
%xls notes: switchblocks: OO blocks B:G; TA blocks H:M; OA blocks N:S;
%switch T:Y
[~,files, ~]=xlsread(spreadsheet, 'A2:A19');
[allOOblocks,~,~]=xlsread(spreadsheet, 'B2:G19');
[allTAblocks,~,~]=xlsread(spreadsheet, 'H2:M19');
[allOAblocks,~,~]=xlsread(spreadsheet, 'N2:S19');
[allswitchblocks,~,~]=xlsread(spreadsheet, 'T2:Y19');
tasktypes = {'odoronly','toneattn','odorattn','switch'}; 
correct = {'correct','incorrect'};
regions = {'RESP'};
pokes = {'Cpoke','Rpoke','Lpoke'};
trialinfo = {'odo_to_cout','cout_to_choice','block','trial','odor','tone','trialonset','blockpcor','RESP_filt'};
regionpokeinfo = [regions, pokes, trialinfo];
trialtypes = {'congruent','incongruent','all'};
% extractdata = 1; %1 for tone ON, 0 for tone OFF
% if extractdata==0
%     tonetype = 'toneoff';
% elseif extractdata==1
%     tonetype = 'toneon';
% end


times = [-10.4 8;...%length of data to pull out/filter (relative to trl# marker)
    -4 2]; %length of data to actually save after filtering (relative to trl# marker)
RESPfs=6.103515625000000e+02;
Pokefs=3.814697265625000e+02;
%LFPfs=3.051757812500000e+03/3;
saveindex = [round((times(2,1)+10.4)*RESPfs) round((times(2,2)+10.4)*RESPfs);... %for LFP data
    round((times(2,1)+10.4)*Pokefs) round((times(2,2)+10.4)*Pokefs)]; %for poke data

%input the channel number within Wave that correspond to each brain region.
%Put the channel you want to analyze for the OT, PFC, and OB first. This
%single channel will be used for fft, the other channel will be put into
%XXchan2 below.
%xls notes: switchblocks T:Y
RESPchan = 1;

for a=1:length(files) %FOR EACH FILE

    %note that for all multi-site LFP data files, the units for the neural data are uV
    filename = files{a};
    data=TDTbin2mat(filename);
    %load(strcat(filename,'_epocs.mat'));
    clear OOblocks TAblocks OAblocks switchblocks OOblocktimes TAblocktimes OAblocktimes switchblocktimes wholeshebang alldata;
    
    %Create variables containing block numbers only for the row
    %corresponding to the particular file being analyzed at this round of
    %the loop, while deleting any cells containing NaNs.
    

    for b=1:length(allOOblocks(a,:))
        if (allOOblocks(a,b))>0
            OOblocks(b)=allOOblocks(a,b);
        else
        end
    end
    for b=1:length(allTAblocks(a,:))
        if (allTAblocks(a,b))>0
            TAblocks(b)=allTAblocks(a,b);
        else
        end
    end
    for b=1:length(allOAblocks(a,:))
        if (allOAblocks(a,b))>0
            OAblocks(b)=allOAblocks(a,b);
        else
        end
    end
    for b=1:length(allswitchblocks(a,:))
        if(allswitchblocks(a,b))>0
            switchblocks(b)=allswitchblocks(a,b);
        else
        end
    end
    if ~exist('switchblocks')% if there were no switch blocks
        switchblocks = NaN;
    else
    end
    
    
    %Creates variables containing the start and end times of each block to be
    %analyzed, based on the user input given above about which blocks to
    %analyze for each trial type.
    OOblocktimes=zeros(length(OOblocks),2);
    for x=1:length(OOblocks)
        OOblocktimes(x, 1)=data.epocs.BlkN.onset(OOblocks(x));
        OOblocktimes(x, 2)=data.epocs.BlkN.offset(OOblocks(x));
    end

    TAblocktimes=zeros(length(TAblocks),2);
    for x=1:length(TAblocks)
        TAblocktimes(x, 1)=data.epocs.BlkN.onset(TAblocks(x));
        TAblocktimes(x, 2)=data.epocs.BlkN.offset(TAblocks(x));
    end

    OAblocktimes=zeros(length(OAblocks),2);
    for x=1:length(OAblocks)
        OAblocktimes(x, 1)=data.epocs.BlkN.onset(OAblocks(x));
        OAblocktimes(x, 2)=data.epocs.BlkN.offset(OAblocks(x));
    end
    
    switchblocktimes=zeros(length(switchblocks),2);
    if ~isnan(switchblocks(1)) %check if there are switch blocks
        for x=1:length(switchblocks)
            switchblocktimes(x,1)=data.epocs.BlkN.onset(switchblocks(x));
            switchblocktimes(x,2)=data.epocs.BlkN.offset(switchblocks(x));
        end
    else
        switchblocktimes = NaN; %if there are no switch blocks
    end

    %Create wholeshebang output, which contains matrix of all trials params. 
    
    wholeshebang = {'trial type','block','trial#','odor','tone','congruent','correct','odo to cout','cout to choice'};
    
    %initialize structure that will contain all spectrum data broken down by trial type, for each brain region.
    %initializing so I can measure the size of it to determine what row to
    %write data to, and it won't break when it doesn't exist yet. yay ugly
    %code!
    alldata.odoronly.correct.congruent.RESP=zeros(150,6106);
    alldata.odoronly.correct.congruent.Cpoke=zeros(150,2291);
    alldata.odoronly.correct.congruent.Rpoke=zeros(150,2291);
    alldata.odoronly.correct.congruent.Lpoke=zeros(150,2291);
    alldata.odoronly.correct.congruent.odo_to_cout=zeros(150,1);
    alldata.odoronly.correct.congruent.cout_to_choice=zeros(150,1);
    alldata.odoronly.correct.congruent.block=zeros(150,1);
    alldata.odoronly.correct.congruent.trial=zeros(150,1);
    alldata.odoronly.correct.congruent.odor=zeros(150,1);
    alldata.odoronly.correct.congruent.tone=zeros(150,1);
    alldata.odoronly.correct.congruent.trialonset=zeros(150,1);
    alldata.odoronly.correct.congruent.blockpcor=zeros(150,1);
    alldata.odoronly.correct.congruent.RESP_filt=zeros(150,6106);    

    alldata.odoronly.correct.incongruent=alldata.odoronly.correct.congruent;
    alldata.odoronly.correct.all=alldata.odoronly.correct.congruent;
    alldata.odoronly.incorrect=alldata.odoronly.correct;
    alldata.toneattn=alldata.odoronly;
    alldata.odorattn=alldata.odoronly;
    alldata.switch=alldata.odoronly;
    
    alldata.all.RESP=zeros(500,6106);
    alldata.all.Cpoke=zeros(500,2291);
    alldata.all.Rpoke=zeros(500,2291);
    alldata.all.Lpoke=zeros(500,2291);
    alldata.all.odo_to_cout=zeros(500,1);
    alldata.all.cout_to_choice=zeros(500,1);
    alldata.all.block=zeros(500,1);
    alldata.all.trial=zeros(500,1);
    alldata.all.odor=zeros(500,1);
    alldata.all.tone=zeros(500,1);
    alldata.all.trialonset=zeros(500,1);
    alldata.all.blockpcor=zeros(500,1);
    alldata.all.RESP_filt=zeros(500,6106);
    
    alldata.all4baseline=alldata.all;
    
    alldata.info.filename=files{a};
    alldata.info.RESPfs=[];
    alldata.info.pokefs=[];
    alldata.info.times=times;
    alldata.info.timesinfo = {'data extracted for filtering (relative to trl onset)';'data saved(relative to trl onset)'};
    alldata.info.filters=[];
    
    %% One block at a time, then one trial at a time, move through each trial,
    %pulling out trial info and data along the way, and writing to
    %wholeshebang.

    for t=1:length(tasktypes)
         if t==1
             blocktimes=OOblocktimes;
         elseif t==2
             blocktimes=TAblocktimes;
         elseif t==3
             blocktimes=OAblocktimes;
         elseif t==4
             if ~isnan(switchblocktimes)%if there are switchblocks
                blocktimes=switchblocktimes;
             else
                 blocktimes = NaN;
             end
         end

        for x=1:size(blocktimes,1) % Get all trial times for odor only blocks, one block at a time
            if ~isnan(blocktimes(1))
            data = TDTbin2mat(filename, 'T1', blocktimes(x,1), 'T2', blocktimes(x,2));

            %this adjusts things in the case of block 1, where there is a pulse
            %of all channels at the very beginning of the file. This removes
            %that "trial onset"
            if data.epocs.TrlN.onset(1)<0.01
                data.epocs.TrlN.onset=data.epocs.TrlN.onset(2:21);
                data.epocs.pCor.data=data.epocs.pCor.data(2);
            else
            end
                for y=1:length(data.epocs.TrlN.onset)                  
                    [row,col]=size(wholeshebang);
                    rownum=row+1;
                    wholeshebang{rownum,1} = tasktypes(t);
                     if t==1
                         wholeshebang{rownum,2} = OOblocks(x);
                     elseif t==2
                         wholeshebang{rownum,2} = TAblocks(x);
                     elseif t==3
                         wholeshebang{rownum,2} = OAblocks(x);
                     elseif t==4
                         wholeshebang{rownum,2} = switchblocks(x);
                     end
                    wholeshebang{rownum,3} = y;

                    onset=data.epocs.TrlN.onset(y);

                    tempOdoB=[];
                    for z=1:length(data.epocs.OdoB.onset) %subtract the trl# time from all OdoB times. If this value is greater than 0, just write 100.
                        tempOdoB(z,1)=data.epocs.OdoB.onset(z,1)-data.epocs.TrlN.onset(y,1);
                        if tempOdoB(z,1)>0
                            tempOdoB(z,1)=10000;
                        end
                        z=z+1;
                    end
                    clear odoval;
                    clear odoindex;
                    %determines the index of the OdoB onset for the trial by finding
                    %the smallest value (meaning it is the closest to the trial onset)
                    %This allows you to then use that index to determine the identity
                    %of the odor and the tone for that trial.
                    [odoval,odoindex] = min(abs(tempOdoB)); 
                    wholeshebang{rownum,4} = data.epocs.OdoB.data(odoindex);
                    wholeshebang{rownum,5} = data.epocs.Tone.data(odoindex);
                    if data.epocs.OdoB.data(odoindex)==4 && data.epocs.Tone.data(odoindex)==0
                        wholeshebang{rownum,6} = 'congruent';
                    elseif data.epocs.OdoB.data(odoindex)==1 && data.epocs.Tone.data(odoindex)==1
                        wholeshebang{rownum,6} = 'congruent';
                    else
                        wholeshebang{rownum,6} = 'incongruent';
                    end

                    %determine if the trial was correct, incorrect, or omit        
                    tempThit = [];
                    for z=1:length(data.epocs.Thit.onset) %subtract the trl# onset from all thit times
                        tempThit(z,1)=data.epocs.Thit.onset(z,1)-data.epocs.TrlN.onset(y,1);
                        z=z+1;
                    end

                    [~,hitindex]=min(abs(tempThit));

                    tempTrFA = [];
                    for z=1:length(data.epocs.TrFA.onset) %subtract the trl# onset from all miss times
                        tempTrFA(z,1)=data.epocs.TrFA.onset(z,1)-data.epocs.TrlN.onset(y,1);
                        z=z+1;
                    end
                    if isempty(tempTrFA) %if there are no incorrect trials, write 100 as a placeholder it doesn't mess up the later code.
                        tempTrFA=100;
                    else
                    end
                    [~,missindex]=min(abs(tempTrFA));

                    tempTrMI = [];
                    for z=1:length(data.epocs.TrMI.onset) %subtract the trl# onset from all omit times
                        tempTrMI(z,1)=data.epocs.TrMI.onset(z,1)-data.epocs.TrlN.onset(y,1);
                        z=z+1;
                    end
                    if isempty(tempTrMI) %if there are no incorrect trials, write 100 as a placeholder it doesn't mess up the later code.
                        tempTrMI=100;
                    else
                    end
                    [~,omitindex]=min(abs(tempTrMI));

                    %Find Cout value associated with this trial
                    tempCout = [];
                    for z=1:length(data.epocs.Cout.onset) %subtract the trl# onset from all cout times
                        tempCout(z,1)=data.epocs.Cout.onset(z,1)-data.epocs.TrlN.onset(y,1);
                        z=z+1;
                    end
                    [~, coutindex]=min(abs(tempCout));

                    %write odo to cout to wholeshebang
                    wholeshebang{rownum,8} = (data.epocs.Cout.onset(coutindex))-(data.epocs.OdoB.onset(odoindex)); %odob to cout

                    %write all min times to array, then determine the smallest min
                    %time, and write correct/incorrect/omit based on the index of the
                    %smallest value. Then calculate the time from cout to hit/miss/omit
                    %and write to wholeshebang.
                    mins=[min(abs(tempThit)), min(abs(tempTrFA)), min(abs(tempTrMI))];
                    [~, trlindex]=min(mins);
                    trloutcomes={'correct','incorrect','omit'};
                    wholeshebang{rownum,7}=trloutcomes{trlindex};
                    if trlindex==1
                        wholeshebang{rownum,9} = (data.epocs.Thit.onset(hitindex))-(data.epocs.Cout.onset(coutindex));            
                    elseif trlindex==2
                        wholeshebang{rownum,9} = (data.epocs.TrFA.onset(missindex))-(data.epocs.Cout.onset(coutindex));
                    else
                        wholeshebang{rownum,9} = (data.epocs.TrMI.onset(omitindex))-(data.epocs.Cout.onset(coutindex));
                    end

                    
                    %import LFP data
                    trialdata = TDTbin2mat(filename, 'T1', data.epocs.TrlN.onset(y)+times(1,1), 'T2', data.epocs.TrlN.onset(y)+times(1,2));

                    % downsample the data by a factor of 3 
                    %trialdata.streams.Wave.data = downsample(trialdata.streams.Wave.data',3)';
                    %trialdata.streams.Wave.fs=trialdata.streams.Wave.fs/3;
    
                    % band pass 2nd order filter 0.5-100 hz
                    %trialdata = TDTdigitalfilter(trialdata,'Wave',[0.5,100]);

                    % band stop 2nd order filter 59-61 hz
                    %trialdata=TDTdigitalfilter(trialdata,'Wave',[59 61],'TYPE','stop');
                           
                    %write the lfp data to a variable for easier referencing below
                    RESP=trialdata.streams.RESP.data(RESPchan,saveindex(1,1):saveindex(1,2));

                    Cpoke = trialdata.streams.CPok.data(saveindex(2,1):saveindex(2,2));
                    Rpoke = trialdata.streams.RPok.data(saveindex(2,1):saveindex(2,2));
                    Lpoke = trialdata.streams.LPok.data(saveindex(2,1):saveindex(2,2));
                    
                    %Filter data again, this time from 0.5-12 Hz, for use
                    %in morlet wavelet analysis after extraction
                    trialdata_theta = TDTdigitalfilter(trialdata,'RESP',[0.5 10]);
                    
                    RESP_filt = trialdata_theta.streams.RESP.data(RESPchan,saveindex(1,1):saveindex(1,2));
                    
                    %ADD SMOOTHING FUNCTION HERE?
                    %
                    %write data to a big structs, which are broken apart by
                    %tone/congruent/correct/trial type/brain region/epoc
                    if ~isempty(data.epocs.pCor.data) %calculate the percent correct for the block if it is an incomplete block
                        blockpcor = data.epocs.pCor.data;
                    else
                        blockpcor = (length(data.epocs.Thit.data)/(length(data.epocs.Thit.data)+length(data.epocs.TrFA.data)+length(data.epocs.TrMI.data)))*100;
                    end
    
                    if strcmp((trloutcomes{trlindex}),'correct') || strcmp((trloutcomes{trlindex}),'incorrect') %skip if trial was omitted        
                        if data.epocs.Tone.data(odoindex)==extractdata;
                            %put the trial data in the proper congruent/incongruent struct for the task type/outcome
                            temprow = max(find(alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).RESP(:,2)))+1;
                            if isempty (temprow)
                                temprow = 1;
                            else
                            end
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).RESP(temprow,1:length(RESP))=RESP;
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).Cpoke(temprow,1:length(Cpoke))=Cpoke;
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).Rpoke(temprow,1:length(Rpoke))=Rpoke; 
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).Lpoke(temprow,1:length(Lpoke))=Lpoke; 
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).odo_to_cout(temprow,1)=wholeshebang{rownum,8};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).cout_to_choice(temprow,1)=wholeshebang{rownum,9};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).block(temprow,1)=wholeshebang{rownum,2};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).trial(temprow,1)=wholeshebang{rownum,3};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).odor(temprow,1)=wholeshebang{rownum,4};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).tone(temprow,1)=wholeshebang{rownum,5};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).trialonset(temprow,1)=data.epocs.TrlN.onset(y);
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).blockpcor(temprow,1)=blockpcor;
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).(wholeshebang{rownum,6}).RESP_filt(temprow,1:length(RESP_filt))=RESP_filt;
                            
                            %now put the trial data in the "all" structure for the task type/outcome
                            temprow = max(find(alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.RESP(:,2)))+1;
                            if isempty (temprow)
                                temprow = 1;
                            else
                            end
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.RESP(temprow,1:length(RESP))=RESP;
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.Cpoke(temprow,1:length(Cpoke))=Cpoke;
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.Rpoke(temprow,1:length(Rpoke))=Rpoke; 
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.Lpoke(temprow,1:length(Lpoke))=Lpoke; 
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.odo_to_cout(temprow,1)=wholeshebang{rownum,8};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.cout_to_choice(temprow,1)=wholeshebang{rownum,9};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.block(temprow,1)=wholeshebang{rownum,2};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.trial(temprow,1)=wholeshebang{rownum,3};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.odor(temprow,1)=wholeshebang{rownum,4};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.tone(temprow,1)=wholeshebang{rownum,5};
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.trialonset(temprow,1)=data.epocs.TrlN.onset(y);
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.blockpcor(temprow,1)=blockpcor;
                            alldata.(tasktypes{t}).(trloutcomes{trlindex}).all.RESP_filt(temprow,1:length(RESP_filt))=RESP_filt;
                            
                            %also write data to a matrix of all trials
                            %regardless of task/trial type, for inspecting
                            %all trials
                            temprow = max(find(alldata.all.RESP(:,2)))+1;
                            if isempty (temprow)
                                temprow = 1;
                            else
                            end
                            alldata.all.RESP(temprow,1:length(RESP))=RESP;
                            alldata.all.Cpoke(temprow,1:length(Cpoke))=Cpoke;
                            alldata.all.Rpoke(temprow,1:length(Rpoke))=Rpoke;
                            alldata.all.Lpoke(temprow,1:length(Lpoke))=Lpoke;
                            alldata.all.odo_to_cout(temprow,1)=wholeshebang{rownum,8};
                            alldata.all.cout_to_choice(temprow,1)=wholeshebang{rownum,9};
                            alldata.all.block(temprow,1)=wholeshebang{rownum,2};
                            alldata.all.trial(temprow,1)=wholeshebang{rownum,3};
                            alldata.all.odor(temprow,1)=wholeshebang{rownum,4};
                            alldata.all.tone(temprow,1)=wholeshebang{rownum,5};
                            alldata.all.trialonset(temprow,1)=data.epocs.TrlN.onset(y);
                            alldata.all.blockpcor(temprow,1)=blockpcor;
                            alldata.all.RESP_filt(temprow,1:length(RESP_filt))=RESP_filt;
                            
                            %also write data to a matrix of all trials
                            %regardless of task/trial type, if the trial is
                            %correct and from a block with >=80% correct,
                            %for computing cross condition baseline
                            if strcmp(trloutcomes{trlindex}, 'correct') && round(blockpcor)>=80
                                temprow = max(find(alldata.all4baseline.RESP(:,2)))+1;
                                if isempty (temprow)
                                    temprow = 1;
                                else
                                end
                                alldata.all4baseline.RESP(temprow,1:length(RESP))=RESP;
                                alldata.all4baseline.Cpoke(temprow,1:length(Cpoke))=Cpoke;
                                alldata.all4baseline.Rpoke(temprow,1:length(Rpoke))=Rpoke;
                                alldata.all4baseline.Lpoke(temprow,1:length(Lpoke))=Lpoke;
                                alldata.all4baseline.odo_to_cout(temprow,1)=wholeshebang{rownum,8};
                                alldata.all4baseline.cout_to_choice(temprow,1)=wholeshebang{rownum,9};
                                alldata.all4baseline.block(temprow,1)=wholeshebang{rownum,2};
                                alldata.all4baseline.trial(temprow,1)=wholeshebang{rownum,3};
                                alldata.all4baseline.odor(temprow,1)=wholeshebang{rownum,4};
                                alldata.all4baseline.tone(temprow,1)=wholeshebang{rownum,5};
                                alldata.all4baseline.trialonset(temprow,1)=data.epocs.TrlN.onset(y);
                                alldata.all4baseline.blockpcor(temprow,1)=blockpcor;
                                alldata.all4baseline.RESP_filt(temprow,1:length(RESP_filt))=RESP_filt;
                            else
                            end
                        else
                        end
                    else
                    end
                end %each trial
            else
            end
        end %each block
     end %each task type
    
    %delete zero rows
    for a=1:length(tasktypes)
        for b=1:length(correct)
            for c=1:length(trialtypes)
                ntrials = max(find(alldata.(tasktypes{a}).(correct{b}).(trialtypes{c}).RESP(:,2)));
                for d=1:length(regionpokeinfo)
                    if ~isempty(ntrials)
                        alldata.(tasktypes{a}).(correct{b}).(trialtypes{c}).(regionpokeinfo{d})(ntrials+1:150,:)=[];
                    elseif isempty(ntrials)
                        alldata.(tasktypes{a}).(correct{b}).(trialtypes{c}).(regionpokeinfo{d})=[];
                    end
                end
            end
        end
    end
    
    %delete zero rows for data.(tonetype).all
    ntrials = max(find(alldata.all.RESP(:,2)));
    for d=1:length(regionpokeinfo)
        alldata.all.(regionpokeinfo{d})(ntrials+1:500,:)=[];
    end
    
    %delete zero rows for data.(tonetype).all4baseline
    ntrials = max(find(alldata.all4baseline.RESP(:,2)));
    for d=1:length(regionpokeinfo)
        alldata.all4baseline.(regionpokeinfo{d})(ntrials+1:500,:)=[];
    end
    
    alldata.info.RESPfs = trialdata.streams.RESP.fs;
    alldata.info.pokefs = trialdata.streams.CPok.fs;
    alldata.info.filters = {'acquired at 0.5-20 Hz, with 60-120-180 notch'};
    alldata.info.RESP_filt = {'0.5-10 hz for morlet convolution'};
    alldata.all.info = {'all trials with tone off'};
    alldata.all4baseline.info = {'all correct trials with tone off and block at criterion (>=80%'};
    
    
    %% save variables
    save(strcat(filename,savefiletag),'alldata','wholeshebang');
end %each file


