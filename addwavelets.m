%add wavelet convolution of OB traces to alldata structs, and detect peaks
%of the convolved signal for use in sniffing analysis

%wavelet convolution code was modified from chapter 12 of the analyzing
%neural time series data book

clear;
addpath(genpath('C:\TDT\TDTMatlabSDK\TDTSDK'));
addpath(genpath('C:\Users\smelluser\Documents\MATLAB\chronux_2_12'));


spreadsheet='allfilenames_sniffing.xlsx';
[~,files, ~]=xlsread(spreadsheet, 'A2:A19');
loadfiletag='toneoff.mat';
savefiletag='toneoff_morlet.mat';

%structure info for cycling through structure below

correct = {'correct','incorrect'};
trialtypes = {'congruent','incongruent','all'};


%% cycle through each session/trial and create the convolved morlet wavelet
%note that for 3 of the 5 rats, the OB trace is inverted prior to
%convolution and peak detection.

%% Wavelet convolution with inverted OBLFP for rats 066, 073, and 074,
% to ensure that peaks detected reflect inhale peaks rather than exhale
% troughs (this is based on previous sniffing analysis done in spike2, 
% where I chose these rats to invert their OB traces for peak detection)

somefiles = [files(5:9);files(16:20)]; %rats 066 and 073

for h = 1:length(somefiles)
    load(strcat(somefiles{h},loadfiletag));
    if contains(strcat(somefiles{h},loadfiletag),'toneon')
        tonetype = 'toneon';
        tasktypes = {'toneattn','odorattn'};
    else
        tonetype = 'toneoff';
        tasktypes = {'odoronly','toneattn','odorattn'};
    end
    % LFP data from OB for one trial
    LFPtime = -3.6:1/alldata.info.LFPfs:2.4;

    f=8; %frequency of morlet wavelet
    n = 3; %number of cycles
    s = n/(2*pi*f); 
    stext = '3/(2*pi*f)'; %for saving in data structure below
    time = -1:1/alldata.info.LFPfs:1;
    sine_wave = exp(1i*2*pi*f.*time);
    gaussian_win = exp(-time.^2./(2*s^2));
    wavelet = sine_wave .* gaussian_win;
    % half of the wavelet size, useful for chopping off edges after convolution.
    halfwaveletsize = ceil(length(wavelet)/2);
    
    for a=1:length(tasktypes)
        for b=1:length(correct)
            for c=1:length(trialtypes)
                %multiply currdata by -1, to invert the trace and detect
                %troughs as peaks
                currdata = (alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBLFP_theta(:,1:6104));
                currdata = currdata.*(-1);
                if ~isempty(currdata)
                    %create matrix where the wavelet convolved data will go
                    alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorlet=zeros(length(currdata(:,1)),6104);
                    alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaks=NaN(length(currdata(:,1)),100);
                    alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaklocs=NaN(length(currdata(:,1)),100);
                    for t=1:length(currdata(:,1))
                        % convolve with data
                        % compute Gaussian
                        n_conv = length(wavelet) + length(currdata) - 1;
                        
                        fft_w = fft(wavelet,n_conv);
                        fft_e = fft(currdata(t,:),n_conv);
                        ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
                        wavelet_conv_data = real(ift(halfwaveletsize:end-halfwaveletsize+1));
                        
                        [peaks,locs] = findpeaks(wavelet_conv_data,LFPtime);
                        
                        alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorlet(t,:)=wavelet_conv_data;
                        alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaks(t,1:(length(peaks)))=peaks;
                        alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaklocs(t,1:(length(locs)))=locs;
                    end
                else
                end
            end
        end
    end %end tasktypes
    
    %repeat the same thing for the alldata struct
    %multiply currdata by -1, to invert the trace and detect
    %troughs as peaks
    currdata = (alldata.(tonetype).all.OBLFP_theta(:,1:6104));
    currdata = currdata.*(-1);
    %create matrix where the wavelet convolved data will go
    alldata.(tonetype).all.OBmorlet=zeros(length(currdata(:,1)),6104);
    alldata.(tonetype).all.OBmorletpeaks=NaN(length(currdata(:,1)),100);
    alldata.(tonetype).all.OBmorletpeaklocs=NaN(length(currdata(:,1)),100);
    for t=1:length(currdata(:,1))
        % convolve with data
        % compute Gaussian
        n_conv = length(wavelet) + length(currdata) - 1;

        fft_w = fft(wavelet,n_conv);
        fft_e = fft(currdata(t,:),n_conv);
        ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
        wavelet_conv_data = real(ift(halfwaveletsize:end-halfwaveletsize+1));

        [peaks,locs] = findpeaks(wavelet_conv_data,LFPtime);

        alldata.(tonetype).all.OBmorlet(t,:)=wavelet_conv_data;
        alldata.(tonetype).all.OBmorletpeaks(t,1:(length(peaks)))=peaks;
        alldata.(tonetype).all.OBmorletpeaklocs(t,1:(length(locs)))=locs;
    end
    
    %save info about morlet wavelet/parameters
    alldata.info.LFPtime = LFPtime;
    alldata.info.wavelet.frequency = f;
    alldata.info.wavelet.ncycles = n;
    alldata.info.wavelet.s = s;
    alldata.info.wavelet.sformula = stext;
    alldata.info.wavelet.time = time;
    alldata.info.wavelet.sine_wave = sine_wave;
    alldata.info.wavelet.gaussian_win = gaussian_win;
    alldata.info.wavelet.wavelet = wavelet;
    alldata.info.wavelet.OBtraceinverted = 'yes';

    %save the file
    save(strcat(somefiles{h},savefiletag),'alldata','wholeshebang');
end %end file


%% Now files that don't require OBLFP inversion

otherfiles = [files(1:4);files(10:15);files(21:23)]; %rats 064, 071, and 074
for h=1:length(otherfiles)
    load(strcat(otherfiles{h},loadfiletag));
    if contains(somefiles{h},'toneon')
        tonetype = 'toneon';
        tasktypes = {'toneattn','odorattn'};
    else
        tonetype = 'toneoff';
        tasktypes = {'odoronly','toneattn','odorattn'};
    end
    
    % LFP data from OB for one trial
    %trialdata = alldata.(tonetype).all.OBLFP(10,1:6104);
    LFPtime = -3.6:1/alldata.info.LFPfs:2.4;

    f=8; %frequency of morlet wavelet
    n=3; %number of cycles
    s = n/(2*pi*f); 
    stext = '3/(2*pi*f)'; %for saving in data structure below
    time = -1:1/alldata.info.LFPfs:1;
    sine_wave = exp(1i*2*pi*f.*time);
    gaussian_win = exp(-time.^2./(2*s^2));
    wavelet = sine_wave .* gaussian_win;
    % half of the wavelet size, useful for chopping off edges after convolution.
    halfwaveletsize = ceil(length(wavelet)/2);
    
    for a=1:length(tasktypes)
        for b=1:length(correct)
            for c=1:length(trialtypes)
                currdata = (alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBLFP_theta(:,1:6104)); 
                if ~isempty(currdata)
                    %create matrix where the wavelet convolved data will go
                    alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorlet=zeros(length(currdata(:,1)),6104);
                    alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaks=NaN(length(currdata(:,1)),100);
                    alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaklocs=NaN(length(currdata(:,1)),100);
                    for t=1:length(currdata(:,1))
                        % convolve with data
                        % compute Gaussian
                        n_conv = length(wavelet) + length(currdata) - 1;
                        
                        fft_w = fft(wavelet,n_conv);
                        fft_e = fft(currdata(t,:),n_conv);
                        ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
                        wavelet_conv_data = real(ift(halfwaveletsize:end-halfwaveletsize+1));
                        
                        [peaks,locs] = findpeaks(wavelet_conv_data,LFPtime);
                        
                        alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorlet(t,:)=wavelet_conv_data;
                        alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaks(t,1:(length(peaks)))=peaks;
                        alldata.(tonetype).(tasktypes{a}).(correct{b}).(trialtypes{c}).OBmorletpeaklocs(t,1:(length(locs)))=locs;
                    end
                else
                end
            end
        end
    end %end tasktypes
    
    %repeat the same thing for the alldata struct
    %multiply currdata by -1, to invert the trace and detect
    %troughs as peaks
    currdata = alldata.(tonetype).all.OBLFP_theta(:,1:6104);
    %create matrix where the wavelet convolved data will go
    alldata.(tonetype).all.OBmorlet=zeros(length(currdata(:,1)),6104);
    alldata.(tonetype).all.OBmorletpeaks=NaN(length(currdata(:,1)),100);
    alldata.(tonetype).all.OBmorletpeaklocs=NaN(length(currdata(:,1)),100);
    for t=1:length(currdata(:,1))
        % convolve with data
        % compute Gaussian
        n_conv = length(wavelet) + length(currdata) - 1;

        fft_w = fft(wavelet,n_conv);
        fft_e = fft(currdata(t,:),n_conv);
        ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
        wavelet_conv_data = real(ift(halfwaveletsize:end-halfwaveletsize+1));

        [peaks,locs] = findpeaks(wavelet_conv_data,LFPtime);

        alldata.(tonetype).all.OBmorlet(t,:)=wavelet_conv_data;
        alldata.(tonetype).all.OBmorletpeaks(t,1:(length(peaks)))=peaks;
        alldata.(tonetype).all.OBmorletpeaklocs(t,1:(length(locs)))=locs;
    end
    
    %save info about morlet wavelet/parameters
    alldata.info.LFPtime = LFPtime;
    alldata.info.wavelet.frequency = f;
    alldata.info.wavelet.ncycles = n;
    alldata.info.wavelet.s = s;
    alldata.info.wavelet.sformula = stext;
    alldata.info.wavelet.time = time;
    alldata.info.wavelet.sine_wave = sine_wave;
    alldata.info.wavelet.gaussian_win = gaussian_win;
    alldata.info.wavelet.wavelet = wavelet;
    alldata.info.wavelet.OBtraceinverted = 'no';

    %save the file
    save(strcat(otherfiles{h},savefiletag),'alldata','wholeshebang');
end
