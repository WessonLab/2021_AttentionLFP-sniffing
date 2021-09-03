function varargout = trialGUI_wavelets_sniffing(varargin)
% TRIALGUI_WAVELETS_SNIFFING MATLAB code for trialGUI_wavelets_sniffing.fig
%      TRIALGUI_WAVELETS_SNIFFING, by itself, creates a new TRIALGUI_WAVELETS_SNIFFING or raises the existing
%      singleton*.
%
%      H = TRIALGUI_WAVELETS_SNIFFING returns the handle to a new TRIALGUI_WAVELETS_SNIFFING or the handle to
%      the existing singleton*.
%
%      TRIALGUI_WAVELETS_SNIFFING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIALGUI_WAVELETS_SNIFFING.M with the given input arguments.
%
%      TRIALGUI_WAVELETS_SNIFFING('Property','Value',...) creates a new TRIALGUI_WAVELETS_SNIFFING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trialGUI_wavelets_sniffing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trialGUI_wavelets_sniffing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trialGUI_wavelets_sniffing

% Last Modified by GUIDE v2.5 16-Apr-2021 15:40:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trialGUI_wavelets_sniffing_OpeningFcn, ...
                   'gui_OutputFcn',  @trialGUI_wavelets_sniffing_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before trialGUI_wavelets_sniffing is made visible.
function trialGUI_wavelets_sniffing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trialGUI_wavelets_sniffing (see VARARGIN)

% Choose default command line output for trialGUI_wavelets_sniffing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trialGUI_wavelets_sniffing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trialGUI_wavelets_sniffing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectfile.
function selectfile_Callback(hObject, eventdata, handles)
% hObject    handle to selectfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%select and load the file
filename = uigetfile('*.mat');
load (filename);
addpath(genpath('C:\Users\smelluser\Documents\MATLAB\chronux_2_12'));
plotregions = {'RESP'};
%display the file name
set (handles.filenamedisplay, 'String', erase(filename,'.mat'));

%determine the xvals for plotting the RESP and poke data
poketimes = (alldata.info.times(2,1):1/alldata.info.pokefs:alldata.info.times(2,2))+0.4;
RESPtimes = (alldata.info.times(2,1):1/alldata.info.RESPfs:alldata.info.times(2,2))+0.4;
if contains(filename,'toneon')
    tonetype = 'toneon';
else
    tonetype = 'toneoff';
end

%determine the number of trials and set the current trial to 1, display the
%current and total number of trials
[ntrials,~] = size(alldata.(tonetype).all.RESP);
if ~contains(filename,'rejected')
    currtrial = 1;
elseif contains(filename,'rejected')
    currtrial = alldata.info.GUIcurrtrial;
end
set (handles.maxtrialnumber,'String',num2str(ntrials));
set (handles.currenttrial,'String',num2str(currtrial));
set (handles.slider1, 'Min',1,'Max',ntrials,'Value',currtrial,'SliderStep',[1/(ntrials-1) 0.1]);

%plot the RESP data and poke data
for a=1:length(plotregions)
    axes(handles.(plotregions{a}));
    plot(RESPtimes,alldata.(tonetype).all.(plotregions{a})(currtrial,1:length(RESPtimes)));
    axis tight;
    xlabel 'Time from odor onset';
    ylabel '\muV';
    title (plotregions{a});
end
axes(handles.centerpokeplot);
plot(poketimes,alldata.(tonetype).all.Cpoke(currtrial,1:length(poketimes)));
axis tight;

%determine the default wavelet frequency/wavelet cycles
set (handles.waveletfrequency,'String',num2str(alldata.info.wavelet.frequency));
set (handles.ncycles,'String',num2str(alldata.info.wavelet.ncycles));

%plot the wavelet convolved signal + detected peaks with the RESP
axes(handles.RESPmorlet);
hold on;
if strcmp(alldata.info.wavelet.RESPtraceinverted,'no')
    morletRESP = plot(RESPtimes,alldata.(tonetype).all.RESP_filt(currtrial,1:length(RESPtimes)),'LineWidth',2,'hittest','off');
else
    morletRESP = plot(RESPtimes,(alldata.(tonetype).all.RESP_filt(currtrial,1:length(RESPtimes)).*(-1)),'LineWidth',2,'hittest','off');
end
morlettrace = plot(RESPtimes,alldata.(tonetype).all.RESPmorlet(currtrial,1:length(RESPtimes)),'LineWidth',2,'hittest','off');
peaks = alldata.(tonetype).all.RESPmorletpeaks(currtrial,:);
peaks = peaks(~isnan(peaks));
locs = alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,:);
locs = locs(~isnan(locs));
allpeaks = plot(locs,peaks,'ro','MarkerFaceColor','r','hittest','off');
axis tight;
xlim ([-2 2]);
xlabel 'Time from odor onset';
ylabel '\muV';
f= get(handles.waveletfrequency,'String');
title(strcat((num2str(f)),'Hz Morlet convolved, with peaks'));

if ~contains(filename,'rejected')
    %create a structure to store the info about rejected trials
    rejected.trials = zeros(ntrials,1);
    rejected.reason = cell(ntrials,1);
    rejected.notes = cell(ntrials,1);

    %create a structure to store info about wavelet params
    alldata.(tonetype).all.RESPmorletfreqs = ones(ntrials,1)*(alldata.info.wavelet.frequency);
    alldata.(tonetype).all.RESPmorletcycles = ones(ntrials,1)*(alldata.info.wavelet.ncycles);
elseif contains(filename,'rejected')
end

%set app data that will be used by other functions
setappdata (0,'alldata',alldata);
setappdata (0,'poketimes',poketimes);
setappdata (0,'RESPtimes',RESPtimes);
setappdata (0,'plotregions',plotregions);
setappdata (0,'rejected',rejected);
setappdata (0,'allpeaks',allpeaks);
setappdata (0,'morletRESP',morletRESP);
setappdata (0,'morlettrace',morlettrace);
setappdata (0,'tonetype',tonetype);
setappdata (0,'wholeshebang',wholeshebang);
setappdata (0,'filename',filename);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

currtrial = round(get(hObject,'Value'));

%import needed app data
alldata = getappdata(0,'alldata');
poketimes = getappdata(0,'poketimes');
RESPtimes = getappdata(0,'RESPtimes');
plotregions = getappdata(0,'plotregions');
rejected = getappdata(0,'rejected');
morletRESP = getappdata(0,'morletRESP');
morlettrace = getappdata(0,'morlettrace');
allpeaks = getappdata(0,'allpeaks');
tonetype = getappdata(0,'tonetype');
set (handles.currenttrial,'String',num2str(currtrial));
set (handles.waveletfrequency,'String',num2str(alldata.(tonetype).all.RESPmorletfreqs(currtrial)));
set (handles.ncycles,'String',num2str(alldata.(tonetype).all.RESPmorletcycles(currtrial)));

%set the keep/reject status and reasons for the trial
if rejected.trials(currtrial)==0
    set(handles.keeptrial,'Value',1);
    set(handles.artifact,'Value',0);
    set(handles.highpower,'Value',0);
    set(handles.other,'Value',0);
    set(handles.notes,'String','Notes');
else
    set (handles.rejecttrial,'Value',1);
    reason = rejected.reason(currtrial);
    nreasons = length(reason{1,1});
    for r=1:nreasons
        set(handles.(reason{1,1}{r}),'Value',1);
    end
    set(handles.notes,'String',rejected.notes(currtrial));
end

%plot the RESP data and poke data
for a=1:length(plotregions)
    axes(handles.(plotregions{a}));
    plot(RESPtimes,alldata.(tonetype).all.(plotregions{a})(round(currtrial),1:length(RESPtimes)));
    axis tight;
    xlabel 'Time from odor onset';
    ylabel '\muV';
    title (plotregions{a});
end
axes(handles.centerpokeplot);
plot(poketimes,alldata.(tonetype).all.Cpoke(round(currtrial),1:length(poketimes)));
axis tight;

%plot the wavelet convolved signal + detected peaks with the RESP
axes(handles.RESPmorlet);
delete(morletRESP);
delete(morlettrace);
delete(allpeaks);
%if the plot unfiltered box is checked, delete the unfiltered trace from
%the plot and uncheck the box
if get(handles.unfiltered,'Value')==1
    set(handles.unfiltered,'Value',0);
    unfilteredRESP=getappdata(0,'unfilteredRESP');
    axes(handles.RESPmorlet);
    delete(unfilteredRESP);
else
end
hold on;

if strcmp(alldata.info.wavelet.RESPtraceinverted,'no')
    morletRESP = plot(RESPtimes,alldata.(tonetype).all.RESP_filt(currtrial,1:length(RESPtimes)),'Color',[0, 0.4470, 0.7410],'LineWidth',2,'hittest','off');
else
    morletRESP = plot(RESPtimes,(alldata.(tonetype).all.RESP_filt(currtrial,1:length(RESPtimes)).*(-1)),'Color',[0, 0.4470, 0.7410],'LineWidth',2,'hittest','off');
end
morlettrace = plot(RESPtimes,alldata.(tonetype).all.RESPmorlet(currtrial,1:length(RESPtimes)),'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'hittest','off');
peaks = alldata.(tonetype).all.RESPmorletpeaks(currtrial,:);
peaks = peaks(~isnan(peaks));
locs = alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,:);
locs = locs(~isnan(locs));
allpeaks = plot(locs,peaks,'ro','MarkerFaceColor','r','hittest','off');
axis tight;
xlim ([-2 2]);
xlabel 'Time from odor onset';
ylabel '\muV';
f= get(handles.waveletfrequency,'String');
title(strcat((num2str(f)),'Hz Morlet convolved, with peaks'));

%set allpeaks to appdata, so it can be accessed by the buttondownfcn
setappdata (0,'allpeaks',allpeaks);
setappdata (0,'morletRESP',morletRESP);
setappdata (0,'morlettrace',morlettrace);
setappdata (0,'allpeaks',allpeaks);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on slider1 and none of its controls.
function slider1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in artifact.
function artifact_Callback(hObject, eventdata, handles)
% hObject    handle to artifact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of artifact


% --- Executes on button press in highpower.
function highpower_Callback(hObject, eventdata, handles)
% hObject    handle to highpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of highpower


% --- Executes on button press in other.
function other_Callback(hObject, eventdata, handles)
% hObject    handle to other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of other



function notes_Callback(hObject, eventdata, handles)
% hObject    handle to notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of notes as text
%        str2double(get(hObject,'String')) returns contents of notes as a double


% --- Executes during object creation, after setting all properties.
function notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rejecttrial.
function rejecttrial_Callback(hObject, eventdata, handles)
% hObject    handle to rejecttrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rejecttrial


% --- Executes on button press in settrialtype.
function settrialtype_Callback(hObject, eventdata, handles)
% hObject    handle to settrialtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%determine which radio button is selected
t = get(handles.keeprejectgroup,'SelectedObject');
x = get(t,'String');

%determine which checkboxes are checked and assign strings to one variable
allreasons = {'artifact','highpower','other'};
for a=1:length(allreasons)
    rejectreasons(a) = get(handles.(allreasons{a}),'Value');
end
rejectreasons = find(rejectreasons);
for d=1:length(rejectreasons)
    reasons{d} = allreasons{d};
end
    
%get any notes from the notes box
notes = get(handles.notes,'String');
if strcmp(notes,'Notes');
    notes = '';
end

%get info about the current trial
currtrial = round(get(handles.slider1,'Value'));

%get the rejected struct from appdata
rejected = getappdata(0,'rejected');

%set the rejected struct
if strcmp(x,'Keep')
    rejected.trials(currtrial,1)=0;
elseif strcmp(x,'Reject')
    rejected.trials(currtrial,1)=1;
    rejected.reason(currtrial,1)={reasons};
    rejected.notes(currtrial,1)={notes};
end
setappdata (0,'rejected',rejected);



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%add the rejected struct to the all data struct and save it
alldata = getappdata(0,'alldata');
rejected = getappdata(0,'rejected');
%wholeshebang isn't modified at all, but should remain in the saved file for future use
wholeshebang = getappdata(0,'wholeshebang'); 

currtrial = str2num(get(handles.currenttrial,'String'));
alldata.info.GUIcurrtrial = currtrial;

filename = getappdata(0,'filename');

if ~contains(filename,'rejected')
    filename = strcat(get(handles.filenamedisplay, 'String'),'_rejected.mat');
    save(filename,'alldata','rejected','wholeshebang');
elseif contains(filename,'rejected')
    filename = get(handles.filenamedisplay, 'String');
    save(filename,'alldata','rejected','wholeshebang');
end


% --- Executes on slider movement.
function morletslider_Callback(hObject, eventdata, handles)
% hObject    handle to morletslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
xval = get(hObject,'Value');
set(handles.RESPmorlet,'Xlim',[xval-1.5 xval+1.5]);

% --- Executes during object creation, after setting all properties.
function morletslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to morletslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function waveletfrequency_Callback(hObject, eventdata, handles)
% hObject    handle to waveletfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of waveletfrequency as text
%        str2double(get(hObject,'String')) returns contents of waveletfrequency as a double


% --- Executes during object creation, after setting all properties.
function waveletfrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to waveletfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setfrequency.
function setfrequency_Callback(hObject, eventdata, handles)
% hObject    handle to setfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = str2num(get(handles.waveletfrequency,'String'));
ncycles = str2num(get(handles.ncycles,'String'));

%get necessary data
alldata = getappdata(0,'alldata');
RESPtimes = getappdata(0,'RESPtimes');
tonetype = getappdata(0,'tonetype');
currtrial = str2num(get(handles.currenttrial,'String'));
if strcmp(alldata.info.wavelet.RESPtraceinverted,'no')
    currdata = alldata.(tonetype).all.RESP_filt(currtrial,1:3663);
elseif strcmp(alldata.info.wavelet.RESPtraceinverted,'yes')
    currdata = (alldata.(tonetype).all.RESP_filt(currtrial,1:3663)).*(-1);
end
morlettrace = getappdata(0,'morlettrace');
allpeaks = getappdata(0,'allpeaks');

%re-do the morlet wavelet with the new frequency
s = ncycles/(2*pi*f); %4.5 is 'n' which can be changed
time = -1:1/alldata.info.RESPfs:1;
sine_wave = exp(1i*2*pi*f.*time);
gaussian_win = exp(-time.^2./(2*s^2));
wavelet = sine_wave .* gaussian_win;
halfwaveletsize = ceil(length(wavelet)/2);

%convolve the trial data
n_conv = length(wavelet) + length(currdata) - 1;
                        
fft_w = fft(wavelet,n_conv);
fft_e = fft(currdata,n_conv);
ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
wavelet_conv_data = real(ift(halfwaveletsize:end-halfwaveletsize+1));

%detect the peaks
[peaks,locs] = findpeaks(wavelet_conv_data,RESPtimes);

%replot the data
axes(handles.RESPmorlet);
delete(morlettrace);
delete(allpeaks);
hold on;
morlettrace = plot(RESPtimes,wavelet_conv_data,'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'hittest','off');
allpeaks = plot(locs,peaks,'ro','MarkerFaceColor','r','hittest','off');
axis tight;
xlim ([-2 2]);
title(strcat((num2str(f)),'Hz Morlet convolved, with peaks'));

%copied the section below from the save new morlet callback, because I kept
%forgetting to click save, resulting in the peaks beign saved but not info
%about the frequency/convolved trace
alldata.(tonetype).all.RESPmorletfreqs(currtrial)= f;
alldata.(tonetype).all.RESPmorletcycles(currtrial)=ncycles;
alldata.(tonetype).all.RESPmorlet(currtrial,1:length(wavelet_conv_data))=wavelet_conv_data;
%first set rows to NaNs, in case there are fewer peaks in the new set
alldata.(tonetype).all.RESPmorletpeaks(currtrial,:) = NaN;
alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,:) = NaN;
%now insert the new peak/loc data
alldata.(tonetype).all.RESPmorletpeaks(currtrial,1:length(allpeaks.YData)) = allpeaks.YData;
alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,1:length(allpeaks.XData)) = allpeaks.XData;

%set the modified wavelet convolved trace, detected peaks/locs to app data,
%so that it can be accessed by the savenewmorlet callback
setappdata(0,'wavelet_conv_data',wavelet_conv_data);
setappdata(0,'morlettrace',morlettrace);
setappdata(0,'allpeaks',allpeaks);
setappdata(0,'alldata',alldata);


% --- Executes on button press in savenewmorlet.
function savenewmorlet_Callback(hObject, eventdata, handles)
% hObject    handle to savenewmorlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currtrial = str2num(get(handles.currenttrial,'String'));
alldata = getappdata(0,'alldata');
f = str2num(get(handles.waveletfrequency,'String'));
ncycles = str2num(get(handles.ncycles,'String'));
wavelet_conv_data = getappdata(0,'wavelet_conv_data');
allpeaks = getappdata(0,'allpeaks');
tonetype = getappdata(0,'tonetype');

%assign the modified morlet frequency/convolved trace/detected peaks to the RESPmorletfreqs matrix
alldata.(tonetype).all.RESPmorletfreqs(currtrial)= f;
alldata.(tonetype).all.RESPmorletcycles(currtrial)=ncycles;
alldata.(tonetype).all.RESPmorlet(currtrial,1:length(wavelet_conv_data))=wavelet_conv_data;
%first set rows to NaNs, in case there are fewer peaks in the new set
alldata.(tonetype).all.RESPmorletpeaks(currtrial,:) = NaN;
alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,:) = NaN;
%now insert the new peak/loc data
alldata.(tonetype).all.RESPmorletpeaks(currtrial,1:length(allpeaks.YData)) = allpeaks.YData;
alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,1:length(allpeaks.XData)) = allpeaks.XData;

setappdata(0,'alldata',alldata);



function ncycles_Callback(hObject, eventdata, handles)
% hObject    handle to ncycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncycles as text
%        str2double(get(hObject,'String')) returns contents of ncycles as a double


% --- Executes during object creation, after setting all properties.
function ncycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function RESPmorlet_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to RESPmorlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mousePos = get(hObject,'CurrentPoint');
mousePos = mousePos(1,1:2);

%get the peaks/locs data for this trial
allpeaks = getappdata(0,'allpeaks');
peakdata = [allpeaks.XData',allpeaks.YData'];

xdifference = peakdata(:,1)-mousePos(1);
[~,nearestidx] = min(abs(xdifference));
bothpts = [peakdata(nearestidx,:);mousePos];
distance = pdist(bothpts,'euclidean');
if distance <= 600
    %edit and replot the peak data
    peakdata(nearestidx,:)=[];
    axes(handles.RESPmorlet);
    delete(allpeaks);
    allpeaks = plot(peakdata(:,1),peakdata(:,2),'ro','MarkerFaceColor','r','hittest','off');
    
    %import and edit alldata struct
    alldata=getappdata(0,'alldata');
    currtrial = str2num(get(handles.currenttrial,'String'));
    tonetype = getappdata(0,'tonetype');
    alldata.(tonetype).all.RESPmorletpeaks(currtrial,:) = NaN;
    alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,:) = NaN;
    %now insert the new peak/loc data
    alldata.(tonetype).all.RESPmorletpeaks(currtrial,1:length(allpeaks.YData)) = allpeaks.YData;
    alldata.(tonetype).all.RESPmorletpeaklocs(currtrial,1:length(allpeaks.XData)) = allpeaks.XData;
    
    %save new allpeaks/alldata to appdata
    setappdata(0,'allpeaks',allpeaks);
    setappdata(0,'alldata',alldata);
else
end


% --- Executes on button press in unfiltered.
function unfiltered_Callback(hObject, eventdata, handles)
% hObject    handle to unfiltered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of unfiltered

alldata = getappdata(0,'alldata');
RESPtimes = getappdata(0,'RESPtimes');
currtrial = str2num(get(handles.currenttrial,'String'));
tonetype = getappdata(0,'tonetype');

if get(hObject,'Value')==1
    %plot the unfiltered OB trace
    axes(handles.RESPmorlet);
    hold on;
    if strcmp(alldata.info.wavelet.RESPtraceinverted,'no')
        unfilteredRESP = plot(RESPtimes,alldata.(tonetype).all.RESP(currtrial,1:length(RESPtimes)),'Color',[0.4660, 0.6740, 0.1880],'hittest','off');
    else
        unfilteredRESP = plot(RESPtimes,(alldata.(tonetype).all.RESP(currtrial,1:length(RESPtimes)).*(-1)),'Color',[0.4660, 0.6740, 0.1880],'hittest','off');
    end
    axis tight;
    xlim ([-2 2]);
    setappdata (0,'unfilteredRESP',unfilteredRESP);
elseif get(hObject,'Value')==0
    unfilteredRESP = getappdata(0,'unfilteredRESP');
    axes(handles.RESPmorlet);
    delete (unfilteredRESP);
end

% 
% % --- Executes on mouse press over axes background.
% function RESPmorlet_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to RESPmorlet (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
