%----Luciferase Assay Analyzer: Graphical User Interface for analyzing luciferase assay data----

    %software: MikroWin version 5.19
    %machine:Berthold Centro XS3 LB960

    %data collection: settings -> measurement -> repeated: enter total time (total duration assay),
        %counting time (how long each measurement takes), cycle time in seconds (after how much time the measurement is repeated) 

    %setting up export: installation -> driver -> export -> RawData Export Driver (double click)
        % -> data layout (by well), kinetic layout (Time/Position), time format
        % (ss), operation mode (no box ticked), Statistics export (general
        % statistics), format (text file), directory (C:\MikroWinData\Transfer)

    % Luciferase data export: file -> export -> Active Export Driver (RawData Export Driver)

    % Luciferase data format: *.txt

    % For each experiment a plate-design (*_pd.xlsx) file should be created that
        % has the same filename as the *.txt file, supplemented with *_pd.xlsx
        % The plate design indicates for each well the experimental condition
        % The format of the *_pd.xlsx file is the following:
            % A1: empty cell
            % A2, A3, A4, ... : A, B, C, ...(alphabet of nr of rows of plate)
            % B1, C1, D1, ... : 1, 2, 3, ...(nr of columns of plate)
            % fill from B2 with the name of the condition in each well

    % How to use the GUI:    
        % open GUI by clicking run in MATLAB
        % load data by clicking 'load data' in GUI and selecting *.txt file
        % 

function varargout = gui_select_molts(varargin)
% GUI_SELECT_MOLTS MATLAB code for gui_select_molts.fig
%      GUI_SELECT_MOLTS, by itself, creates a new GUI_SELECT_MOLTS or raises the existing
%      singleton*.
%
%      H = GUI_SELECT_MOLTS returns the handle to a new GUI_SELECT_MOLTS or the handle to
%      the existing singleton*.
%
%      GUI_SELECT_MOLTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SELECT_MOLTS.M with the given input arguments.
%
%      GUI_SELECT_MOLTS('Property','Value',...) creates a new GUI_SELECT_MOLTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_select_molts_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_select_molts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_select_molts

% Last Modified by GUIDE v2.5 10-Feb-2020 15:52:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_select_molts_OpeningFcn, ...
    'gui_OutputFcn',  @gui_select_molts_OutputFcn, ...
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


% --- Executes just before gui_select_molts is made visible.
function gui_select_molts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_select_molts (see VARARGIN)

% Choose default command line output for gui_select_molts
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = gui_select_molts_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning off
%Load luciferase data
[filename,filepath]=uigetfile('*.*','Select Luciferase Data ');

%remove extensions from file name
idx = strfind(filename,'.');
idx=idx(end);%in case the filename contains more then 1 .
ext=filename(idx+1:end);
filename=filename(1:idx-1);

%put file name in handles
handles.filename=filename;handles.filepath=filepath;

% there are two possible filetypes:
% *.txt is the raw data as exported by MikroWin2010 software using settings
    % described above
% *.mat if the dataset has already been opened in the GUI

if exist([filepath,filename,'.mat'],'file')
    handles=load_data(handles);
else
    if strcmpi(ext,'xlsx')
        [X,samples] = xlsread([filepath,filename]);
        handles.sampleRate = 1/12;%set manually in hours
        first_sample = 2;
        handles.samples = samples(first_sample:end);
        handles.X = X(:,first_sample:end);
    elseif strcmpi(ext,'txt')
        formatSpec = '%s%[^\n\r]';
        fileID = fopen([filepath,filename,'.',ext],'r');
        delimiter = '';
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
        dataArray =dataArray{1:end-1};
        [~, ~, PD] = xlsread([filepath,filename,'_pd.xlsx'],'Sheet1');
        
        %find line where data starts, this is specific for the *.txt file from
            %the Mikrowin2010 software. The lines of code might break if the format of the
            %exported txt file changes.
        i=1;
        while isempty(strfind(dataArray{i},'Time/Pos'))
            i=i+1;
        end
        c=1;
        for j=i+1:length(dataArray)
            C = textscan(dataArray{j},'%s');
            tmp=cell2mat(C{1}(1));
            col=find(strcmp(tmp(1),PD(:,1)));
            row=find(str2double(tmp(2:end))==cell2mat(PD(1,:)));
            
            if ~isnan(PD{col,row})
                handles.samples(c)=PD(col,row);
                handles.X(:,c)=cellfun(@str2double,C{1}(2:end));
                c=c+1;
            end
        end
       
        %find sample rate. The lines of code might break if the format of the
            %exported txt file changes.
        i=1;
        while isempty(strfind(dataArray{i},'Cycle Time [s]'))
            i=i+1;
        end
        ind1 = strfind(dataArray{i},'Cycle Time [s]')+length('Cycle Time [s]');
        ind2 = strfind(dataArray{i},'Measurement')-1;
        tmp = dataArray{i}(ind1:ind2);
        tmp(strfind(tmp,' '))=[];tmp(strfind(tmp,'.'))=[];
        handles.sampleRate = str2double(tmp)/36000;
    end
    
    handles.Xc = [];
    
        
    %set all molt-larval labels to 0
    handles.y = zeros(size(handles.X));
    for i=1:size(handles.X,2)
        %detect hatch: the first time point (starting from TP4 -> to avoid
            %edge effects) that exceeds the mean+5*stdev of the first 20 TPs
        try;handles.Hatch(i)=find(1:length(handles.X(:,i))>3 & handles.X(:,i)'>(mean(handles.X(1:20,i))+5*std(handles.X(1:20,i))),1,'first')-1;end
        %trend correct data: According to Olmedo et al., Genetics, 2015 (Methods and Fig S2C,D)
        handles.Xc(:,i) = calc_trend_corrected_data(handles.X(:,i));
    end
    
    handles.Molt = ones(4,2,length(handles.samples));
    handles.annotated = logical(zeros(size(handles.samples)));% annotated is a logical, indicating whether a sample has been annotated (1) or not (0)
    handles.valid = logical(ones(size(handles.samples)));% valid is a logical, indicating whether a sample is valid (1) or not (0) 
end

%plot luminescence data and set values for dropdown menus
handles.uniq_samples = unique(handles.samples);
    set(handles.Strain1Selection,'String',handles.uniq_samples);
    set(handles.Strain2Selection,'String',handles.uniq_samples);
    set(handles.Strain1Selection,'Value',1);
    set(handles.Strain2Selection,'Value',2);
    set(handles.strain1Lumi,'String',handles.uniq_samples);
    set(handles.strain2Lumi,'String',handles.uniq_samples);
    set(handles.strain1Lumi,'Value',1);
    set(handles.strain2Lumi,'Value',2);

%plot data
handles.current_sample=1;
handles.ymaxlim = 1.1*max(handles.X(:,1));
plot_data(handles)

%update and save handles
guidata(hObject, handles);
save_data(handles)

% --- Executes on button press in PreviousSample.
function PreviousSample_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.current_sample > 1
    handles.current_sample = handles.current_sample - 1;
    handles.ValidSample.Value = handles.valid(handles.current_sample);
    %plot data
    plot_data(handles)
end

%set M1-M4 to valid for previous sample
handles.M1valid.Value=1;
handles.M2valid.Value=1;
handles.M3valid.Value=1;
handles.M4valid.Value=1;

%update and save handles
guidata(hObject, handles);
save_data(handles)

% --- Executes on button press in NextSample.
function NextSample_Callback(hObject, eventdata, handles)
% hObject    handle to NextSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.current_sample < length(handles.samples)
    handles.current_sample = handles.current_sample +1;
    handles.ValidSample.Value = handles.valid(handles.current_sample);
    %plot data
    plot_data(handles)
    c = handles.current_sample;
    try
        if handles.valid(c)==1 && ~handles.annotated(c)
            [handles.y(:,c), handles.Molt(:,:,c)]= detect_molts(handles.Xc(:,c),handles);
        end
    end
     
    %plot data
    plot_data(handles)
end

%set M1-M4 to valid for next sample
handles.M1valid.Value=1;
handles.M2valid.Value=1;
handles.M3valid.Value=1;
handles.M4valid.Value=1;

%update and save handles
guidata(hObject, handles);
save_data(handles)

% --- Executes on button press in AnnotateMolts.
function AnnotateMolts_Callback(hObject, eventdata, handles,m)
% hObject    handle to AnnotateMolts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%select start and end of the 4 molts
c = handles.current_sample;

if handles.valid(c)
    %select molts which where not valid
        [x, ~,button] = ginput(2);
        x=x/handles.sampleRate;
        if ~any(button==3)
            ind = ceil(x(1)):floor(x(2));
            handles.Molt(m,:,c)=[ind(1) ind(end)];
            handles.y(ind,c)= 1;
        end
        plot_data(handles)
    idx = find(handles.Molt(:,1,c)~=0);
        handles.Molt(1:length(idx),:,c) = sort(handles.Molt(idx,:,c));
end

handles.annotated(c)=true;

%update and save handles
save_data(handles)
guidata(hObject, handles);

% --- Executes on button press in AnnotateHatch.
function AnnotateHatch_Callback(hObject, eventdata, handles)
% hObject    handle to AnnotateHatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%manual annotation of hatch
c = handles.current_sample;

%select hatch
[x, ~,button] = ginput(1);
if any(button==3)
else
    handles.Hatch(c)=round(x*(1/handles.sampleRate));
    ind=1:handles.Hatch(c);
    plot_data(handles)
end

guidata(hObject, handles);
save_data(handles)

% --- Executes on button press in ValidSample.
function ValidSample_Callback(hObject, eventdata, handles)
% hObject    handle to ValidSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%indicate whether sample is valid or not and should be included in
    %boxplots, heatmaps and statistics or nor
c = handles.current_sample;
handles.valid(c) = handles.ValidSample.Value;
if handles.ValidSample.Value==0
    handles.y(:,c)=0;
    handles.Hatch(c)=0;
else
    try;handles.Hatch(c)=find(1:length(handles.X(:,c))>3 & handles.X(:,c)'>(mean(handles.X(1:20,c))+5*std(handles.X(1:20,c))),1,'first')-1;end%end
end
plot_data(handles)

%update handles
guidata(hObject, handles);
save_data(handles)


function save_data(handles)
X=handles.X;
Xc=handles.Xc;
y=handles.y;
valid=handles.valid;
current_sample=handles.current_sample;
samples=handles.samples;
Molt=handles.Molt;
Hatch=handles.Hatch;
sampleRate=handles.sampleRate;
annotated=handles.annotated;
save([handles.filepath,handles.filename,'.mat'],'X','Xc','y','valid','current_sample','samples','Molt','Hatch','sampleRate','annotated')

function handles=load_data(handles)
load([handles.filepath,handles.filename,'.mat'])
handles.X=X;
handles.Xc=Xc;
handles.y=y;
handles.valid=valid;
handles.current_sample=current_sample;
handles.samples=samples;
handles.Molt=Molt;
handles.Hatch=Hatch;
handles.sampleRate=sampleRate;
handles.annotated=annotated;


% --- Executes on button press in ExportToExcel.
function ExportToExcel_Callback(hObject, eventdata, handles)
% hObject    handle to ExportToExcel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %generates csv file with tab-delimited data of Hatch, MoltN_entry and
    %MoltN_exit of valid samples

j=1;
results={};
for i=1:length(handles.samples)
    if handles.valid(i)==1
        results{j,1}=handles.samples{i};
        results{j,2}=handles.Hatch(i)*handles.sampleRate;
        results{j,3}=handles.Molt(1,1,i)*handles.sampleRate;
        results{j,4}=handles.Molt(1,2,i)*handles.sampleRate;
        results{j,5}=handles.Molt(2,1,i)*handles.sampleRate;
        results{j,6}=handles.Molt(2,2,i)*handles.sampleRate;
        results{j,7}=handles.Molt(3,1,i)*handles.sampleRate;
        results{j,8}=handles.Molt(3,2,i)*handles.sampleRate;
        results{j,9}=handles.Molt(4,1,i)*handles.sampleRate;
        results{j,10}=handles.Molt(4,2,i)*handles.sampleRate;
        j=j+1;
    end
end
fileID = fopen([handles.filepath,handles.filename(1:end),'_molt_detection_results.csv'],'w');
formatSpec = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n';
formatSpec2 = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
[nrows,ncols] = size(results);
labels = {'condition','hatch','M1_entry','M1_exit','M2_entry','M2_exit','M3_entry','M3_exit','M4_entry','M4_exit'};
table = [labels;results];
[nrowTable,ncolTable] = size(table);

fprintf(fileID,formatSpec2,table{1,:})
for row = 2:nrowTable
    fprintf(fileID,formatSpec,table{row,:});
end

fclose(fileID);




% --- Executes on button press in MakeHeatmap.
function MakeHeatmap_Callback(hObject, eventdata, handles)
% hObject    handle to MakeHeatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%size of heatmap depends on max duration
max_duration=max(squeeze(handles.Molt(4,2,:))-handles.Hatch')+60;
if max(max_duration+handles.Hatch) >  size(handles.X,1)
    max_duration =  min(size(handles.X,1)-handles.Hatch');
end

%pop-up on how data should be plotted and sorted in the heatmap
    %include hatch
    j=1;tmp_data_bin=[];tmp_data=[];
    choice_include_hatch = questdlg('Include Hatch?','Heatmaps','Yes','No','No');
    switch choice_include_hatch
        case 'No'
    for i=1:length(handles.samples)
        tmp_data_bin(j,:)=handles.y(handles.Hatch(i)+1:handles.Hatch(i)+max_duration,i)';
        tmp_data(j,:)=handles.Xc(handles.Hatch(i)+1:handles.Hatch(i)+max_duration,i)';
        tmp_data_raw(j,:)=handles.X(handles.Hatch(i)+1:handles.Hatch(i)+max_duration,i)';
        j=j+1; 
    end
    t=linspace(1,max_duration*handles.sampleRate);
        case 'Yes'
         for i=1:length(handles.samples)
        tmp_data_bin(j,:)=handles.y(:,i)';
        tmp_data(j,:)=handles.Xc(:,i)';
        tmp_data_raw(j,:)=handles.X(:,i)';
        j=j+1; 
        end
        t = linspace(1,size(handles.X,1)*handles.sampleRate);
    end

tmp=unique(handles.samples);

    %each strain in seperate heatmap or all strains together in one heatmap
        %and sort on which molt
    choice_together_or_seperately = questdlg('What type of heatmap do you want to plot?','Heatmaps','Strains separately','All strains together','Strains separately');
    molt2sortOn = questdlg('Sort on wich molt?','Select molt','1','2','3','1'); 
    molt2sortOn = str2double(molt2sortOn);
    switch choice_together_or_seperately
        
        %strains seperately
         case 'Strains separately'
    for i=1:length(tmp)
        ind=find(handles.valid==1 & strcmp(handles.samples,tmp{i}));
        if ~isempty(ind)
            switch choice_include_hatch
                case 'No'
                    [~,ind2]=sort(squeeze(handles.Molt(1,1,ind))-handles.Hatch(ind)');
                case 'Yes'
                    [~,ind2]=sort(squeeze(handles.Molt(molt2sortOn,1,ind)));
            end
        %make heatmaps
    figure;imagesc(t,1:j-1,tmp_data_bin(ind(ind2),:));xlabel('Time (h)');set(gca,'YTickLabel','');title(tmp{i});pause(0.1)
    figure;imagesc(t,1:j,log(tmp_data(ind(ind2),:)));xlabel('Time (h)');set(gca,'YTickLabel','');title(tmp{i});colorbar;colormap(gray);pause(0.1)
        end
    end
    
        %all strains together
        case 'All strains together'
    ind=find(handles.valid==1);
    switch choice_include_hatch
        case 'No'
            [~,ind2]=sort(squeeze(handles.Molt(1,1,ind))-handles.Hatch(ind)');
        case 'Yes'
            [~,ind2]=sort(squeeze(handles.Molt(molt2sortOn,1,ind)));
    end

        %make heatmaps
    figure;imagesc(t,1:j-1,tmp_data_bin(ind(ind2),:));xlabel('Time (h)');set(gca,'YTickLabel','');pause(0.1)
    figure;imagesc(t,1:j,log(tmp_data(ind(ind2),:)));xlabel('Time (h)');set(gca,'YTickLabel','');colorbar;colormap(gray);pause(0.1)
 end



function Xc = calc_trend_corrected_data(X)
%trend correct the data according to Olmedo et al 2015
windowSize=72;  %sliding window of 72
Xma = filter(ones(1,windowSize)/windowSize,1,X);
windowSize=1;
Xf = filter(ones(1,windowSize)/windowSize,1,X);
Xc = Xf./Xma;

function [y, Molt] = detect_molts(X,handles)
y=zeros(size(X));
dy=([0;diff(X<0.75)]);
ind=find(dy~=0);
valid=ones(size(ind));
for i=1:length(ind)
    N1=1/handles.sampleRate;
    if(dy(ind(i)))==1 %potential onset
        if ind(i)+N1 >length(dy)
            valid(i)=0;
        elseif any(dy(ind(i)+1:ind(i)+N1)~=0)
            valid(i)=0;
        end
    else %potential end
        if ind(i)-N1 <1 || ind(i)+N1 >length(dy)%offsett should be preceded by 1 hour, extra:
            valid(i)=0;
        elseif any(dy(ind(i)-N1:ind(i)-1)~=0)
            valid(i)=0;
        end
    end
end
ind(valid==0)=[];
if ~isempty(ind)
if dy(ind(1))==-1;ind(1)=[];end
if dy(ind(end))==1;ind(end)=[];end
onset=ind(dy(ind)==1);
offset=ind(dy(ind)==-1)-1;
durations=offset-onset;
if length(onset)>4
    [~,tmp]=sort(durations,'descend');
    onset=onset(tmp(1:4));
    offset=offset(tmp(1:4));
end
Molt=zeros(4,2);
for i=1:length(onset)
    y(onset(i):offset(i))=1;
    Molt(i,1:2)=[onset(i) offset(i)];
end
if any(Molt(:,1)==0)
        ind_nz = find(Molt(:,1)~=0); % ind_nz == nonzero molts
        ind_z = find(Molt(:,1)==0);
        [~,idx]= sort(Molt(ind_nz,1));
        Molt = [Molt(ind_nz(idx),:); Molt(ind_z,:)];
    else
        [~,idx]= sort(Molt(:,1));
        Molt = Molt(idx,:); 
    end
else
   Molt=zeros(4,2); 
end


% --- Executes on button press in DetectMolts.
function DetectMolts_Callback(hObject, eventdata, handles)
% hObject    handle to DetectMolts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c = handles.current_sample;
if handles.valid(c)==1
    [handles.y(:,c) handles.Molt(:,:,c)]= detect_molts(handles.Xc(:,c),handles);
    %if 4 molts are detected, sort them. if less then 4 molts are found set
        %last molts to zero
    handles.annotated(c)=false;
end
plot_data(handles)
%set M1-M4 to valid for next sample
handles.M1valid.Value=1;
handles.M2valid.Value=1;
handles.M3valid.Value=1;
handles.M4valid.Value=1;
%save data
guidata(hObject, handles);
save_data(handles)


% --- Executes on button press in M1valid.
function M1valid_Callback(hObject, eventdata, handles)
c = handles.current_sample;

if ~any(handles.Molt(:,1,c)==0)
    handles.Molt(1,:,c)=0;
else
    idx = find(handles.Molt(:,1,c)~=0);
    handles.Molt(idx(1),:,c)=0;
end

handles.y(cumsum(abs(([0; diff(handles.y(:,c))])))==1,c)=0;
guidata(hObject, handles);
save_data(handles)
plot_data(handles)
AnnotateMolts_Callback(hObject, eventdata, handles,1);
handles.M1valid.Value=1;


% --- Executes on button press in M2valid.
function M2valid_Callback(hObject, eventdata, handles)
c = handles.current_sample;
if ~any(handles.Molt(:,1,c)==0)
    handles.Molt(2,:,c)=0;
else
    idx = find(handles.Molt(:,1,c)~=0);
    if length(idx)>1
    handles.Molt(idx(2),:,c)=0;
    end
end

handles.y(cumsum(abs(([0; diff(handles.y(:,c))])))==3,c)=0;
guidata(hObject, handles);
save_data(handles)
plot_data(handles)
AnnotateMolts_Callback(hObject, eventdata, handles,2);
handles.M2valid.Value=1;


% --- Executes on button press in M3valid.
function M3valid_Callback(hObject, eventdata, handles)
c = handles.current_sample;
if ~any(handles.Molt(:,1,c)==0)
    handles.Molt(3,:,c)=0;
else
    idx = find(handles.Molt(:,1,c)~=0);
    if length(idx)>2
    handles.Molt(idx(3),:,c)=0;
    end
end

handles.y(cumsum(abs(([0; diff(handles.y(:,c))])))==5,c)=0;
guidata(hObject, handles);
save_data(handles)
plot_data(handles)
AnnotateMolts_Callback(hObject, eventdata, handles,3);
handles.M3valid.Value=1;


% --- Executes on button press in M4valid.
function M4valid_Callback(hObject, eventdata, handles)
c = handles.current_sample;
if ~any(handles.Molt(:,1,c)==0)
    handles.Molt(4,:,c)=0;
else
    idx = find(handles.Molt(:,1,c)~=0);
    if length(idx)>3
    handles.Molt(idx(4),:,c)=0;
    end
end
handles.y(cumsum(abs(([0; diff(handles.y(:,c))])))==7,c)=0;
guidata(hObject, handles);
save_data(handles)
plot_data(handles)
AnnotateMolts_Callback(hObject, eventdata, handles,4);
handles.M4valid.Value=1;


function plot_data(handles)
c = handles.current_sample;
%plot data
handles.axes1;
t=[1:length(handles.X(:,c))]*handles.sampleRate;
hold off;plot(t,handles.X(:,c));title([handles.samples{c},'  ',num2str(c),'/',num2str(size(handles.y,2))])
ind = find(handles.y(:,c)==1);
hold on;plot(t(ind),handles.X(ind,c),'g.');
handles.current_sample = c;
ind=1:handles.Hatch(c);
hold on; plot(t(ind),handles.X(ind,c),'r.');
%axis([1 t(end) 0 3000]);
%axis([1 t(end) 0 1.1*max(handles.X(:,c))]);
axis([1 t(end) 0 handles.ymaxlim]);


function GoTo_Callback(hObject, eventdata, handles)
% hObject    handle to GoTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%enter a sample number and go to this sample immediately
if ~isnan(str2double(eventdata.Source.String))
    handles.current_sample=str2double(eventdata.Source.String);
    handles.ValidSample.Value = handles.valid(handles.current_sample);
    plot_data(handles)
    guidata(hObject, handles);
    save_data(handles)
else
    disp([(eventdata.Source.String),' is not a number!'])
    eventdata.Source.String=num2str(handles.current_sample);
end
% Hints: get(hObject,'String') returns contents of GoTo as text
%        str2double(get(hObject,'String')) returns contents of GoTo as a double


% --- Executes during object creation, after setting all properties.
function GoTo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GoTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%to use the rightarrow and leftarrow on the keyboard to go to next and
    %previous sample respectively
if strcmp(eventdata.Key,'rightarrow')
    NextSample_Callback(hObject, eventdata, handles)
elseif strcmp(eventdata.Key,'leftarrow')
    PreviousSample_Callback(hObject, eventdata, handles)
elseif strcmp(eventdata.Key,'escape')
    
end


% --- Executes on button press in MakeBoxplots.
function MakeBoxplots_Callback(hObject, eventdata, handles)
% hObject    handle to MakeBoxplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% let pop-up appear to choose strains to be plotted in the boxplot
    uniq_samples_stable = unique(handles.samples, 'stable'); %order is as in _pd (handles.uniq_samples is alphabetically/numerically ordered)
    Prompt{1}= 'Strains  to be plotted: ';
    Formats(1).type = 'list';
    Formats(1).style = 'listbox';
    Formats(1).items = uniq_samples_stable;
    Formats(1).format = 'text';
    Formats(1).limits = [0 length(Formats(1).items)];
    Formats(1).size = [200 length(Formats(1).items)*20];
    [cfg,Cancelled] = inputsdlg(Prompt,'Select Strains',Formats); %cfg=number indicating selected strains


% let popup appear to adapt the order of the strains to be plotted
    %enter a 1 for the strain to be plotted first in boxplot, a 2 for the
    %strain to be plotted second, a 3 for the third, and so on
    
    strains_name = uniq_samples_stable(cfg{1}); %strains to be plotted

    for i=1:length(strains_name)
        Prompt_order{i} = strains_name{i};
        Formats_order(i).type = 'list';
        Formats_order(i).style = 'popupmenu';
        Formats_order(i).items = num2cell(1:length(strains_name));
        Formats_order(i).format = 'text';
        Formats_order(i).limits = [0 length(1:length(strains_name))];
        Formats_order(i).size = [200 length(1:length(strains_name))*20];
    end
        strains_orderID = inputsdlg(Prompt_order','Strain order',Formats_order)'; %indicates for each strain on which spot you want it to be: 1 for first, 2 for second
        
        %transform spot into column position
        for i = 1:length(strains_orderID)
            strains_orderPos(strains_orderID{i}) = i;
        end 
        %order strains
        strains_name_ordered = strains_name(strains_orderPos);
    
%sort handles.samples by strain order
       samples_orderID = (zeros(size(handles.samples)));
    for i=1:length(strains_orderID)
       samples_orderID(strcmp(handles.samples,strains_name_ordered(i)))=i; %couple the sample idx to the order
    end
        [samples_orderID_sorted,idx_order] = sort(samples_orderID);
        handles.samples_ordered = handles.samples(idx_order);
        handles.valid_ordered = handles.valid(idx_order);

%sort handles.Molt and handles.Hatch
    handles.Molt_ordered = handles.Molt(:,:,idx_order);
    handles.Hatch_ordered = handles.Hatch(:,idx_order);
    
%determine the Molt, Larval stage and Intermolt durations
    M=squeeze(handles.Molt_ordered(:,2,:)-handles.Molt_ordered(:,1,:))';
    L(:,1)=squeeze(handles.Molt_ordered(1,2,:))-handles.Hatch_ordered';
    L(:,2:4)=squeeze(handles.Molt_ordered(2:4,2,:)-handles.Molt_ordered(1:3,2,:))';
    I(:,1)=squeeze(handles.Molt_ordered(1,1,:))-handles.Hatch_ordered';
    I(:,2:4)=squeeze(handles.Molt_ordered(2:4,1,:)-handles.Molt_ordered(1:3,2,:))';
     
%define boxplot positions, width of boxplot is 0.5
    % and distance between the center of 2 groups of boxplots is 0.25 
    box_positions = [];
    width = length(cfg{1})+0.5;
    is_even = rem(length(cfg{1}),2)==0;
    for i=1:4;
        if is_even
            shift = (i-1)*width + 0.5;
            box_positions = cat(2,box_positions, (1:length(cfg{1}))-0.5+shift);
        else
            shift = (i-1)*width + 1;
            box_positions = cat(2,box_positions, (1:length(cfg{1}))-0.5+shift); 
        end
    end


%Boxplots according to defined order
   
    %use selected strains and valid sampled
    ind = logical(zeros(size(handles.samples_ordered))); 
    for i=1:length(strains_name_ordered)
       ind(strcmp(handles.samples_ordered,strains_name_ordered(i)) & handles.valid_ordered)=1;
    end

    %generate x-labels for te boxes
    labels=repmat({'',' ','  ','   '},sum(ind),1);
    labels=labels(:)';
    final_labels = strcat(repmat(handles.samples_ordered(ind),1,4),labels);
    clear labels
    if is_even
        tick_positions = box_positions(ceil(length(cfg{1})/2)+[0:3]*length(cfg{1}))+0.5;
    else
        tick_positions = box_positions(ceil(length(cfg{1})/2)+[0:3]*length(cfg{1}));
    end
    
    %generate boxplots
    tmp_M = M(ind,:)*handles.sampleRate;
    figure;boxplot(tmp_M(:),final_labels,'ColorGroup',repmat(handles.samples_ordered(ind),1,4),...
       'positions',box_positions,'LabelOrientation','inline');title('Molts');ylabel('Time (hrs)');axis([0.5 box_positions(end)+0.5 0 1.25*max(tmp_M(:))])
    hold on;for i=1:4;text(tick_positions(i),0.15*min(tmp_M(:)),['M',num2str(i)],'HorizontalAlignment','Center');end
   
    tmp_L = L(ind,:)*handles.sampleRate;
    figure;boxplot(tmp_L(:),final_labels,'ColorGroup',repmat(handles.samples_ordered(ind),1,4),...
       'positions',box_positions,'LabelOrientation','inline');title('Larval stages');ylabel('Time (hrs)');axis([0.5 box_positions(end)+0.5 0 1.25*max(tmp_L(:))]);
    hold on;for i=1:4;text(tick_positions(i),0.25*min(tmp_L(:)),['L',num2str(i)],'HorizontalAlignment','Center');end
    
    tmp_I = I(ind,:)*handles.sampleRate;
    figure;boxplot(tmp_I(:),final_labels,'ColorGroup',repmat(handles.samples_ordered(ind),1,4),...
       'positions',box_positions,'LabelOrientation','inline');title('Intermolts');ylabel('Time (hrs)');axis([0.5 box_positions(end)+0.5 0 1.25*max(tmp_I(:))]);
    hold on;for i=1:4;text(tick_positions(i),0.25*min(tmp_I(:)),['I',num2str(i)],'HorizontalAlignment','Center');end
  
 
% --- Executes on button press in Statistics.
function Statistics_Callback(hObject, eventdata, handles)
% hObject    handle to Statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %determine the Molt, Larval stage and Intermolt durations
    M=squeeze(handles.Molt(:,2,:)-handles.Molt(:,1,:))';
    L(:,1)=squeeze(handles.Molt(1,2,:))-handles.Hatch';
    L(:,2:4)=squeeze(handles.Molt(2:4,2,:)-handles.Molt(1:3,2,:))';
    I(:,1)=squeeze(handles.Molt(1,1,:))-handles.Hatch';
    I(:,2:4)=squeeze(handles.Molt(2:4,1,:)-handles.Molt(1:3,2,:))';
    
    %statistics (t-test) of boxplots 
    y=zeros(size(handles.samples));
    strain1 = handles.Strain1Selection.String{handles.Strain1Selection.Value};
    strain2 = handles.Strain2Selection.String{handles.Strain2Selection.Value};
    y(strcmp(handles.samples,strain1))=1;
    y(strcmp(handles.samples,strain2))=2;
 
    %disp(['Please note that the type of statistical test depends on research design, distribution of data and type of variables',char(10),'Do not blindly use the statistics below',char(10),strain1,' vs. ',strain2,' (strain1 vs strain2, ','dT=strain2-strain1, ','percentage=dT/Tstrain1)'])
    %ind=y~=0 & handles.valid;
    
        %molts
        delta_time_median_molt = (median(M(y==2 & handles.valid,:))-median(M(y==1 & handles.valid,:)))*handles.sampleRate;
        delta_time_mean_molt = (mean(M(y==2 & handles.valid,:))-mean(M(y==1 & handles.valid,:)))*handles.sampleRate;
        perc_median_molt =(delta_time_median_molt(:)'./(median(M(y==1 & handles.valid,:))*handles.sampleRate))*100;
        perc_mean_molt =(delta_time_mean_molt(:)'./(mean(M(y==1 & handles.valid,:))*handles.sampleRate))*100;
            %Wilcoxon test
            %p = ranksum(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i));
            
        sign_ttest_molt={};
        sign_vartest_molt={};
        for i=1:4
            [h,p_molt] = ttest2(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i));
            [l,v_molt] = vartest2(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i)); %v: p-value of variance
            p_ttest_molt(i) = p_molt;
            p_vartest_molt(i) = v_molt;
            if p_molt<0.001
                sign_ttest_molt{i} = '***';
            elseif p_molt==0.001
                sign_ttest_molt{i} = '***';
            elseif (p_molt<=0.01) && (p_molt>0.001)
                sign_ttest_molt{i} = '**';
            elseif p_molt>0.05
                sign_ttest_molt{i} = 'ns';
            elseif (p_molt<=0.05) && (p_molt>0.01)
                sign_ttest_molt{i} = '*';
            end
            
            if v_molt<0.001
                sign_vartest_molt{i} = '***';
            elseif v_molt==0.001
                sign_vartest_molt{i} = '***';
            elseif (v_molt<=0.01) && (v_molt>0.001)
                sign_vartest_molt{i} = '**';
            elseif v_molt>0.05
                sign_vartest_molt{i} = 'ns';
            elseif (v_molt<=0.05) && (v_molt>0.01)
                sign_vartest_molt{i} = '*';
            end
        end
        
        %intermolts
        delta_time_median_intermolt = (median(I(y==2 & handles.valid,:))-median(I(y==1 & handles.valid,:)))*handles.sampleRate;
        delta_time_mean_intermolt = (mean(I(y==2 & handles.valid,:))-mean(I(y==1 & handles.valid,:)))*handles.sampleRate;
        perc_median_intermolt =(delta_time_median_intermolt(:)'./(median(I(y==1 & handles.valid,:))*handles.sampleRate))*100;
        perc_mean_intermolt =(delta_time_mean_intermolt(:)'./(mean(I(y==1 & handles.valid,:))*handles.sampleRate))*100;
                     
        sign_ttest_intermolt={};
        sign_vartest_intermolt={};
        for i=1:4
            [h,p_intermolt] = ttest2(I(y==1 & handles.valid,i),I(y==2 & handles.valid,i));
            [l,v_intermolt] = vartest2(I(y==1 & handles.valid,i),I(y==2 & handles.valid,i)); %v: p-value of variance
            p_ttest_intermolt(i) = p_intermolt;
            p_vartest_intermolt(i) = v_intermolt;
            if p_intermolt<0.001
                sign_ttest_intermolt{i} = '***';
            elseif p_intermolt==0.001
                sign_ttest_intermolt{i} = '***';
            elseif (p_intermolt<=0.01) && (p_intermolt>0.001)
                sign_ttest_intermolt{i} = '**';
            elseif p_intermolt>0.05
                sign_ttest_intermolt{i} = 'ns';
            elseif (p_intermolt<=0.05) && (p_intermolt>0.01)
                sign_ttest_intermolt{i} = '*';
            end
            
            if v_intermolt<0.001
                sign_vartest_intermolt{i} = '***';
            elseif v_intermolt==0.001
                sign_vartest_intermolt{i} = '***';
            elseif (v_intermolt<=0.01) && (v_intermolt>0.001)
                sign_vartest_intermolt{i} = '**';
            elseif v_intermolt>0.05
                sign_vartest_intermolt{i} = 'ns';
            elseif (v_intermolt<=0.05) && (v_intermolt>0.01)
                sign_vartest_intermolt{i} = '*';
            end
        end
        
         %larvalstage
        delta_time_median_larvalstage = (median(L(y==2 & handles.valid,:))-median(L(y==1 & handles.valid,:)))*handles.sampleRate;
        delta_time_mean_larvalstage = (mean(L(y==2 & handles.valid,:))-mean(L(y==1 & handles.valid,:)))*handles.sampleRate;
        perc_median_larvalstage =(delta_time_median_larvalstage(:)'./(median(L(y==1 & handles.valid,:))*handles.sampleRate))*100;
        perc_mean_larvalstage =(delta_time_mean_larvalstage(:)'./(mean(L(y==1 & handles.valid,:))*handles.sampleRate))*100;
                     
        sign_ttest_larvalstage={};
        sign_vartest_larvalstage={};
        for i=1:4
            [h,p_larvalstage] = ttest2(L(y==1 & handles.valid,i),L(y==2 & handles.valid,i));
            [l,v_larvalstage] = vartest2(L(y==1 & handles.valid,i),L(y==2 & handles.valid,i)); %v: p-value of variance
            p_ttest_larvalstage(i) = p_larvalstage;
            p_vartest_larvalstage(i) = v_larvalstage;
            if p_larvalstage<0.001
                sign_ttest_larvalstage{i} = '***';
            elseif p_larvalstage==0.001
                sign_ttest_larvalstage{i} = '***';
            elseif (p_larvalstage<=0.01) && (p_larvalstage>0.001)
                sign_ttest_larvalstage{i} = '**';
            elseif p_larvalstage>0.05
                sign_ttest_larvalstage{i} = 'ns';
            elseif (p_larvalstage<=0.05) && (p_larvalstage>0.01)
                sign_ttest_larvalstage{i} = '*';
            end
            
            if v_larvalstage<0.001
                sign_vartest_larvalstage{i} = '***';
            elseif v_larvalstage==0.001
                sign_vartest_larvalstage{i} = '***';
            elseif (v_larvalstage<=0.01) && (v_larvalstage>0.001)
                sign_vartest_larvalstage{i} = '**';
            elseif v_larvalstage>0.05
                sign_vartest_larvalstage{i} = 'ns';
            elseif (v_larvalstage<=0.05) && (v_larvalstage>0.01)
                sign_vartest_larvalstage{i} = '*';
            end
        end
        
        
        
        %for i=1:4
            %Wilcoxon test
            %p = ranksum(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i)); 

            %two-sample t-test
            %[h,p_molt] = ttest2(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i));
            %[l,v_molt] = vartest2(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i)); %v: p-value of variance
            %if p<0.001
                %disp(['Molt ',num2str(i),': dT_median = ',num2str(delta_time_median(i),'%.2f'),' hrs (',num2str(perc_median(i), '%.1f'),'%), dT_mean = ',num2str(delta_time_mean(i),'%.2f'),' hrs (',num2str(perc_mean(i), '%.1f'),'%), ','t-test (p < 0.001) ***, ', 'F-test (p = ',num2str(v,'%.3f'),')'])
            %elseif p==0.001
                %disp(['Molt ',num2str(i),': dT_median = ',num2str(delta_time_median(i),'%.2f'),' hrs (',num2str(perc_median(i), '%.1f'),'%), dT_mean = ',num2str(delta_time_mean(i),'%.2f'),' hrs (',num2str(perc_mean(i), '%.1f'),'%), ','t-test (p = ',num2str(p,'%.3f'), ') ***, ', 'F-test (p = ',num2str(v,'%.3f'),')'])
            %elseif (p<=0.01) && (p>0.001)
                %disp(['Molt ',num2str(i),': dT_median = ',num2str(delta_time_median(i),'%.2f'),' hrs (',num2str(perc_median(i), '%.1f'),'%), dT_mean = ',num2str(delta_time_mean(i),'%.2f'),' hrs (',num2str(perc_mean(i), '%.1f'),'%), ','t-test (p = ',num2str(p,'%.3f'), ') **, ', 'F-test (p = ',num2str(v,'%.3f'),')'])
            %elseif p>0.05
                %disp(['Molt ',num2str(i),': dT_median = ',num2str(delta_time_median(i),'%.2f'),' hrs (',num2str(perc_median(i), '%.1f'),'%), dT_mean = ',num2str(delta_time_mean(i),'%.2f'),' hrs (',num2str(perc_mean(i), '%.1f'),'%), ','t-test (p = ',num2str(p,'%.3f'), ') ns, ', 'F-test (p = ',num2str(v,'%.3f'),')'])
            %elseif (p<=0.05) && (p>0.01)
                %disp(['Molt ',num2str(i),': dT_median = ',num2str(delta_time_median(i),'%.2f'),' hrs (',num2str(perc_median(i), '%.1f'),'%), dT_mean = ',num2str(delta_time_mean(i),'%.2f'),' hrs (',num2str(perc_mean(i), '%.1f'),'%), ','t-test )p = ',num2str(p,'%.3f'), ') *, ', 'F-test (p = ',num2str(v,'%.3f'),')'])
            %end
        %end

       
        


labels = {'DevStage','dT_median(hrs)','dT_median(%)','dT_mean(hrs)','dT_mean(%)','pvalue_t-test','sign_t-test','pvalue_var-test','sign_var-test'};
DevStage = {'Molt1','Molt2','Molt3','Molt4','Intermolt1','Intermolt2','Intermolt3','Intermolt4','LarvalStage1','LarvalStage2','LarvalStage3','Larvalstage4'};
dT_median = num2cell([delta_time_median_molt delta_time_median_intermolt delta_time_median_larvalstage]);
dT_median_perc = num2cell([perc_median_molt perc_median_intermolt perc_median_larvalstage]);
dT_mean = num2cell([delta_time_mean_molt delta_time_mean_intermolt delta_time_mean_larvalstage]);
dT_mean_perc = num2cell([perc_mean_molt perc_mean_intermolt perc_mean_larvalstage]);
pval_ttest = num2cell([p_ttest_molt p_ttest_intermolt p_ttest_larvalstage]);
sign_ttest = ([sign_ttest_molt sign_ttest_intermolt sign_ttest_larvalstage]);
pval_vartest = num2cell([p_vartest_molt p_vartest_intermolt p_vartest_larvalstage]);
sign_vartest = ([sign_vartest_molt sign_vartest_intermolt sign_vartest_larvalstage]);
results = [DevStage;dT_median;dT_median_perc;dT_mean;dT_mean_perc;pval_ttest;sign_ttest;pval_vartest;sign_vartest]';

fileID = fopen([handles.filepath,handles.filename(1:end),'_DevStages_stats_',strain1,'_vs_',strain2,'.csv'],'w');
formatSpec = '%s\t%.6f\t%.6f\t%.6f\t%.6f\t%e\t%s\t%e\t%s\n';
formatSpec2 = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'; % t=tab seperated, %s=string

table = [labels;results];
[nrowTable,ncolTable] = size(table);

fprintf(fileID,formatSpec2,table{1,:})
for row = 2:nrowTable
    fprintf(fileID,formatSpec,table{row,:});
end
fclose(fileID);
            
            
           


% --- Executes on button press in CheckAnnotation.
function CheckAnnotation_Callback(hObject, eventdata, handles)
% hObject    handle to CheckAnnotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%this function checks whether annotation of molts occured correctly, i.e. whether all durations are > 0 
%wrong annotation may occur when use of GUI buttons did not happen in
    %correct order, or when not all molts are annotated
Growth_durations = [squeeze(handles.Molt(1,1,:))'-handles.Hatch ;squeeze(handles.Molt(2:4,1,:)-handles.Molt(1:3,2,:)) ];
if any(any(Growth_durations(:,handles.valid)<0))
    [~,c] = find(any(Growth_durations<0) & handles.valid);
   disp(['Error, Sample ', num2str(c(1)),' wrongly annotated'])
    handles.current_sample=c(1);
    plot_data(handles)
    guidata(hObject, handles);
    save_data(handles)
else
valid = handles.valid;
disp('all samples are correct')
end


% --- Executes on selection change in Strain1Selection.
function Strain1Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Strain1Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Strain1Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Strain1Selection


% --- Executes during object creation, after setting all properties.
function Strain1Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Strain1Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Strain2Selection.
function Strain2Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Strain2Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Strain2Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Strain2Selection


% --- Executes during object creation, after setting all properties.
function Strain2Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Strain2Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in SetHatchToZero.
function SetHatchToZero_Callback(hObject, eventdata, handles)
% hObject    handle to SetHatchToZero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of SetHatchToZero

% to set the hatch to zero or to annoatate the hatch for ALL samples
if handles.SetHatchToZero.Value==1 %if box is checked, the hatch of all samples is set to zero
    handles.Hatch(:)=0;
else %if box is unchecked, annotate the hatch for valid samples only
        for i=1:size(handles.X,2)
            if handles.valid(i)==1
            try;handles.Hatch(i)=find(1:length(handles.X(:,i))>3 & handles.X(:,i)'>(mean(handles.X(1:20,i))+5*std(handles.X(1:20,i))),1,'first')-1;end
            end 
        end
end
plot_data(handles)
guidata(hObject, handles);
save_data(handles)


function yMax_Callback(hObject, eventdata, handles)
% hObject    handle to yMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%function to adjust the max value of the y-axis by typing a number
%without a number the ylim will be adjusted according to the max value of
    %the data for the current sample

c = handles.current_sample;
if ~isnan(str2double(eventdata.Source.String))
     %handles.axis_lim = [1 t(end) 0 str2double(eventdata.Source.String)];
     handles.ymaxlim = str2double(eventdata.Source.String);
else
     %handles.axis_lim = [1 t(end) 0 1.1*max(handles.X(:,c))];
     handles.ymaxlim = 1.1*max(handles.X(:,c));
end
guidata(hObject, handles);
plot_data(handles)
 
% Hints: get(hObject,'String') returns contents of yMax as text
%        str2double(get(hObject,'String')) returns contents of yMax as a double


% --- Executes during object creation, after setting all properties.
function yMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LumiLevels.
function LumiLevels_Callback(hObject, eventdata, handles)
% hObject    handle to LumiLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% first let pop-up appear to choose strains for the boxplot
uniq_samples_stable = unique(handles.samples, 'stable'); %order is as in _pd (handles.uniq_samples is alphabetically/numerically ordered)
    Prompt{1}= 'Strains  to be plotted: ';
    Formats(1).type = 'list';
    Formats(1).style = 'listbox';
    Formats(1).items = uniq_samples_stable;
    Formats(1).format = 'text';
    Formats(1).limits = [0 length(Formats(1).items)];
    Formats(1).size = [200 length(Formats(1).items)*20];
    [cfg,Cancelled] = inputsdlg(Prompt,'Select Strains',Formats); %cfg=number indicating selected strains

% let popup appear to adapt the order of the strains to be plotted
    %enter a 1 for the strain to be plotted first in boxplot, a 2 for the
    %strain to be plotted second, a 3 for the third, and so on
    
    strains_name = uniq_samples_stable(cfg{1}); %strains to be plotted

    for i=1:length(strains_name)
        Prompt_order{i} = strains_name{i};
        Formats_order(i).type = 'list';
        Formats_order(i).style = 'popupmenu';
        Formats_order(i).items = num2cell(1:length(strains_name));
        Formats_order(i).format = 'text';
        Formats_order(i).limits = [0 length(1:length(strains_name))];
        Formats_order(i).size = [200 length(1:length(strains_name))*20];
    end
        strains_orderID = inputsdlg(Prompt_order','Strain order',Formats_order)'; %indicates for each strain on which spot you want it to be: 1 for first, 2 for second
        
        %transform spot into column position
        for i = 1:length(strains_orderID)
            strains_orderPos(strains_orderID{i}) = i;
        end 
        %order strains
        strains_name_ordered = strains_name(strains_orderPos);
    
%sort handles.samples by strain order
       samples_orderID = (zeros(size(handles.samples)));
    for i=1:length(strains_orderID)
       samples_orderID(strcmp(handles.samples,strains_name_ordered(i)))=i; %couple the sample idx to the order
    end
        [samples_orderID_sorted,idx_order] = sort(samples_orderID);
        samples_ordered = handles.samples(idx_order);
        valid_ordered = handles.valid(idx_order);

%sort handles.Molt and handles.Hatch
    Molt_ordered = handles.Molt(:,:,idx_order);
    X_ordered = handles.X(:,idx_order);

%determine the median luminescence during the molt
lumi=[];
j=1;
    for i=1:length(valid_ordered)
        if valid_ordered(i)==1
            lumi(j,1) = log10(median(X_ordered(Molt_ordered(1,1,i):Molt_ordered(1,2,i),i)));
            lumi(j,2) = log10(median(X_ordered(Molt_ordered(2,1,i):Molt_ordered(2,2,i),i)));
            lumi(j,3) = log10(median(X_ordered(Molt_ordered(3,1,i):Molt_ordered(3,2,i),i)));
            lumi(j,4) = log10(median(X_ordered(Molt_ordered(4,1,i):Molt_ordered(4,2,i),i)));
        end
        j=j+1;
    end

%define boxplot positions, width of boxplot is 0.5
    % and distance between the center of 2 groups of boxplots is 0.25 
box_positions = [];
width = length(cfg{1})+0.5;
is_even = rem(length(cfg{1}),2)==0;
for i=1:4;
    if is_even
        shift = (i-1)*width + 0.5;
        box_positions = cat(2,box_positions, (1:length(cfg{1}))-0.5+shift);
    else
        shift = (i-1)*width + 1;
        box_positions = cat(2,box_positions, (1:length(cfg{1}))-0.5+shift); 
    end
end

    %use selected strains and valid sampled
    ind = logical(zeros(size(samples_ordered))); 
    for i=1:length(strains_name_ordered)
       ind(strcmp(samples_ordered,strains_name_ordered(i)) & valid_ordered)=1;
    end

    %generate x-labels for te boxes
    labels=repmat({'',' ','  ','   '},sum(ind),1);
    labels=labels(:)';
    final_labels = strcat(repmat(samples_ordered(ind),1,4),labels);
    clear labels
    if is_even
        tick_positions = box_positions(ceil(length(cfg{1})/2)+[0:3]*length(cfg{1}))+0.5;
    else
        tick_positions = box_positions(ceil(length(cfg{1})/2)+[0:3]*length(cfg{1}));
    end
    
    %generate boxplots
    tmp_Lumi = lumi(ind,:);
    figure;boxplot(tmp_Lumi(:),final_labels,'ColorGroup',repmat(samples_ordered(ind),1,4),...
       'positions',box_positions,'LabelOrientation','inline');title('Luminescence during Molt');ylabel('Median luminescence (log10)');axis([0.5 box_positions(end)+0.5 0 1.25*max(tmp_Lumi(:))])
    hold on;for i=1:4;text(tick_positions(i),0.15*min(tmp_Lumi(:)),['M',num2str(i)],'HorizontalAlignment','Center');end


% --- Executes on button press in StatLuminescence.
function StatLuminescence_Callback(hObject, eventdata, handles)
% hObject    handle to StatLuminescence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
    

    y=zeros(size(handles.samples));
    strain1 = handles.Strain1Selection.String{handles.Strain1Selection.Value};
    strain2 = handles.Strain2Selection.String{handles.Strain2Selection.Value};
    y(strcmp(handles.samples,strain1))=1;
    y(strcmp(handles.samples,strain2))=2;
    
    lumi=[];
    j=1;
    for i=1:length(handles.valid)
        if handles.valid(i)==1
            lumi(j,1) = log10(median(handles.X(handles.Molt(1,1,i):handles.Molt(1,2,i),i)));
            lumi(j,2) = log10(median(handles.X(handles.Molt(2,1,i):handles.Molt(2,2,i),i)));
            lumi(j,3) = log10(median(handles.X(handles.Molt(3,1,i):handles.Molt(3,2,i),i)));
            lumi(j,4) = log10(median(handles.X(handles.Molt(4,1,i):handles.Molt(4,2,i),i)));
        end
        j=j+1;
    end
    
    %molts
        delta_median_lumi = (median(lumi(y==2 & handles.valid,:))-median(lumi(y==1 & handles.valid,:)));
        delta_mean_lumi = (mean(lumi(y==2 & handles.valid,:))-mean(lumi(y==1 & handles.valid,:)));
        perc_median_lumi =(delta_median_lumi(:)'./(median(lumi(y==1 & handles.valid,:))))*100;
        perc_mean_lumi =(delta_mean_lumi(:)'./(mean(lumi(y==1 & handles.valid,:))))*100;
            %Wilcoxon test
            %p = ranksum(M(y==1 & handles.valid,i),M(y==2 & handles.valid,i));
            
        sign_ttest_lumi={};
        sign_vartest_lumi={};
        for i=1:4
            [h,p_lumi] = ttest2(lumi(y==1 & handles.valid,i),lumi(y==2 & handles.valid,i));
            [l,v_lumi] = vartest2(lumi(y==1 & handles.valid,i),lumi(y==2 & handles.valid,i)); %v: p-value of variance
            p_ttest_lumi(i) = p_lumi;
            p_vartest_lumi(i) = v_lumi;
            if p_lumi<0.001
                sign_ttest_lumi{i} = '***';
            elseif p_lumi==0.001
                sign_ttest_lumi{i} = '***';
            elseif (p_lumi<=0.01) && (p_lumi>0.001)
                sign_ttest_lumi{i} = '**';
            elseif p_lumi>0.05
                sign_ttest_lumi{i} = 'ns';
            elseif (p_lumi<=0.05) && (p_lumi>0.01)
                sign_ttest_lumi{i} = '*';
            end
            
            if v_lumi<0.001
                sign_vartest_lumi{i} = '***';
            elseif v_lumi==0.001
                sign_vartest_lumi{i} = '***';
            elseif (v_lumi<=0.01) && (v_lumi>0.001)
                sign_vartest_lumi{i} = '**';
            elseif v_lumi>0.05
                sign_vartest_lumi{i} = 'ns';
            elseif (v_lumi<=0.05) && (v_lumi>0.01)
                sign_vartest_lumi{i} = '*';
            end
        end
        
labels = {'DevStage','FC_medLumiMolt_median(log10)','FC_medLumiMolt_median(%)','FC_medLumiMolt_mean(log10)','FC_medLumiMolt_mean(%)','pvalue_t-test','sign_t-test','pvalue_var-test','sign_var-test'};
DevStage = {'Molt1','Molt2','Molt3','Molt4'};
d_lumi_median = num2cell(delta_median_lumi);
d_lumi_median_perc = num2cell(perc_median_lumi);
d_lumi_mean = num2cell(delta_mean_lumi);
d_lumi_mean_perc = num2cell(perc_mean_lumi);
pval_ttest = num2cell(p_ttest_lumi);
sign_ttest = sign_ttest_lumi;
pval_vartest = num2cell(p_vartest_lumi);
sign_vartest = sign_vartest_lumi;
results = [DevStage;d_lumi_median;d_lumi_median_perc;d_lumi_mean;d_lumi_mean_perc;pval_ttest;sign_ttest;pval_vartest;sign_vartest]';

fileID = fopen([handles.filepath,handles.filename(1:end),'_log10medLumiMolt_stats_',strain1,'_vs_',strain2,'.csv'],'w');
formatSpec = '%s\t%.6f\t%.6f\t%.6f\t%.6f\t%e\t%s\t%e\t%s\n';
formatSpec2 = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'; % t=tab seperated, %s=string

table = [labels;results];
[nrowTable,ncolTable] = size(table);

fprintf(fileID,formatSpec2,table{1,:})
for row = 2:nrowTable
    fprintf(fileID,formatSpec,table{row,:});
end
fclose(fileID);
    

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in strain1Lumi.
function strain1Lumi_Callback(hObject, eventdata, handles)
% hObject    handle to strain1Lumi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function strain1Lumi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strain1Lumi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in strain2Lumi.
function strain2Lumi_Callback(hObject, eventdata, handles)
% hObject    handle to strain2Lumi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function strain2Lumi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strain2Lumi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
