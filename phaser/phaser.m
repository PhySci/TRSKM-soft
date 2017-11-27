function varargout = phaser(varargin)
% PHASER M-file for phaser.fig
%      PHASER, by itself, creates a new PHASER or raises the existing
%      singleton*.
%
%      H = PHASER returns the handle to a new PHASER or the handle to
%      the existing singleton*.
%
%      PHASER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASER.M with the given input arguments.
%
%      PHASER('Property','Value',...) creates a new PHASER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phaser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phaser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phaser

% Last Modified by GUIDE v2.5 23-Mar-2016 16:42:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phaser_OpeningFcn, ...
                   'gui_OutputFcn',  @phaser_OutputFcn, ...
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


% --- Executes just before phaser is made visible.
function phaser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phaser (see VARARGIN)

% Choose default command line output for phaser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



function varargout = phaser_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
  global ch1 ch2 distance ch1Correct ch2Correct rMax rMaxCorrect
  cursorPoint = get(handles.axes3, 'CurrentPoint');
  curX = cursorPoint(1,1);
  curY = cursorPoint(1,2);
  if (abs(curX) < rMax && abs(curY) < rMax)
    set(handles.shiftX,'String',sprintf('%-5.2e',curX));
    set(handles.shiftY,'String',sprintf('%-5.2e',curY));
    ch1Correct = ch1 - curX;
    ch2Correct = ch2 - curY;
    rMaxCorrect = max(sqrt(ch1Correct.*ch1Correct+ch2Correct.*ch2Correct));
    
    plot(handles.axes2,distance,ch1Correct,distance,ch2Correct);
    legend(handles.axes2,{'Chanel1', 'Chanel2'});
    
    plot(handles.axes4,ch1Correct,ch2Correct);
    xlim(handles.axes4,[-1.1*rMaxCorrect 1.1*rMaxCorrect]);
    ylim(handles.axes4,[-1.1*rMaxCorrect 1.1*rMaxCorrect]); 
    xlabel(handles.axes4,'Chanel 1');
    ylabel(handles.axes4,'Chanel 2');
    set(handles.angleInd,'String','0');
  end;

function angleVal_Callback(hObject, eventdata, handles)
% hObject    handle to angleVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angleVal as text
%        str2double(get(hObject,'String')) returns contents of angleVal as a double

% --- Executes during object creation, after setting all properties.
function angleVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angleVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end


% --- Executes on button press in rotateBtn.
function rotateBtn_Callback(hObject, eventdata, handles)
% hObject    handle to rotateBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  global rawSignal distance correctSignal Monitors rMax rMaxCorrect 
  angle = str2num(get(handles.angleVal,'String'));
  
  totalAngle = str2num(get(handles.angleInd,'String'));
  set(handles.angleInd,'String',num2str(totalAngle+angle));
  
  rotMatx = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];
  inp = [real(correctSignal) imag(correctSignal)];
  out = (rotMatx*inp.').';
  correctSignal = out(:,1) + i*out(:,2);
  
  plotCorrectSignal(handles)
  
% Callback for load h5 file
function LoadH5_Callback(hObject, eventdata, handles)
    global rawSignal distance correctSignal Monitors rMax fullName timeScale freqScale    
    % open file dialog 
    [fName,fPath,~] = uigetfile({'*.h5';'*.*'});
    if (fName == 0)
        return
    else
        fullName = fullfile(fPath,fName);
        set(handles.InputFile,'String',fullName);    
        % read array of results
        res = h5read(fullName,'/Signal').';
        
        % check if it was backward movement
        if res(1,1)>res(end,1)
            res = flipud(res);
        end
        
        distance = res(:,1);
        rawSignal =  res(:,4)+ i*res(:,5);
        correctSignal =  res(:,4)+ i*res(:,5);
        rMax = max(abs(correctSignal ));
        Monitors = res(:,2)+i*res(:,3);
        
        timeScale = 2*(distance - min(distance))/3e2; % m
        dt = abs(mean(diff(timeScale)));
        freqScale = linspace(-0.5/dt,0.5/dt,size(timeScale,1)).';  
        
        plotRawSignal(handles);
        plotCorrectSignal(handles);
        
    end
    
function SaveH5_Callback(hObject, eventdata, handles)
    global distance correctSignal fullName
    [fName,fPath,~] = uiputfile('./*.h5','Save as...',fullName);
    
    saveName = fullfile(fPath,fName);
    
    res = [distance, real(correctSignal), imag(correctSignal)];
    path = '/correctedSignal';
    try
      h5create(saveName,path,[301 3],'ChunkSize',[301 3]);  
    catch
    end
    h5write(saveName,path,res);

% --------------------------------------------------------------------
function Main_Callback(hObject, eventdata, handles)
% hObject    handle to Main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function plotRawSignal(handles)
    global rawSignal distance correctSignal Monitors rMax timeScale freqScale 
    % plot time-resolved signal
    plot(handles.axes1,timeScale,real(rawSignal),timeScale,imag(rawSignal));
    legend(handles.axes1,'Chanel 1','Chanel 2');
    
    % plot complex plot
    plot(handles.axes3,real(rawSignal),imag(rawSignal));
    xlim(handles.axes3,[-1.1*rMax 1.1*rMax]);
    ylim(handles.axes3,[-1.1*rMax 1.1*rMax]);
    xlabel(handles.axes3,'Chanel 1');
    ylabel(handles.axes3,'Chanel 2');
    
    % calculate and plot FFT
    FFT1 = abs(fftshift(fft(real(rawSignal))));
    FFT2 = abs(fftshift(fft(imag(rawSignal))));
    plot(handles.axes5,freqScale,FFT1,freqScale,FFT2);
    xlim(handles.axes5,[0 25]);
    xlabel(handles.axes5,'Frequency (GHz)');
    legend(handles.axes5,'Chanel 1','Chanel 2');
    
function plotCorrectSignal(handles)
    global distance correctSignal Monitors rMax freqScale timeScale
    
    plot(handles.axes2,timeScale,real(correctSignal),timeScale,imag(correctSignal));
    legend(handles.axes2,{'Chanel1', 'Chanel2'});
    
    % complex plot
    plot(handles.axes4,real(correctSignal),imag(correctSignal));
    xlim(handles.axes4,[-1.1*rMax 1.1*rMax]);
    ylim(handles.axes4,[-1.1*rMax 1.1*rMax]);
    xlabel(handles.axes4,'Chanel 1');
    ylabel(handles.axes4,'Chanel 2');
    
    % calculate and plot FFT
    FFT1 = abs(fftshift(fft(real(correctSignal))));
    FFT2 = abs(fftshift(fft(imag(correctSignal))));
    plot(handles.axes6,freqScale,FFT1,freqScale,FFT2);
    xlim(handles.axes6,[0 25]);
    xlabel(handles.axes6,'Frequency (GHz)');
    legend(handles.axes6,'Chanel 1','Chanel 2');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global correctSignal
    FFT1 = abs(fftshift(fft(real(correctSignal))));
    num2clip(FFT1);
