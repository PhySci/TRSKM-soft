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

% Last Modified by GUIDE v2.5 06-Sep-2015 19:12:47

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



% --- Outputs from this function are returned to the command line.
function varargout = phaser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function loadAct_Callback(hObject, eventdata, handles)
% hObject    handle to loadAct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  global ch1 ch2 distance rMax fName fPath
  [fName, fPath] = uigetfile('*.*');
  
  if (fName == 0)
    return
  end
  fullName = strcat(fPath,fName);
  set(handles.InputFile,'String',strcat('Input file: ',fullName));
  M = dlmread(fullName);
  distance = M(:,1);
  ch1 = M(:,4);
  ch2 = M(:,5);
  rMax = max(sqrt(ch1.*ch1+ch2.*ch2));
  
  plot(handles.axes1,distance,ch1,distance,ch2);
  legend(handles.axes1,{'Chanel1', 'Chanel2'});
 
  plot(handles.axes3,ch1,ch2);
  xlim(handles.axes3,[-1.1*rMax 1.1*rMax]);
  ylim(handles.axes3,[-1.1*rMax 1.1*rMax]);

  
  xlabel(handles.axes3,'Chanel 1');
  ylabel(handles.axes3,'Chanel 2');
  xlabel(handles.axes4,'Chanel 1');
  ylabel(handles.axes4,'Chanel 2');


% --------------------------------------------------------------------
function saveAct_Callback(hObject, eventdata, handles)
  global distance ch1Correct ch2Correct fName fPath
  [~, name, ext] = fileparts(fName);
  defaultName = strcat(fPath,name,'-corr',ext);
  [saveName,savePath] = uiputfile('*.txt','Save result as...',defaultName)
  dlmwrite(strcat(savePath,saveName), [distance ch1Correct ch2Correct], '\t');

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
  global distance ch1Correct ch2Correct rMax rMaxCorrect 
  angle = str2num(get(handles.angleVal,'String'));
  
  totalAngle = str2num(get(handles.angleInd,'String'));
  set(handles.angleInd,'String',num2str(totalAngle+angle));
  
  rotMatx = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];
  inp = [ch1Correct ch2Correct];
  out = (rotMatx*inp.').';
  ch1Correct = out(:,1);
  ch2Correct = out(:,2);
  
  plot(handles.axes4,ch1Correct,ch2Correct);
  xlim(handles.axes4,[-1.1*rMaxCorrect 1.1*rMaxCorrect]);
  ylim(handles.axes4,[-1.1*rMaxCorrect 1.1*rMaxCorrect]);
  xlabel(handles.axes4,'Chanel 1');
  ylabel(handles.axes4,'Chanel 2');
  
  plot(handles.axes2,distance,ch1Correct,distance,ch2Correct);
  legend(handles.axes2,{'Chanel1', 'Chanel2'});
