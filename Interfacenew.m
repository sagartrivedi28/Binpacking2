function varargout = RFL_Interface2(varargin)
% RFL_INTERFACE2 MATLAB code for RFL_Interface2.fig
%      RFL_INTERFACE2, by itself, creates a new RFL_INTERFACE2 or raises the existing
%      singleton*.
%
%      H = RFL_INTERFACE2 returns the handle to a new RFL_INTERFACE2 or the handle to
%      the existing singleton*.
%
%      RFL_INTERFACE2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RFL_INTERFACE2.M with the given input arguments.
%
%      RFL_INTERFACE2('Property','Value',...) creates a new RFL_INTERFACE2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RFL_Interface2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RFL_Interface2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RFL_Interface2

% Last Modified by GUIDE v2.5 23-Aug-2016 17:17:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RFL_Interface2_OpeningFcn, ...
                   'gui_OutputFcn',  @RFL_Interface2_OutputFcn, ...
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


% --- Executes just before RFL_Interface2 is made visible.
function RFL_Interface2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RFL_Interface2 (see VARARGIN)

% Choose default command line output for RFL_Interface2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RFL_Interface2 wait for user response (see UIRESUME)
% uiwait(handles.RFL);


% --- Outputs from this function are returned to the command line.
function varargout = RFL_Interface2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_fname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fname as text
%        str2double(get(hObject,'String')) returns contents of edit_fname as a double


% --- Executes during object creation, after setting all properties.
function edit_fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Fload.
function pushbutton_Fload_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Fload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.fname,handles.pname] = uigetfile('*.xls;*.xlsx','Select The Product List File');
handles.xlspath = fullfile(handles.pname,handles.fname);
set(handles.edit_fname,'String',handles.xlspath);
guidata(hObject,handles);


function edit_Maxno_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Maxno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Maxno as text
%        str2double(get(hObject,'String')) returns contents of edit_Maxno as a double


% --- Executes during object creation, after setting all properties.
function edit_Maxno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Maxno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_grid_Callback(hObject, eventdata, handles)
% hObject    handle to edit_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_grid as text
%        str2double(get(hObject,'String')) returns contents of edit_grid as a double


% --- Executes during object creation, after setting all properties.
function edit_grid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_wsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wsize as text
%        str2double(get(hObject,'String')) returns contents of edit_wsize as a double


% --- Executes during object creation, after setting all properties.
function edit_wsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ArrangeDice.
function pushbutton_ArrangeDice_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ArrangeDice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); cla;
axes(handles.axes2); cla;
Rswitch = get(handles.RFlag,'Value');
MPW_flag = get(handles.MPW_Flag,'Value');

handles.Maxno = str2num(get(handles.edit_Maxno,'String'));
handles.gsize = str2num(get(handles.edit_grid,'String'));
handles.gsize = handles.gsize*10^-3;
handles.wsize = str2num(get(handles.edit_wsize,'String'));
handles.slwidth = str2num(get(handles.edit_slwidth,'String'));
if length(handles.slwidth)==1
    handles.slwidth = [handles.slwidth handles.slwidth];
end;
handles.Rxmax = str2num(get(handles.edit_Rxmax,'String'));
handles.Rymax = str2num(get(handles.edit_Rymax,'String'));
if MPW_flag
    
    [~,~,Dice_data1] = xlsread(handles.xlspath,1,['B2:I' num2str(handles.Maxno+1)]);    
    % B = ID, C=X-Size, D=Y-Size, E=X-Loc, F=Y-Loc, G=DVol
    FDice_data = Dice_data1;
    Fcor = [cell2mat(Dice_data1(:,4))-cell2mat(Dice_data1(:,2))/2-handles.slwidth(1)/2, cell2mat(Dice_data1(:,5))-cell2mat(Dice_data1(:,3))/2-handles.slwidth(2)/2,...
            cell2mat(Dice_data1(:,4))+cell2mat(Dice_data1(:,2))/2+handles.slwidth(1)/2, cell2mat(Dice_data1(:,5))+cell2mat(Dice_data1(:,3))/2+handles.slwidth(2)/2];
    Fvol = cell2mat(Dice_data1(:,8));
    Fid = (1:size(Dice_data1,1))';
    Ddim = [cell2mat(Dice_data1(:,2))+handles.slwidth(1) cell2mat(Dice_data1(:,3))+handles.slwidth(2)];
    Ddim = ceil(Ddim/handles.gsize)*handles.gsize;
    FDice_data = Dice_data1;
    
    Gpara.Rmax = [handles.Rxmax handles.Rymax]*10^3;
    Gpara.Sline = handles.slwidth;
    Gpara.Rno=ones(handles.Maxno,1);
    for i=2:size(Dice_data1,1)
         if strcmp(Dice_data1{i-1},Dice_data1{i})
              Gpara.Rno(i)=Gpara.Rno(i-1);
         else
              Gpara.Rno(i)=Gpara.Rno(i-1)+1;
         end
    end
    [Fsets,Fwafer] = MPW_DicePacking(Fcor,Ddim(:,1),Ddim(:,2),Fvol,Gpara,handles.axes1);
    
    Lcor = [Fcor(:,1)+Ddim(:,1)/2,Fcor(:,2)+Ddim(:,2)/2];
    FDdim = [cell2mat(Dice_data1(:,2)), cell2mat(Dice_data1(:,3))];
    Rnew = [max(Fcor(:,3))-min(Fcor(:,1)) max(Fcor(:,4))-min(Fcor(:,2))]/10^3;
    Dstruct = DieOffset(Rnew(1),Rnew(2),handles.wsize);
    GDPW = sum(Dstruct.map_in(:));    
    Plotdata([],Dstruct,[],handles.axes2);

else

    [~,~,Dice_data1] = xlsread(handles.xlspath,1,['B2:K' num2str(handles.Maxno+1)]);
    
    Fwafer = Inf;
    GDPW = 1;
    VOLmin = cell2mat(Dice_data1(:,4)); VOLmin = min(VOLmin(VOLmin~=0));
    Dsize = cell2mat(Dice_data1(:,2)).*cell2mat(Dice_data1(:,3));  %  Use size factor to repeate
    RVOL = floor(cell2mat(Dice_data1(:,4))/VOLmin);

    if ~Rswitch
        RVOL(:)=1;
    end;
    
    for ar_limit=0:max(RVOL)-1
        Rauto = zeros(size(RVOL));
        for ar=0:sum((RVOL(RVOL(:)>0)-1))
            Rflag=true;
            Dice_data = Dice_data1;
            NDice_data={};
            Rno = [];
        
            for r = 1:size(Dice_data,1)
                Dice_data{r,9}=Dice_data{r,9}+Rauto(r);
                Dice_data{r,4} = Dice_data{r,4}/Dice_data{r,9};
                NDice_data = [NDice_data; repmat(Dice_data(r,:),Dice_data{r,9},1)];
                Rno = [Rno; r*ones(Dice_data{r,9},1)];
            end;
            Dice_data = NDice_data;

            Did = (1:size(Dice_data,1))';
            Ddim = [cell2mat(Dice_data(:,2))+handles.slwidth(1) cell2mat(Dice_data(:,3))+handles.slwidth(2)];
            Ddim = ceil(Ddim/handles.gsize)*handles.gsize;

            Dvol = cell2mat(Dice_data(:,4));
            Dreq = cell2mat(Dice_data(:,5:10));

            Gpara.Rmax = [handles.Rxmax handles.Rymax]*10^3;
            Gpara.Sline = handles.slwidth;
            Gpara.Rno = Rno;        
            Gpara.wsize = []; %handles.wsize;
            Gpara.TBflag = false;
            if get(handles.radiobutton_testline,'Value')
                Gpara.TBflag = true;
                Gpara.TBcenter = handles.TBcenter;
                Gpara.TBcorner = handles.TBcorner;
            end;            
        
            Fstruct = DicePackingR6(Did,Ddim(:,1),Ddim(:,2),Dvol,Dreq,Gpara);
            Tcor = Fstruct.Fcor;
            Twafer = Fstruct.Fwafer;
        
            if sum(Twafer)==0
                break;
            else
                Rnew1 = [max(Tcor(:,3))-min(Tcor(:,1)) max(Tcor(:,4))-min(Tcor(:,2))]/10^3;
                Dstruct = DieOffset(Rnew1(1),Rnew1(2),handles.wsize);
                GDPW1 = sum(Dstruct.map_in(:));
            end;

            if sum(ceil(Twafer/GDPW1))<sum(ceil(Fwafer/GDPW))            
                Plotdata(Fstruct,Dstruct,handles.axes1,handles.axes2);
                Fcor= Fstruct.Fcor;          
                Fsets= Fstruct.Fsets;            
                Fwafer= Fstruct.Fwafer;            
                FDdim = [Fstruct.Ddimx-handles.slwidth(1) Fstruct.Ddimy-handles.slwidth(2)];
                Fid = Did;
                Fvol = Dvol;
                FDice_data=Dice_data;
                Rshift = [max(Fcor(:,3))-min(Fcor(:,1)) max(Fcor(:,4))-min(Fcor(:,2))]/2;
                Lcor = [Fcor(:,1)+Fstruct.Ddimx/2-Rshift(1),Fcor(:,2)+Fstruct.Ddimy/2-Rshift(2)];
                GDPW = sum(Dstruct.map_in(:)); 
                Rnew = Rnew1;
            end;

            for r = 1:size(Dice_data1,1)
                if ((Dice_data1{r,9}+Rauto(r))<RVOL(r)) && Rflag && (RVOL(r)-Rauto(r))==max(RVOL-Rauto) && Rauto(r)<ar_limit
                    Rauto(r)=Rauto(r)+1;
                    Rflag=false;
                end;
            end;
            if Rflag
               break; 
            end
            
        end;
    end

end;

% Options for Stepper/Normal Reticle
% RFL import function (Group Info + Same Cut-set Info) with decided coordinates

SetID = cell(length(Fid),1);
XLset = SetID;
XLFsets = Fsets';
MDPW = zeros(length(Fid),1);
for l=1:length(Fsets)
    for i=1:length(Fsets{l})
        SetID{Fsets{l}(i)}=[SetID{Fsets{l}(i)} l];
        MDPW(Fsets{l}(i)) = MDPW(Fsets{l}(i))+ceil(Fwafer(l)/GDPW)*GDPW;         % Provides # dices = GDPW x all cut-sets
        %Lcor(Fsets{l}(i),:) = [Fcor(Fsets{l}(i),1)-Rshift(1)+FDdim(Fsets{l}(i),1)/2,Fcor(Fsets{l}(i),2)-Rshift(2)+FDdim(Fsets{l}(i),2)/2];
        XLset{Fsets{l}(i)} = num2str(SetID{Fsets{l}(i)});
        if Fvol(Fsets{l}(i))==0
            XLFsets{l}(i)=0;
        end;
    end;
    XLFsets{l} = num2str(XLFsets{l}(XLFsets{l}~=0));
end;
 XLset(Fvol==0)={'0'};
 MDPW(Fvol==0)=0;
% GDPW = (cell2mat(cellfun(@length,SetID,'UniformOutput', false)).*sum(map_in(:))).*(Dvol~=0);

handles.xlspathout = fullfile(handles.pname,[handles.fname '_out.xlsx']);

Myxlwrite(handles.xlspathout,{'Dice No' 'Dice ID' 'Dice Xcor' 'Dice Ycor' 'Set ID' 'Dice-Yield' 'Volume'},1,'A1:G1');
Myxlwrite(handles.xlspathout,Fid,1,['A2:A' num2str(length(Fid)+1)]);
Myxlwrite(handles.xlspathout,Dice_data(:,1),1,['B2:B' num2str(length(Fid)+1)]);
Myxlwrite(handles.xlspathout,Lcor,1,['C2:D' num2str(length(Fid)+1)]);
Myxlwrite(handles.xlspathout,XLset,1,['E2:E' num2str(length(Fid)+1)]);
Myxlwrite(handles.xlspathout,MDPW,1,['F2:F' num2str(length(Fid)+1)]); 
Myxlwrite(handles.xlspathout,Fvol,1,['G2:G' num2str(length(Fid)+1)]);

% Wafer Count: Calculation need to refine
Myxlwrite(handles.xlspathout,{'Set No' 'Dice No ' 'Wafer/Set'},1,'I1:K1');
Myxlwrite(handles.xlspathout,[1:length(Fsets)]',1,['I2:I' num2str(length(Fsets)+1)]);
Myxlwrite(handles.xlspathout,XLFsets,1,['J2:J' num2str(length(Fsets)+1)]);
Myxlwrite(handles.xlspathout,ceil(Fwafer/GDPW),1,['K2:K' num2str(length(Fsets)+1)]);
Myxlwrite(handles.xlspathout,{'Total Wafer =',sum(ceil(Fwafer/GDPW))},1,['J' num2str(length(Fsets)+2) ':K' num2str(length(Fsets)+2)]);

Myxlwrite(handles.xlspathout,FDdim,1,['M2:N' num2str(length(Fid)+1)]);

% Myxlwrite(handles.xlspath,ceil(Dvol./GDPW),2,['H2:H' num2str(length(Did)+1)]);
% Myxlwrite(handles.xlspath,{'Total Wafer =', num2str(sum(ceil(Dvol(GDPW~=0)./GDPW(GDPW~=0))-1)+length(Fsets))},2,'J2:K2');

set(handles.edit_bsize,'String',num2str(Rnew));
set(handles.edit_bshift,'String',num2str([Dstruct.xoffset,Dstruct.yoffset]));
set(handles.edit_shotno,'String',num2str(GDPW));

guidata(hObject,handles);





function edit_bsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bsize as text
%        str2double(get(hObject,'String')) returns contents of edit_bsize as a double


% --- Executes during object creation, after setting all properties.
function edit_bsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bshift_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bshift as text
%        str2double(get(hObject,'String')) returns contents of edit_bshift as a double


% --- Executes during object creation, after setting all properties.
function edit_bshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_shotno_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shotno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shotno as text
%        str2double(get(hObject,'String')) returns contents of edit_shotno as a double


% --- Executes during object creation, after setting all properties.
function edit_shotno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shotno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when RFL is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to RFL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_slwidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slwidth as text
%        str2double(get(hObject,'String')) returns contents of edit_slwidth as a double


% --- Executes during object creation, after setting all properties.
function edit_slwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Rymax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Rymax as text
%        str2double(get(hObject,'String')) returns contents of edit_Rymax as a double


% --- Executes during object creation, after setting all properties.
function edit_Rymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Rxmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rxmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Rxmax as text
%        str2double(get(hObject,'String')) returns contents of edit_Rxmax as a double


% --- Executes during object creation, after setting all properties.
function edit_Rxmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rxmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RFlag.
function RFlag_Callback(hObject, eventdata, handles)
% hObject    handle to RFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RFlag


% --- Executes on button press in MPW_Flag.
function MPW_Flag_Callback(hObject, eventdata, handles)
% hObject    handle to MPW_Flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MPW_Flag


% --- Executes on button press in radiobutton_testline.
function radiobutton_testline_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_testline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton_testline,'Value')
    prompt = {'TL Width :','TL Height :','# TL Column :', '# TL Row :', 'Center TL Width :', 'Center TL Height :', '# TL Blocks :'};
    dlg_title = 'Enter Testline Details';
    num_lines = 1;
%     defaultans = {'20','hsv'};
    TBdata = inputdlg(prompt,dlg_title,num_lines);    
    
    handles.TBcorner(1) = str2double(TBdata{1});
    handles.TBcorner(2) = str2double(TBdata{2});
    handles.TBcorner(3) = str2double(TBdata{3});
    handles.TBcorner(4) = str2double(TBdata{4});
    
    handles.TBcenter(1) = str2double(TBdata{5});
    handles.TBceter(2) = str2double(TBdata{6});
    handles.TBcenter(3) = str2double(TBdata{7});
    
end;

guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton_testline


