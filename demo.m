function varargout = demo(varargin)
% DEMO MATLAB code for demo.fig
%      DEMO, by itself, creates a new DEMO or raises the existing
%      singleton*.
%
%      H = DEMO returns the handle to a new DEMO or the handle to
%      the existing singleton*.
%
%      DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMO.M with the given input arguments.
%
%      DEMO('Property','Value',...) creates a new DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before demo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help demo

% Last Modified by GUIDE v2.5 08-May-2018 08:46:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo_OpeningFcn, ...
                   'gui_OutputFcn',  @demo_OutputFcn, ...
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


% --- Executes just before demo is made visible.
function demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demo (see VARARGIN)

% Choose default command line output for demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile('test.wav','��ѡ�������ļ���');
x=audioread([pathname filename]);
save temp x


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load temp.mat;
sound(x,8000);
axes(handles.axes1);
plot(x);

axes(handles.axes2);
x = double(x);
x = x / max(abs(x));

%��������
FrameLen = 256;
FrameInc = 80;

amp1 = 10;%�����������
amp2 = 1;%������С����
zcr1 = 10;%�����ʴ�����ֵ
zcr2 = 5;%������С����ֵ

maxsilence = 30;  %���������
minlen  = 15;    % ���ʱ������
status  = 0;%0�����Σ�1���ȶΣ�2�����Σ�3������
count   = 0;%�����������ȡ�
silence = 0;%��������

%���������
tmp1  = enframe(x(1:end-1), FrameLen, FrameInc);
tmp2  = enframe(x(2:end)  , FrameLen, FrameInc);
signs = (tmp1.*tmp2)<0;
diffs = (tmp1 -tmp2)>0.02;
zcr   = sum(signs.*diffs, 2);

%�����ʱ����
amp = sum(abs(enframe(filter([1 -0.9375], 1, x), FrameLen, FrameInc)), 2);

%������������
amp1 = min(amp1, max(amp)/4);
amp2 = min(amp2, max(amp)/8);

%��ʼ�˵���
x1 = 0; 
x2 = 0;
for n=1:length(zcr)
   goto = 0;
   switch status
   case {0,1}                   % 0 = ����, 1 = ���ܿ�ʼ
      if amp(n) > amp1          % ȷ�Ž���������
         x1 = max(n-count-1,1);%���ñ��
         status  = 2;
         silence = 0;
         count   = count + 1;
      elseif amp(n) > amp2 | ... % ���ܴ���������
             zcr(n) > zcr2
         status = 1;
         count  = count + 1;
      else                       % ����״̬
         status  = 0;
         count   = 0;
      end
   case 2,                       % 2 = ������
      if amp(n) > amp2 | ...     % ������������
         zcr(n) > zcr2
         count = count + 1;
      else                       % ����������
         silence = silence+1;
         if silence < maxsilence % ����������������δ����
            count  = count + 1;
         elseif count < minlen   % ��������̫�̣���Ϊ������
            status  = 0;
            silence = 0;
            count   = 0;
         else                    % ��������
            status  = 3;%3�Ļ��ǽ����Ρ�
         end
      end
   case 3,
      break;
   end
end   

count = count-silence/2;
x2 = x1 + count -1;

plot(x)
ylabel('Speech');
line([x1*FrameInc x1*FrameInc], [-1 1], 'Color', 'red');
line([x2*FrameInc x2*FrameInc], [-1 1], 'Color', 'red');


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load temp.mat;
load mylabel.mat;
[x1 x2] = vad(x);
x=0.2*x/max(x);
m=mfcc(x);
m=m(x1-2:x2-2,:);
m=m';
m = vqlbg(m,4);
m = m';
m = reshape(m(1:4,:)',[1,96]);
m = [1 m];
elm_predict(m);
load elm_output.mat;
a=sprintf('%s', label(output,:));
set(handles.edit1,'string',a);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
recorder = audiorecorder(8000,16,1);
disp('Start speaking.')
recordblocking(recorder, 3);
disp('End of Recording.');
filename = './1.wav';
y=getaudiodata(recorder);
audiowrite(filename, y, 8000);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname=sprintf('1.wav');
x=audioread(fname);
sound(x,8000);
axes(handles.axes1);
plot(x);

axes(handles.axes2);
x = double(x);
x = x / max(abs(x));

%��������
FrameLen = 256;
FrameInc = 80;

amp1 = 10;%�����������
amp2 = 1;%������С����
zcr1 = 10;%�����ʴ�����ֵ
zcr2 = 5;%������С����ֵ

maxsilence = 30;  %���������
minlen  = 15;    % ���ʱ������
status  = 0;%0�����Σ�1���ȶΣ�2�����Σ�3������
count   = 0;%�����������ȡ�
silence = 0;%��������

%���������
tmp1  = enframe(x(1:end-1), FrameLen, FrameInc);
tmp2  = enframe(x(2:end)  , FrameLen, FrameInc);
signs = (tmp1.*tmp2)<0;
diffs = (tmp1 -tmp2)>0.02;
zcr   = sum(signs.*diffs, 2);

%�����ʱ����
amp = sum(abs(enframe(filter([1 -0.9375], 1, x), FrameLen, FrameInc)), 2);

%������������
amp1 = min(amp1, max(amp)/4);
amp2 = min(amp2, max(amp)/8);

%��ʼ�˵���
x1 = 0; 
x2 = 0;
for n=1:length(zcr)
   goto = 0;
   switch status
   case {0,1}                   % 0 = ����, 1 = ���ܿ�ʼ
      if amp(n) > amp1          % ȷ�Ž���������
         x1 = max(n-count-1,1);%���ñ��
         status  = 2;
         silence = 0;
         count   = count + 1;
      elseif amp(n) > amp2 | ... % ���ܴ���������
             zcr(n) > zcr2
         status = 1;
         count  = count + 1;
      else                       % ����״̬
         status  = 0;
         count   = 0;
      end
   case 2,                       % 2 = ������
      if amp(n) > amp2 | ...     % ������������
         zcr(n) > zcr2
         count = count + 1;
      else                       % ����������
         silence = silence+1;
         if silence < maxsilence % ����������������δ����
            count  = count + 1;
         elseif count < minlen   % ��������̫�̣���Ϊ������
            status  = 0;
            silence = 0;
            count   = 0;
         else                    % ��������
            status  = 3;%3�Ļ��ǽ����Ρ�
         end
      end
   case 3,
      break;
   end
end   

count = count-silence/2;
x2 = x1 + count -1;

plot(x)
ylabel('Speech');
line([x1*FrameInc x1*FrameInc], [-1 1], 'Color', 'red');
line([x2*FrameInc x2*FrameInc], [-1 1], 'Color', 'red');

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load mylabel.mat;
fname=sprintf('test.wav');
x=audioread(fname);
[x1 x2] = vad(x);
x=0.2*x/max(x);
m=mfcc(x);
m=m(x1-2:x2-2,:);
m=m';
m = vqlbg(m,4);
m = m';
m = reshape(m(1:4,:)',[1,96]);
m = [1 m];
elm_predict(m);
%elm_kernel(dataset, m, 1, 1, 'RBF_kernel',100);
load elm_output.mat;
a=sprintf('%s', label(output,:));
set(handles.edit2,'string',a);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname = sprintf('1.wav');
[y fs]=audioread(fname);
fs=8000;
n=length(y);  %ѡȡ�任�ĵ���
y_p=fft(y,n);
f=fs*(0:n/2-1)/n;
%�Լ���������źŽ���ȥ�����
fp=1500;fc=1700;As=100;Ap=1;

wc=2*pi*fc/fs;wp=2*pi*fp/fs;
wdel=wc-wp;
beta=0.112*(As-8.7);
N=ceil((As-8)/2.285/wdel);
wn=kaiser(N+1,beta);
ws=(wp+wc)/2/pi;
b=fir1(N,ws,wn);

%��ǰδ��ͨ�˲�����ƽ׶� ������Ϊȥ��
x=fftfilt(b,y);
X=fft(x,n);
axes(handles.axes3);
plot(x);
sound(x,8000);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname=sprintf('1.wav');
[wavin_t,fs]=audioread(fname);
a=2;      %��������
b=0.01;     %���油������
c=0;        %c=0ʱ���������������п�����c=1ʱ�����п�������
fs=8000;
wav_length=length(wavin_t);

frame_len=320;step_len=160;
frame_num=ceil((wav_length-step_len)/step_len);
wavin=zeros(1,frame_num*frame_len);
wavin(1:wav_length)=wavin_t(:);
inframe=zeros(frame_len,frame_num);
for i=1:frame_num;
    inframe(:,i)=wavin(((i-1)*step_len+1):((i-1)*step_len+frame_len));
end;
%inframe=(ENFRAME(wavin,frame_len,step_len))';   %��֡
%frame_num=size(inframe,2);          %��֡��
window=hamming(frame_len);          %���庺����

%�ֱ��ÿ֡fft�����ֵ�������-----------------------------------------------
for i=1:frame_num;
    fft_frame(:,i)=fft(window.*inframe(:,i));
    abs_frame(:,i)=abs(fft_frame(:,i));
    ang_frame(:,i)=angle(fft_frame(:,i));
end;

%ÿ������֡ƽ��-------------------------------------------------------------
abs_frame_f=abs_frame;
for i=2:frame_num-1;
    abs_frame_f(:,i)=mean(abs_frame(:,(i-1):(i+1)),2);
end;
abs_frame=abs_frame_f;

%���������-----------------------------------------------------------------
%������ÿһԪ��Ϊ��
%g(k)=(Py(k)-a*Pn(k))/Py(k)
%Py��Pn�ֱ�Ϊ���������������Ĺ����׹��ƣ�����MATLAB���Դ���pmtm����������
%�ɸ�����Ҫ����a�Ĵ�С�����õ����õ�Ч��

%�öര�׷�����ÿһ֡���ݽ��й����׹���
for i=1:frame_num;
    per_PSD(:,i)=pmtm(inframe(:,i),3,frame_len,'twosided');
end;

%�Թ����׵�ÿ������֡����ƽ��
per_PSD_f=per_PSD;
for i=2:frame_num-1;
    per_PSD_f(:,i)=mean(per_PSD(:,(i-1):(i+1)),2);
end;
per_PSD=per_PSD_f;

%ȡǰ20֡��Ϊ����֡��ȡ��ƽ����Ϊ�����Ĺ����׹���
noise_PSD=mean(per_PSD(:,1:20),2);

%���������
for k=1:frame_num;
    g(:,k)=(per_PSD(:,k)-a*noise_PSD)./per_PSD(:,k);
end;

%�����油����ֵ������С�ڸ���ֵ������ϵ��������ֵ�����棬�����ɼ�����������
spec_floor=b*noise_PSD./per_PSD(:,k);
spec_floor=spec_floor(:,ones(1,frame_num));
[I,J]=find(g<spec_floor);
gf=g;
gf(sub2ind(size(gf),I,J))=spec_floor(sub2ind(size(gf),I,J));
if c==0;
    g=gf;
else g=gf.^0.5;
end;

%�׼�----------------------------------------------------------------------
sub_frame=g.*abs_frame;

% % ������֡˥��------------------------------------------------------------
% T=20*log10(mean(sub_frame./(abs_noise*ones(1,frame_num))));
% T_noise=mean(T(:,1:20),2);
% c=10^(-2/3);            %˥��ϵ��Ϊ 10^(-1.5)
% noise_frame=find(T<T_noise);
% sub_frame(:,noise_frame)=c*sub_frame(:,noise_frame);

%�������źŻ�ԭ��ʱ��-------------------------------------------------------
wavout=zeros(1,(frame_num-1)*step_len+frame_len);
j=sqrt(-1);
i=1;
for t=1:step_len:((frame_num-1)*step_len+1);
    wavout(:,t:(t+frame_len-1))=wavout(:,t:(t+frame_len-1))+real(ifft(sub_frame(:,i).*exp(j*ang_frame(:,i))))';
    i=i+1;
end;

sound(wavout,8000);
axes(handles.axes4);
plot(wavout);
audiowrite('test.wav',wavout,fs);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load temp.mat;
y=x;
fs=8000;
n=length(y);  %ѡȡ�任�ĵ���
y_p=fft(y,n);
f=fs*(0:n/2-1)/n;
%�Լ���������źŽ���ȥ�����
fp=1500;fc=1700;As=100;Ap=1;

wc=2*pi*fc/fs;wp=2*pi*fp/fs;
wdel=wc-wp;
beta=0.112*(As-8.7);
N=ceil((As-8)/2.285/wdel);
wn=kaiser(N+1,beta);
ws=(wp+wc)/2/pi;
b=fir1(N,ws,wn);

%��ǰδ��ͨ�˲�����ƽ׶� ������Ϊȥ��
x=fftfilt(b,y);
X=fft(x,n);
axes(handles.axes3);
plot(x);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load temp.mat;

a=2;      %��������
b=0.01;     %���油������
c=0;        %c=0ʱ���������������п�����c=1ʱ�����п�������
wavin_t=x;
fs=8000;
wav_length=length(wavin_t);

frame_len=320;step_len=160;
frame_num=ceil((wav_length-step_len)/step_len);
wavin=zeros(1,frame_num*frame_len);
wavin(1:wav_length)=wavin_t(:);
inframe=zeros(frame_len,frame_num);
for i=1:frame_num;
    inframe(:,i)=wavin(((i-1)*step_len+1):((i-1)*step_len+frame_len));
end;
%inframe=(ENFRAME(wavin,frame_len,step_len))';   %��֡
%frame_num=size(inframe,2);          %��֡��
window=hamming(frame_len);          %���庺����

%�ֱ��ÿ֡fft�����ֵ�������-----------------------------------------------
for i=1:frame_num;
    fft_frame(:,i)=fft(window.*inframe(:,i));
    abs_frame(:,i)=abs(fft_frame(:,i));
    ang_frame(:,i)=angle(fft_frame(:,i));
end;

%ÿ������֡ƽ��-------------------------------------------------------------
abs_frame_f=abs_frame;
for i=2:frame_num-1;
    abs_frame_f(:,i)=mean(abs_frame(:,(i-1):(i+1)),2);
end;
abs_frame=abs_frame_f;

%���������-----------------------------------------------------------------
%������ÿһԪ��Ϊ��
%g(k)=(Py(k)-a*Pn(k))/Py(k)
%Py��Pn�ֱ�Ϊ���������������Ĺ����׹��ƣ�����MATLAB���Դ���pmtm����������
%�ɸ�����Ҫ����a�Ĵ�С�����õ����õ�Ч��

%�öര�׷�����ÿһ֡���ݽ��й����׹���
for i=1:frame_num;
    per_PSD(:,i)=pmtm(inframe(:,i),3,frame_len,'twosided');
end;

%�Թ����׵�ÿ������֡����ƽ��
per_PSD_f=per_PSD;
for i=2:frame_num-1;
    per_PSD_f(:,i)=mean(per_PSD(:,(i-1):(i+1)),2);
end;
per_PSD=per_PSD_f;

%ȡǰ20֡��Ϊ����֡��ȡ��ƽ����Ϊ�����Ĺ����׹���
noise_PSD=mean(per_PSD(:,1:20),2);

%���������
for k=1:frame_num;
    g(:,k)=(per_PSD(:,k)-a*noise_PSD)./per_PSD(:,k);
end;

%�����油����ֵ������С�ڸ���ֵ������ϵ��������ֵ�����棬�����ɼ�����������
spec_floor=b*noise_PSD./per_PSD(:,k);
spec_floor=spec_floor(:,ones(1,frame_num));
[I,J]=find(g<spec_floor);
gf=g;
gf(sub2ind(size(gf),I,J))=spec_floor(sub2ind(size(gf),I,J));
if c==0;
    g=gf;
else g=gf.^0.5;
end;

%�׼�----------------------------------------------------------------------
sub_frame=g.*abs_frame;

% % ������֡˥��------------------------------------------------------------
% T=20*log10(mean(sub_frame./(abs_noise*ones(1,frame_num))));
% T_noise=mean(T(:,1:20),2);
% c=10^(-2/3);            %˥��ϵ��Ϊ 10^(-1.5)
% noise_frame=find(T<T_noise);
% sub_frame(:,noise_frame)=c*sub_frame(:,noise_frame);

%�������źŻ�ԭ��ʱ��-------------------------------------------------------
wavout=zeros(1,(frame_num-1)*step_len+frame_len);
j=sqrt(-1);
i=1;
for t=1:step_len:((frame_num-1)*step_len+1);
    wavout(:,t:(t+frame_len-1))=wavout(:,t:(t+frame_len-1))+real(ifft(sub_frame(:,i).*exp(j*ang_frame(:,i))))';
    i=i+1;
end;

sound(wavout,8000);
axes(handles.axes4);
plot(wavout);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton2.
function pushbutton2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
