%函数功能：    滤波  
%时间：        2018.5.5
%作者：        吴宁旭
%输入：        无
%输出：        无
clear;clc
[filename,pathname]=uigetfile('滤波.wav','请选择语音文件：');
[y fs]=audioread([pathname filename]);
%sound(y,fs);
n=length(y);  %选取变换的点数
y_p=fft(y,n);
f=fs*(0:n/2-1)/n;
figure(1);
subplot(2,1,1);
plot(y);
title('原始语音信号采样后时域波形');
xlabel('时间轴');
ylabel('幅值');
subplot(2,1,2);
plot(f, abs(y_p(1:n/2)));
xlabel('频率Hz');
ylabel('频率幅值');

%对加噪的语音信号进行去噪程序
fp=1500;fc=1700;As=100;Ap=1;

wc=2*pi*fc/fs;wp=2*pi*fp/fs;
wdel=wc-wp;
beta=0.112*(As-8.7);
N=ceil((As-8)/2.285/wdel);
wn=kaiser(N+1,beta);
ws=(wp+wc)/2/pi;
b=fir1(N,ws,wn);
figure(3);
freqz(b,1);

%此前未低通滤波器设计阶段 接下来为去噪
x=fftfilt(b,y);
X=fft(x,n);
figure(4);
subplot(2,2,1);plot(abs(y_p));
title('滤波前信号频谱');
subplot(2,2,2);plot(abs(X));
title('滤波后信号频谱');
subplot(2,2,3);plot(y);
title('滤波前信号的波形');
subplot(2,2,4);plot(x);
title('滤波后信号的波形');
sound(x,fs);

















