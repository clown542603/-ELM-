%�������ܣ�    �˲�  
%ʱ�䣺        2018.5.5
%���ߣ�        ������
%���룺        ��
%�����        ��
clear;clc
[filename,pathname]=uigetfile('�˲�.wav','��ѡ�������ļ���');
[y fs]=audioread([pathname filename]);
%sound(y,fs);
n=length(y);  %ѡȡ�任�ĵ���
y_p=fft(y,n);
f=fs*(0:n/2-1)/n;
figure(1);
subplot(2,1,1);
plot(y);
title('ԭʼ�����źŲ�����ʱ����');
xlabel('ʱ����');
ylabel('��ֵ');
subplot(2,1,2);
plot(f, abs(y_p(1:n/2)));
xlabel('Ƶ��Hz');
ylabel('Ƶ�ʷ�ֵ');

%�Լ���������źŽ���ȥ�����
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

%��ǰδ��ͨ�˲�����ƽ׶� ������Ϊȥ��
x=fftfilt(b,y);
X=fft(x,n);
figure(4);
subplot(2,2,1);plot(abs(y_p));
title('�˲�ǰ�ź�Ƶ��');
subplot(2,2,2);plot(abs(X));
title('�˲����ź�Ƶ��');
subplot(2,2,3);plot(y);
title('�˲�ǰ�źŵĲ���');
subplot(2,2,4);plot(x);
title('�˲����źŵĲ���');
sound(x,fs);

















