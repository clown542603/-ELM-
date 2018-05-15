%���ܣ�   ��ʱ��������ʶ��Ч��
%ʱ�䣺   2018.4.20
%���ߣ�   ������
%���룺   ��
%�����   ��
clc;clear;
load mylabel.mat;
load alldata.mat;
load traindata.mat;
load predictdata.mat;

%¼��
recorder = audiorecorder(8000,16,1);
disp('Start speaking.')
recordblocking(recorder, 3);
disp('End of Recording.');
filename = './1.wav';
y=getaudiodata(recorder);
audiowrite(filename, y, 8000);

fname=sprintf('./trainning/%d.wav', 1);
disp(fname);
x=audioread(fname);
[x1 x2] = vad(x);
x=0.2*x/max(x);
m=mfcc(x);
m=m(x1-2:x2-2,:);
a=m;
m=m';
m = vqlbg(m,4);
m = m';
subplot(2,1,1);
plot(a);
m = reshape(m(1:4,:)',[1,96]);
m = [1 m];

elm_train(dataset, 1, 1000, 'sig');
elm_predict(m);
%elm_kernel(dataset, m, 1, 1, 'RBF_kernel',100);
load elm_output.mat;
fprintf('ʶ��Ϊ%s\n', label(output,:));
%��ʬ���� ����ѿ� ������һ ��ˮһս ���ٳ��� ͵�컻�� �Ի����� ���ɹ��� ׳־���� ��������

