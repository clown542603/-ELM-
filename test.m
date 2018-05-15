%功能：   临时测试语音识别效果
%时间：   2018.4.20
%作者：   吴宁旭
%输入：   无
%输出：   无
clc;clear;
load mylabel.mat;
load alldata.mat;
load traindata.mat;
load predictdata.mat;

%录音
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
fprintf('识别为%s\n', label(output,:));
%行尸走肉 金蝉脱壳 百里挑一 背水一战 兵临城下 偷天换日 卧虎藏龙 八仙过海 壮志凌云 长生不死

