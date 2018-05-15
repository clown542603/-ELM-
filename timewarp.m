%函数功能：  提取所有训练模板的mfcc特征参数，并利用VQ-LBG算法规整到四帧，然后规整到一帧，添加标号
%时间：      2018.4.20
%作者：      吴宁旭
%输入：      training文件夹下面的音频文件
%输出：      无
function  t=timewarp()
clear; clc;

dataset=[];
for i=1:100
   fname=sprintf('./trainning/%d.wav', i);
   disp(fname);
   x=audioread(fname);
   [x1 x2] = vad(x);
   x=0.2*x/max(x);
   m=mfcc(x);
   m=m(x1-2:x2-2,:);
   m=m';
   m = vqlbg(m,4);
   dataset = [dataset; m'];
end

%save alldata1 dataset

temp=[];
i=1;
while i < 400
    temp = [temp; reshape(dataset(i:i+3,:)', [1,96])];
    i = i + 4;
end

dataset = temp;

lables = [];

for i=1:10
    lables = [lables; i*ones(10,1)];
end

dataset = [lables, dataset];

save alldata dataset