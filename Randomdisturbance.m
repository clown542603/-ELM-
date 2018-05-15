%函数功能：  将十个样本随机打乱，前九个作为训练样本，最后一个保存为测试样本，并保存
%时间：      2018.4.20
%作者：      吴宁旭
%输入：      无
%输出：      无

function Randomdisturbance()
clear;clc;
load alldata.mat;

%将矩阵随机打乱
i=1;
temp = [];
while i < 100
    A = dataset(i:i+9,:);
    rowrank = randperm(size(A,1));
    B=A(rowrank,:);
    temp = [temp; B];
    i=i+10;
end

dataset = temp;

i=1;
trainset = [];
temp = dataset;
while i < 100
    trainset = [trainset; temp(i:i+8, :)];
    i = i+10;
end

i=1;
predictset = [];
temp = dataset;
while i < 100
    predictset = [predictset; temp(i+9, :)];
    i = i+10;
end

save predictdata predictset
save traindata trainset