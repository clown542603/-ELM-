%�������ܣ�  ��ʮ������������ң�ǰ�Ÿ���Ϊѵ�����������һ������Ϊ����������������
%ʱ�䣺      2018.4.20
%���ߣ�      ������
%���룺      ��
%�����      ��

function Randomdisturbance()
clear;clc;
load alldata.mat;

%�������������
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