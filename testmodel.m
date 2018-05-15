%�������ܣ�    ���������ڵ��ʶ���ʵ�Ӱ��
%ʱ�䣺        2018.4.26
%���ߣ�        ������
%���룺        ��
%�����        ��

function trainmodel()
clear;clc;
warning off
load traindata.mat;
load predictdata.mat;

R = [];

for i = 50:50:2000
    [TrainingTime TrainingAccuracy] = elm_train(trainset, 1, i, 'sig');
    [TestingTime TestingAccuracy] = elm_predict(predictset);
    R = [R; TrainingAccuracy TestingAccuracy];
end

figure
plot(50:50:2000,R(:,2),'b:o')
xlabel('��������Ԫ����')
ylabel('���Լ�Ԥ����ȷ�ʣ�%�� ')
title('��������Ԫ������ ELM ���ܵ�Ӱ��')




