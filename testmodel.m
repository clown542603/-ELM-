%函数功能：    测试隐含节点对识别率的影响
%时间：        2018.4.26
%作者：        吴宁旭
%输入：        无
%输出：        无

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
xlabel('隐含层神经元个数')
ylabel('测试集预测正确率（%） ')
title('隐含层神经元个数对 ELM 性能的影响')




