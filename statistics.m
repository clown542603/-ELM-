%函数功能：    统计10次平均识别准确率
%时间：        2018.5.5
%作者：        吴宁旭
%输入：        无
%输出：        无

function accuracy = statistics()
clear;clc;
load traindata.mat;
load predictdata.mat;

sum = 0;
R = [];
for i=1:10
    Randomdisturbance();
    [TrainingTime TrainingAccuracy] = elm_train(trainset, 1, 800, 'sig');
    [TestingTime TestingAccuracy] = elm_predict(predictset);
    sum = sum + TestingAccuracy;
    R = [R; TrainingAccuracy TestingAccuracy];
end

subplot(2,1,1);
plot(1:10,R(:,1),'b:o')
xlabel('x')
ylabel('训练集预测正确率 ')
title('训练集正确率')
subplot(2,1,2);
plot(1:10,R(:,2),'b:o')
xlabel('x')
ylabel('预测集预测正确率 ')
title('预测集正确率')
accuracy = sum/10;

str = sprintf('十次测试的平均识别率为%d%%', accuracy*100);
disp(str);