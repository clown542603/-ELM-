%�������ܣ�    ͳ��10��ƽ��ʶ��׼ȷ��
%ʱ�䣺        2018.5.5
%���ߣ�        ������
%���룺        ��
%�����        ��

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
ylabel('ѵ����Ԥ����ȷ�� ')
title('ѵ������ȷ��')
subplot(2,1,2);
plot(1:10,R(:,2),'b:o')
xlabel('x')
ylabel('Ԥ�⼯Ԥ����ȷ�� ')
title('Ԥ�⼯��ȷ��')
accuracy = sum/10;

str = sprintf('ʮ�β��Ե�ƽ��ʶ����Ϊ%d%%', accuracy*100);
disp(str);