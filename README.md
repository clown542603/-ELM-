# -ELM-

由于提取mfcc特征参数受录音时间长度和端点检测影响，导致特征参数维度不太一致，
本文采用矢量量化对mfcc特征参数进行规整（降维），然后使用十折交叉验证，在无噪声
的录音文件中识别能力较好，有93左右，但是在实时录音时识别率较差，该模型泛化能力
还是太差了，前面加了一段语音增强，还是不太理想，后期打算提取用EMD优化的mfcc特征参数，
还有采用动态时间规整（dtw）看实时录音的情况，看泛化能力会不会有所增强。