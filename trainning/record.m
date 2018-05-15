recorder = audiorecorder(8000,16,1);
disp('Start speaking.')
recordblocking(recorder, 4);
disp('End of Recording.');
filename = './100.wav';
y=getaudiodata(recorder);
audiowrite(filename, y, 8000);