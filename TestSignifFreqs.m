% TestSignifFreq
rng('default'); % For debugging purposes

t=0:.001:1; 
a=1.25;
x=0;
for k=1:.01:100
    x=x+sin(k*2*pi*t+rand(1)*2*pi-pi)/k^a;
end
% We should now have some version of 1/f^a noise in x.

% Create stronger signals at 10 & 20 Hz
x=x+50*sin(10*2*pi*t)/10^a+50*sin(20*2*pi*t)/20^a;
% Draw power of signal vs. frequency
[f,m]=MyFFT(x,t,[]);
% run 1/f^a noise 
[freqRanges,pVals,tValsSum]=FindSignifFreqsFFT(x,t,[3 30],.1,2,1:2);


