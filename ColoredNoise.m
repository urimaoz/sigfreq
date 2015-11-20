function res = ColoredNoise (alpha,nSamples,nChannels)
% function res = ColoredNoise (alpha,nSamples,nChannels)
% 
% Generate colored noise of the form 1/f^'alpha': 'nSamples' over
% 'nChannels' 
hColoredNoise=dsp.ColoredNoise(alpha,nSamples,nChannels); 
res=step(hColoredNoise); 
end