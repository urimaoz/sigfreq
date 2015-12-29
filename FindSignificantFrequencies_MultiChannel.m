function [freqBands,parameters] = FindSignificantFrequencies_MultiChannel(data,parameters);
% Expects data in the form nTrials x nSamples x nChannels.
% Runs through each channel and returns the frequency bands with power
% greater than expected by a 1/f^alpha noise spectrum. See documentation
% inside the function FindSignificantFrequencies for more details,
% including how to enter a parameter structure. Or omit paramters and you
% will be walked through the creation of a parameter struct as part of the
% function call.

if ~exist('parameters','var')||isempty(parameters)
    parameters = createFindSigFreqsParameterStruct;
end
parameters.plotDataAndSpectra = 0;
parameters.plotFitOfSpectra = 0;
parameters.plotNoiseSpectra = 0;
parameters.plotLineNoiseRemoval = 0;

[nTrials,nSamples,nChannels] = size(data);
freqBands = cell(1,nChannels);

rmi(-1);
for ch = 1:nChannels
    freqBands{ch} = findSignificantFrequencies(data(:,:,ch),parameters);
    rmi(sprintf(';  finished step %d of %d',ch,nChannels));
end
%%
f = parameters.frequencyRange(1):.01:parameters.frequencyRange(2);
isSig = zeros(nChannels,length(f));
for ch = 1:nChannels
    for fb = 1:size(freqBands{ch},1)
        isSig(ch,f>=freqBands{ch}(fb,1) & f<=freqBands{ch}(fb,2)) = 1;
    end
end
figure; imagesc(f,1:nChannels,isSig)
