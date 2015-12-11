function [freqBands,parameters] = FindSignificantFrequencies(data,parameters)

% SYNTAX: freqBands = findSignificantFrequencies(data,parameters)
%
% Given 'data' (nTrials x nSamples) and parameters (optional), this
% function calculates the frequency ranges in the data that differ
% significantly from noise with a power spectrum that matches the data's
% 1/f^alpha best fit spectrum.
%
% This is done in the following steps:
%
% 1. Compute the spectrum of each trial within a requested frequency
% range*. This uses the Chronux toolbox, available at www.chronux.org.
% Please download it if you don't already have it, and make sure it's in
% your path. You may also ask for Chronux to remove line noise from you
% data prior to computing spectra*.
%
% 2. Normalize the spectra to have identical integrals (optional*)
%
% 3. Peform a robust fit of the spectra to the power function
%    Power(f) = c*f^(-alpha)
%
% 4. For nNoiseIterations
%      4a. Generate nTrials of noise that matches the power spectrum found
%      in step 3.
%      4b. Use either the Maris & Oostenveld Method or Principal Components
%      Method* to define frequency bands in the data that differ from the
%      generated noise.
% 5. Return the frequencies that were shown to be significant in at least
% softIntersectionThreshold* percent of the nNoiseIterations.
%
% * Parameters may be passed in as a structure, or if omitted or left
% empty, the user will be prompted to select parameters. To create a
% parameter structure, please use createFindSigFreqsParameterStruct (and
% read the documentation there). Parameters that may be defined by the user
% include:
% *frequencyRange (2-element vector, [fMin, fMax]; default: [5,55])
% *samplingFrequency (should be specified in Hz; default: 1000)
% *tapers (should be a 3 element array, [W T p], where W is the bandwidth
% (in Hz), T is the duration (in seconds) over which tapers are computed,
% and p is a natural number such that Chronux will use 2TW-p tapers. Note
% that if p==1, this achieves the best estimate of the spectrum but can be
% computationally intensive if T and W are large. For more information
% about choosing tapers, please read the Chronux manual or another
% reference on multi-taper analysis). Default: [3 2 3]
% *removeLineNoise (0 or 1; Default: 0)
% *normalizeSpectra (0 or 1: Default: 1)
% *nNoiseIterations (Default: 20)
% *softIntersectionThreshold (A real number greater than 0 and less than or
% equal to 1; Default 0.65)
% divergFunc (Default: 'ttest');
% divergTh (Default: 2);
% neighborMat (Default: ones(size(cond1,3),size(cond1,3)));
% neighborMatTh (Default: 1);
% nMinClustLen (Default: 10);
% nMaxGapUnify (Default: round(nMinClustLen/3));
% pMax4signif (Default: .5);
% bUnifyPosNeg (Default: false);
% nIter (Default: 1000);
%
% The function is a spinoff of the Maris & Oostenvald (2007) method.
%
% (c) Uri Maoz and Emily Mankin, 2015

%% Get Parameters
if ~exist('parameters','var') || isempty(parameters)
    parameters = createFindSigFreqsParameterStruct;
end

%% Compute the spectrum of each trial within a requested frequency range.
[nTrials,nSamples] = size(data);
% transpose the data so that it is in samples by trials
data = transpose(data);


% create a chronux parameters structure, cp
cp = struct('tapers',parameters.tapers,'Fs',parameters.samplingFrequency,...
    'fpass',parameters.frequencyRange,'trialave',0);

if parameters.removeLineNoise
    T = parameters.tapers(2);
    data=rmlinesmovingwinc(data,[T T/2],10,cp,[],...
        parameters.plotLineNoiseRemoval,parameters.lineNoiseFrequencies);
end

[spectra,freq] = mtspectrumc(data,cp);

%
if parameters.plotDataAndSpectra
    if parameters.normalizeSpectra
        nSubplots = 3;
    else
        nSubplots = 2;
    end
    figure('units','normalized','position',[.3,.55,.35,.25]);
    subplot(1,nSubplots,1); plot(data); axis square; title('Data per trial')
    subplot(1,nSubplots,2); plot(freq,spectra); axis square; title('Spectrum per trial')
end


% Normalize the spectra to have identical integrals (optional)

if parameters.normalizeSpectra
    integrals = sum(spectra);
    spectra = spectra./repmat(integrals,size(spectra,1),1);
    if parameters.plotDataAndSpectra
        subplot(1,nSubplots,3); plot(freq,spectra); axis square; title('Normalized spectrum per trial')
    end
end

spectra = transpose(spectra); % nTrials by nFrequencies
data = transpose(data); %nTrials by nSamples

%% Peform a robust fit of the spectra to the power function: Power(f) = c*f^(-alpha)

logBase=2;
[c,alpha,logF,logFreq] = getFit(freq,spectra,logBase);


if parameters.plotFitOfSpectra
    figure;
    nSubplots = 2;
    subplot(1,nSubplots,1); plot(logFreq',logF'); axis square; title('Log-Log Spectra with Best Fit')
    hold on; plot(logFreq(1,:),mean(logF),'r','linewidth',3);
    plot(logFreq(1,[1 end]),c+logFreq(1,[1 end])*(-alpha),'k','linewidth',3);
    xlabel('Log Frequency'); ylabel('Log Power')
    
    subplot(1,nSubplots,2); plot(freq,spectra); axis square; title('Spectra with Best Fit')
    hold on; plot(freq,mean(spectra),'r','linewidth',3);
    plot(freq,(logBase^c)*freq.^(-alpha),'k','linewidth',3);
    xlabel('Frequency'); ylabel('Power')
end



%% For nNoiseIterations
if parameters.plotNoiseSpectra
    f = figure; a = floor(sqrt(parameters.nNoiseIterations));
    b = ceil(sqrt(parameters.nNoiseIterations));
    while a*b<parameters.nNoiseIterations
        b = b+1;
    end
end
%%
for n = parameters.nNoiseIterations:-1:1
    % a. Generate nTrials of noise that matches the power spectrum found in step 3.
    nz =  ColoredNoise(alpha,nSamples,nTrials) * logBase^c;
    noiseSpectra = mtspectrumc(nz,cp);
    
    % scale the noise spectra to have power that matches the real data
    [c2,alpha2,logF,logFreq] = getFit(freq,noiseSpectra',logBase);
    noiseSpectra = noiseSpectra*logBase^(c-c2);
%     if parameters.normalizeSpectra
%         integrals = sum(noiseSpectra);
%         noiseSpectra = noiseSpectra./repmat(integrals,size(noiseSpectra,1),1);
%     else
%         integrals1 = sum(spectra');
%         integrals2 = sum(noiseSpectra);
%         scaleFactor = integrals1./integrals2;
%         noiseSpectra = noiseSpectra.*repmat(scaleFactor,size(noiseSpectra,1),1);
%         noiseSpectra = noiseSpectra';
%     end
    if parameters.plotNoiseSpectra
        ax = subplot(a,b,n,'parent',f);
        plot(freq,noiseSpectra); ch1 = get(ax,'children');
        hold on; plot(freq, spectra); ch2 = get(ax,'children');
        arrayfun(@(x)set(x,'color',(rand(1))*[0 0 1]),ch2)
        arrayfun(@(x)set(x,'color',(rand(1))*[1 0 0]),ch1)
        plot(freq,mean(noiseSpectra,2),'r','linewidth',3)
        plot(freq,mean(spectra),'c','linewidth',3);
    end
    
    noiseSpectra = noiseSpectra'; % should now be nTrials by nFreqs
   
    %%
%     % b. Use  the Maris & Oostenveld Method to define frequency bands in the
%     % data that differ from the generated noise.
% DiffCondsSignif() expects (#trials x #channels x #timestamps (x
% #frequencies  <optional>). So make 'spectra' 3D with #channels=1
    spectra3D=reshape(spectra,size(spectra,1),1,size(spectra,2));
    noiseSpectra3D=reshape(noiseSpectra,...
        size(noiseSpectra,1),1,size(noiseSpectra,2));
    [clustMask{n}, pVals{n}, clustMaskSgnf{n}, pValsSignif{n}] = ...
        DiffCondsSignif (spectra3D, noiseSpectra3D, parameters);
end

%% Return the frequencies that were shown to be significant in at least softIntersectionThreshold percent of the nNoiseIterations.
end

function [c,alpha,logF,logFreq] = getFit(freq,spectra,logBase)

logFreq=log(freq)/log(logBase); logFreq = repmat(logFreq,size(spectra,1),1);
logF=log(spectra)/log(logBase);
% % warnState=warning('off','stats:statrobustfit:IterationLimit');
% This approach should be replaced by non-linear power fitting using
% the fit() function @@@
b=robustfit(logFreq(:),logF(:));
% % warning(warnState);

c = b(1);
alpha = min(-b(2),2);
end