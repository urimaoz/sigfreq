function [freqBands,parameters] = FindSignificantFrequencies(data,parameters)

% SYNTAX: freqBands = FindSignificantFrequencies(data,parameters)
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


% create a chronux parameters structure, chronuxParams
chronuxParams = struct('tapers',parameters.tapers,'Fs',parameters.samplingFrequency,...
    'fpass',parameters.frequencyRange,'trialave',0);

if (parameters.removeLineNoise)
    T = parameters.tapers(2);
    data=rmlinesmovingwinc(data,[T T/2],10,chronuxParams,[],...
        parameters.plotLineNoiseRemoval,parameters.lineNoiseFrequencies);
end

[spectra,freq] = mtspectrumc(data,chronuxParams);

% Smooth the spectrum
if (parameters.plotSpectraAndSmoothed)
    spectraSmooth=SmoothSpectrum(spectra,parameters);
%    [[We need to add plotting code here. But how do we decide which trials
%    to plot? Do we get it as a parameter from the user? Do we plot the
%    first 10 or so trials? Do we plot all the trials (risking having each
%    subolot be very small)?]]
    spectra=spectraSmooth;
else
    spectra=SmoothSpectrum(spectra,parameters);
end

% Plot data and spectra (optional)
if (parameters.plotDataAndSpectra)
    if (parameters.normalizeSpectra)
        nSubplots = 3;
    else
        nSubplots = 2;
    end
    if (strcmpi(get(0,'DefaultFigureWindowStyle'),'normal'))
        figure('units','normalized','position',[.3,.55,.35,.25]);
    else
        figure('units','normalized');
    end
    subplot(1,nSubplots,1); plot(data); axis square; title('Data per trial')
    subplot(1,nSubplots,2); plot(freq,spectra); axis square; title('Spectrum per trial')
end


% Normalize the spectra to have identical integrals (optional)
if (parameters.normalizeSpectra)
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
[c,alpha,logF,logFreq] = getFit(freq,spectra,logBase,...
    parameters.bNonlinearRegression,1);


if (parameters.plotFitOfSpectra)
    figure;
    if (~parameters.bNonlinearRegression)
        nSubplots = 2;
        subplot(1,nSubplots,1);
        plot(logFreq',logF'); axis square; title('Log-Log Spectra with Best Fit')
        hold on; plot(logFreq(1,:),mean(logF),'r','linewidth',3);
        plot(logFreq(1,[1 end]),c+logFreq(1,[1 end])*(-alpha),'k','linewidth',3);
        xlabel('Log Frequency'); ylabel('Log Power')
        
        subplot(1,nSubplots,2);
    end
    
    plot(freq,spectra); axis square; title('Spectra with Best Fit')
    hold on; plot(freq,mean(spectra),'r','linewidth',3);
    if (parameters.bNonlinearRegression)
        plot(freq,c*freq.^(-alpha),'k','linewidth',3);
    else
        plot(freq,(logBase^c)*freq.^(-alpha),'k','linewidth',3);
    end
    xlabel('Frequency'); ylabel('Power')
end

% Get subplot dimensions
if (parameters.plotNoiseSpectra)
    a = floor(sqrt(parameters.nNoiseIterations));
    b = ceil(sqrt(parameters.nNoiseIterations));
    while (a*b<parameters.nNoiseIterations)
        b = b+1;
    end
else
    a = []; b = [];
end

%% Iteratively Choose a Subset of Data on which to Test For Significant Bands

if parameters.topPercentile == 100
    % If using all data, no need to do this iteratively. Run once and stop
    freqBands = getFreqBandsFromData(data,spectra,parameters,a,b,c,chronuxParams);
    keepGoing = 0;
else
    keepGoing = 1;
    freqBands = [];
%     runAtleastOnce = 0;
    if parameters.bNonlinearRegression
        bestFitLine = c*freq.^(-alpha);
        deviationBetween = mean(spectra) - bestFitLine;
    else
        bestLogFitLine = c+logFreq(1,:)*(-alpha);
        meanLogSpectra = mean(logF);
        deviationBetween = meanLogSpectra - bestLogFitLine;
    end
    
    
    
end

while keepGoing
    % Find the frequency where the mean deviation from the best-fit line is
    % maximal (currently, not in absolute value, since we're only looking
    % for positive bumps)
    
    [m,i] = nanmax(deviationBetween);
    
    
    % if we are at a point where the maximal devation between the untested
    % mean spectra and the best fit line is less than zero, then all bumps
    % have been tested already. Stop.
    if  m<0
        keepGoing = 0;
    else
        % for each trial, consider the value of that trial's spectrum at
        % the frequency that most deviates from the mean. Take the trials
        % that represent the top X% of these values, where X is defined by
        % parameters.topPercentile
        thisThreshold = prctile(spectra(:,i),100-parameters.topPercentile);
        trialsToKeep = spectra(:,i)>=thisThreshold;
        theseData = data(trialsToKeep,:);
        theseSpectra = spectra(trialsToKeep,:);
        theseFreqBands = getFreqBandsFromData(theseData,theseSpectra,parameters,a,b,c,chronuxParams);
        freqBands = [freqBands;theseFreqBands];
%         runAtleastOnce = 1;
        sigInds = false(size(freq));
        for tfb = 1:size(theseFreqBands)
            sigInds(freq>=theseFreqBands(tfb,1) & freq<=theseFreqBands(tfb,2)) = 1;
        end
        if ~sigInds(i)
            %then the pre-determined likeliest to be significant frequency
            %isn't significant, so stop;
            keepGoing = 0;
        else
            % eliminate the already significant frequency bands from
            % consideration for possible best frequency range and then
            % continue
            deviationBetween(sigInds) = NaN;
        end
        
    end
end

freqBands = unifyFreqBands(freqBands);



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function freqBands = getFreqBandsFromData(data,spectra,parameters,a,b,c,cp)
        [nTrials,nSamples] = size(data);
        if (parameters.plotNoiseSpectra)
            f = figure;
        end
        for n = parameters.nNoiseIterations:-1:1
            % a. Generate nTrials of noise that matches the power spectrum found in step 3.
            nz =  ColoredNoise(min(alpha,2),nSamples,nTrials) * logBase^c;
            nzAmps = sum(abs(nz));
            dataAmps = sum(abs(data'));
            scaleFactor = repmat(dataAmps./nzAmps,nSamples,1);
            nz = nz.*scaleFactor;
            noiseSpectra = mtspectrumc(nz,cp);
            % Smooth noise spectra same as data spectra
            noiseSpectra = SmoothSpectrum(noiseSpectra,parameters);
            % scale the noise spectra to have power that matches the real data
            [~,c,logF,logFreq] = getFit(freq,noiseSpectra',logBase,...
                parameters.bNonlinearRegression,0);
%             noiseSpectra = noiseSpectra*logBase^(c-c2);
            if (parameters.normalizeSpectra)
                integrals = sum(noiseSpectra);
                noiseSpectra = noiseSpectra./repmat(integrals,size(noiseSpectra,1),1);
            end
            if (parameters.plotNoiseSpectra)
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
            [~, ~, clustMaskSgnf{n}, ~] = ...
                DiffCondsSignif (spectra3D, noiseSpectra3D, parameters);
        end
        
        %% Return the frequencies that were shown to be significant in at least softIntersectionThreshold percent of the nNoiseIterations.
        sigClusts = cell2mat(clustMaskSgnf');
        %         sigClusts = sum(sigClusts~=0)/parameters.nNoiseIterations>parameters.softIntersectionThreshold;
        % changed to >0 because for now we are insisting on positive bumps.
        sigClusts = sum(sigClusts>0)/parameters.nNoiseIterations>parameters.softIntersectionThreshold;
        freqBands = freq(continuousRunsOfTrue(sigClusts));
        freqBands(freqBands(:,1)==freqBands(:,2),:)=[];
        
    end

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [c,alpha,logF,logFreq] = getFit(freq,spectra,logBase,...
    bNonlinearFit,dataNotNoise)

logFreq=log(freq)/log(logBase); logFreq = repmat(logFreq,size(spectra,1),1);
logF=log(spectra)/log(logBase);

if (bNonlinearFit)
    if dataNotNoise
        options=fitoptions('power1','Robust','On');
    else
        options = fitoptions('power1','Robust','Off');
    end
    freq=repmat(freq,size(spectra,1),1);
    fitObj=fit(freq(:),spectra(:),'power1',options);
    c = fitObj.a;
    alpha = -fitObj.b;
else
    % % warnState=warning('off','stats:statrobustfit:IterationLimit');
    B=robustfit(logFreq(:),logF(:));
    % % warning(warnState);
    c = B(1);
    alpha = -B(2);
end


if dataNotNoise
    fprintf('Data; alpha = %0.3g\n',alpha)
else
    fprintf('Noise; alpha = %0.3g\n',alpha)
end
% if alpha>2
% %     if giveWarning
% %     warning(['Alpha was computed as ',num2str(alpha),'. If alpha is far from 2, the fitting will not be appropriate. Please check carefully']);
% %     end
%     alpha = 2;
% end
end



function freqBands = unifyFreqBands(freqBands)
if size(freqBands,1)>1
    [~,i] = sort(freqBands(:,1));
    freqBands = freqBands(i,:);
    counter = 1;
    while counter<size(freqBands,1)
        if freqBands(counter+1,1)>freqBands(counter,2)
            counter = counter+1;
        else
            toUnify = freqBands(:,1)<=freqBands(counter,2);
            maxVal = max(freqBands(toUnify,2));
            freqBands(counter,2) = maxVal;
            toUnify(1:counter) = 0;
            freqBands(toUnify,:) = [];
        end
    end
end

end