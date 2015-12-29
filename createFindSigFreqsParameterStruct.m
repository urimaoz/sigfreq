function parameters = createFindSigFreqsParameterStruct ...
    (useDefaults, varargin)
% function parameters = createFindSigFreqsParameterStruct ...
%     (useDefaults, varargin)
% 
% This function creates a parameter struct that is used throughout the
% locating place field code. If this function is called without input
% arguments, it will interactively ask questions of the user and set
% parameters appropriately. If default values are generally acceptable to
% the user, however, useDefaults can be entered as true, and then any value
% that should be overwritten may be passed in as a parameter-value pair.
% For example, to change only the frequency range and samplingFrequency,
% you could enter: 
% p = createFindSigFreqsParameterStruct(1,'frequencyRange',[10,30],
%   'samplingFrequency',2000); 


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
% *methodForCalculatingSignificance ('MarisOostenveld' or 'PCA': Default: MarisOostenveld)
% *softIntersectionThreshold (A real number greater than 0 and less than or
% equal to 1; Default 0.65)
% * If using Maris/Oostenveld, you may also define:
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
% (c) Emily Mankin, 2015

p.frequencyRange = [5,55]; % What frequency range, in Hz, should the spectra be calculted over? [fMin, fMax]
p.samplingFrequency = 1000; % What was the sampling frequency of the data (in Hz)?
p.tapers = [3 2 3]; %3 element array of taper parameters [bandWidth,duration,p] where the number of tapers computed will be 2TW-p
p.removeLineNoise = 0; % Should Chronux attempt to remove line noise (e.g. 60 Hz or 50 Hz electrical noise) prior to computing spectra
p.lineNoiseFrequencies = [60 120 180]; % At what frequencies should line noise be removed?
p.normalizeSpectra = 1; % Should spectra be normalized to have the same total power so that each trial contributes equal weight to the fit
p.nNoiseIterations = 20; % How many times should new noise be generated for comparing clusters to noise
p.methodForCalculatingSignificance = 'MarisOostenveld'; %'MarisOostenveld' or 'PCA'
p.softIntersectionThreshold = 0.65; % Out of the noise iterations, how often must a freqeuncy band be included to count as significant in the final output?
p.divergFunc = @GetTStat2Robust; % What method should be used to calculate the difference between data and noise at each frequency?
p.divergTh = 2; % threshold on divergence above which a value belongs to cluster
p.neighborMat = []; % What points in the data count as neighbors/How do you define distance between data points?
p.neighborMatTh = .5; % What is the threshold on the distance function to count as a neighbor?
p.nMinClustLen = 10; % What is the minimum number of continguous points of significance to count as a cluster worth examining?
p.nMaxGapUnify = round(p.nMinClustLen/3); % At what point is the gap between two clusters small enough that they should be combined into a single cluster?
p.Max4signif = 0.5; % I don't actually know what this is for. I'd suggest using the default. But let's ask Uri!
p.bUnifyPosNeg = false; % If two clusters are adjacent and significant but have differing signs on their divergence function, should they be joined?
p.nIter = 1000; % How many iterations of bootstrapping should be performed to identify significant clusters?
p.bNonlinearRegression = true; % Use non-linear regression for the fit (true) or linear regression in log space (false), i.e., log-log fitting
p.nSmoothingSpan = 0.1; % Smoothing span for the spectrum. If < 1 signifies ratio of data to use, if >1 signifies #samples to use
p.topPercentile = 50; % What percent of trials should be used at a time in the iterative approach to find bands? Must be between 5 and 100. 100 indicates no iteration.
p.plotDataAndSpectra = 1;
p.plotFitOfSpectra = 1;
p.plotNoiseSpectra = 1;
p.plotLineNoiseRemoval = 0;

if nargin==0||~useDefaults
    cprintf('*string','You may hit ''Enter'' to accept the default value for any parameter.\n')
    p = updateP(p,'frequencyRange','What frequency range, in Hz, should the spectra be calculted over? [fMin, fMax]',@(x)length(x)==2 && x(1)<x(2) && x(1)>=0);
    p = updateP(p,'samplingFrequency','What was the sampling frequency of the data (in Hz)?',@(x)length(x)==1 && x>0);
    p = updateP(p,'tapers','3 element array of taper parameters [bandWidth,duration,p] where the number of tapers computed will be 2TW-p',...
        @(x)length(x)==3 && all(x>0) && x(3)>=1 && x(3)==round(x(3)));
    p = updateP(p,'removeLineNoise','Should Chronux attempt to remove line noise (e.g. 60 Hz or 50 Hz electrical noise) prior to computing spectra');
    if p.removeLineNoise
        p = updateP(p,'lineNoiseFrequencies','At what frequencies should line noise be removed?');
        p = updateP(p,'plotLineNoiseRemoval','Would you like to see the pre/post line noise spectra?');
    end
    p = updateP(p,'normalizeSpectra','Should spectra be normalized to have the same total power so that each trial contributes equal weight to the fit');
    p = updateP(p,'nNoiseIterations','How many times should new noise be generated for comparing clusters to noise',@(x)length(x)==1);
    p = updateP(p,'softIntersectionThreshold','Out of the noise iterations, how often must a freqeuncy band be included to count as significant in the final output?',...
        @(x)length(x)==1 && x>0 && x<=1);
    p = updateP(p,'divergFunc','What method should be used to calculate the difference between data and noise at each frequency? (Please enter as a function handle)',...
        @(x)isa(x,'function_handle'));
    p = updateP(p,'divergTh','Threshold on divergence above which a value belongs to cluster');
    p = updateP(p,'neighborMat','What points in the data count as neighbors/How do you define distance between data points? If left empty, all points will be considered neighbors');
    p = updateP(p,'neighborMatTh','What is the threshold on the distance function to count as a neighbor?');
    p = updateP(p,'nMinClustLen','What is the minimum number of continguous points of significance to count as a cluster worth examining?');
    p = updateP(p,'nMaxGapUnify','At what point is the gap between two clusters small enough that they should be combined into a single cluster?');
    p = updateP(p,'Max4signif','I don''t actually know what this is for. I''d suggest using the default. But let''s ask Uri!');
    p = updateP(p,'bUnifyPosNeg','If two clusters are adjacent and significant but have differing signs on their divergence function, should they be joined?');
    p = updateP(p,'nIter','How many iterations of bootstrapping should be performed to identify significant clusters?');
    p = updateP(p,'bNonlinearRegression ','Use nonlinear regression (true) or linear regression in log space (false) to calculate alpha for 1/(f^alpha)?');
    p = updateP(p,'nSmoothingSpan ','What smoothing span to use for spectrum? Value < 1 signifies ratio of data to use, >1 signifies #samples to use');
    p = updateP(p,'topPercentile','What percent of trials should be used at a time in the iterative approach to find bands? Must be between 5 and 100. 100 indicates no iteration.',@(x)x>=5 & x <=100);
end

p = parseParameters(p,varargin{:});
parameters = p;

% ========================== Local subfunctions ===========================
function p = updateP(p,field,question,assertionFunction)
% Ask user to update field 
disp(field)
spaces = regexp(question,'\s');
cutoffs = [0 spaces(logical(diff(ceil(spaces/75)))) length(question)];
for i=1:length(cutoffs)-1
    disp(['     ',question(cutoffs(i)+1:cutoffs(i+1))])
end
if isnumeric(p.(field))
    answer = input(['Default: ',num2str(p.(field)),'\n']);
elseif isa(p.(field),'function_handle')
    answer = input(['Default: ',func2str(p.(field)),'\n']);
else
    answer = input(['Default: ',p.(field),'\n']);
end

if ~isempty(answer) && exist('assertionFunction','var')&&~isempty(assertionFunction)
    tf = assertionFunction(answer);
    if ~tf
        cprintf('*string','That answer is not acceptable, please try again.\n')
        p = updateP(p,field,question,assertionFunction);
    else
        p.(field)=answer;
    end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function p = parseParameters(p,varargin)
% SYNTAX: p = parseParameters(p,varargin)
% Takes a struct, p, with default values and overwrites any values that are
% passed in as part of parameter-value pairs inside varargin. Parameter
% names passed in must be in the struct p already, which prevents users
% from mistyping an input argument name and not realizing that their new
% value did not overwrite a default value.
%
% Emily Mankin

assert(mod(length(varargin),2)==0,'Variable arguments must be passed in as pairs')

for i=1:2:length(varargin)
    if ~ischar(varargin{i})
        error('Optional Arguments must be of the form ''parameterName'', parameterValue');
    elseif isfield(p,varargin{i})
        eval(['p.',varargin{i},' = varargin{i+1};'])
    else
        error(['"', varargin{i},'"',' is not a valid input parameter.']);
    end
end