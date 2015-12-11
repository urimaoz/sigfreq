rng default; 
addpath ../../Replicablility_Majed/Replicable-Difference/
addpath ../../SignifFreqExtract_Emily/
IsChronuxInpath;

alpha=-1; nSamples=1000; nTrials=50; nChannels=1;
% cond1=rand(nTrials,nChannels,nSamples);
% cond2=rand(nTrials,nChannels,nSamples);
% cond2(:,:,450:550)=cond2(:,:,450:550)+2;
cond1=permute(reshape(...
    ColoredNoise(-alpha,nSamples,nTrials*nChannels),...
    [nSamples,nTrials,nChannels]),[2 3 1]); 
cond2=permute(reshape(...
    ColoredNoise(-alpha,nSamples,nTrials*nChannels),...
    [nSamples,nTrials,nChannels]),[2 3 1]); 

% We want to add a bump in the 10-20 Hz range to condition 2
sDataFile='TestDiffCondsSignifData.mat';
if (exist(sDataFile,'file'))
    load(sDataFile);
else
    t=linspace(0,1,nSamples);
    rmi(-1);
    for c=1:nChannels
        for tr=1:nTrials
            x=0;
            for k=20:.01:25
                x=x+sin(k*2*pi*t+rand(1)*2*pi-pi);
            end
            cond2(tr,c,:)=cond2(tr,c,:)+reshape(x,1,1,[])./28;
            rmi(sprintf('; c=%d, tr=%d',c,tr));
        end
    end
    save(sDataFile,'cond2');
end

params=createFindSigFreqsParameterStruct(1,...
    'frequencyRange',[5 120],...
    'divergFunc',@dPrime,...
    'divergTh',2,...
    'neighborMat',ones(nChannels,nChannels),...
    'neighborMatTh',2,...
    'nMinClustLen',0,...
    'nMaxGapUnify',5,...
    'bUnifyPosNeg',false,...
    'nIter',100);
FindSignificantFrequencies(squeeze(cond2(:,1,:)),params);