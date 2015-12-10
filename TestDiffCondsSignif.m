rng default; 
addpath ../../Replicablility_Majed/Replicable-Difference/
addpath ../../SignifFreqExtract_Emily/

alpha=-1; nSamples=1000; nTrials=50; nChannels=10;
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
sDataFile='TestDiffCondsSignifData';
if (exist(sDataFile,'file'))
    load(sDataFile);
else
    t=linspace(0,1,nSamples);
    rmi(-1);
    for c=1:nChannels
        for tr=1:nTrials
            x=0;
            for k=10:.01:20
                x=x+sin(k*2*pi*t+rand(1)*2*pi-pi);
            end
            cond2(tr,c,:)=cond2(tr,c,:)+reshape(x,1,1,[]);
            rmi(sprintf('; c=%d, tr=%d',c,tr));
        end
    end
    save(sDataFile,'cond2');
end

params=struct(...
    'divergFunc',@dPrime,...
    'divergTh',1,...
    'neighborMat',ones(nChannels,nChannels),...
    'neighborMatTh',2,...
    'nMinClustLen',0,...
    'nMaxGapUnify',5,...
    'pMax4signif',0.05,...
    'bUnifyPosNeg',false,...
    'nIter',100);
DiffCondsSignif(cond1,cond2,params);