function [clustMask, pVals, clustMaskSgnf, pValsSignif] = ...
    DiffCondsSignif (cond1, cond2, params)
% function [clustMask, pVals, clustMaskSgnf, pValsSignif] = ...
%     DiffCondsSignif (cond1, cond2, params)
% 
% The function computes clusters on differences between 'cond1' and 'cond2'
% being at least 'params.divergTh' apart in a t test.
% Inputs:
% (1&2) matrices of values for conditions 1 and 2, 'cond1' and 'cond2'
% (#trials or #subjects x #channels x #timestamps (x #frequencies; 
% <optional>)), 
% (3) a parameters structure, params, consisting of:
% * divergFunc: Some divergence function between the conditions. The
% function returns a higher value the more divergent cond1&2 are. The sign
% typically indicates which of the conditions is larger, if relevant.
% divergFunc(cond1,cond2) operates on the first dimention of cond1&2.
% * divergTh: threshold on divergence above which a value belongs to cluster
% * d: distance matrix between channels
% * dTh: threshold on distance matrix for channels to be neighbors 
% * nMinClustLen: minimal cluster length allowed for initial clustering
% * nMaxGapUnify: maximal gap between clusters to unify them into one 
% * pMax4signif: maximal p value to keep in 'clustMaskSgnf' & 'pValsSignif'
% * bUnifyPosNeg: allow unification of positive and negative clusters?
% * nIter: number of bootstrapping iterations 
%  Default values for these parameters set in GetParams() below
% 
% The function returns a mask the size of the second and onward dimensions
% of cond1 (or cond2; the first dimension is the one compared in divergFunc 
% and hence lost). The values in 'clustMask' are k for cluster #k
% (negative for clusters with t values < 0; 0 for no cluster). 'pVals' are
% the pValues of all the clusters according to the bootstrapping
% distribution. 'clustMaskSgnf' and 'pValsSignif' are only the significant
% clusters and p values, as defined 
% 
% This is an offshoot of Maris & Oostenveld (J Neurosci Meth 2007). The
% code is meant to overcome the multiple-comparison problem -- e.g., how to
% go from local t-test significance to overall significance, when the
% bonferroni correction is too conservative. 
% 
% Uri Maoz, Caltech, Created 9/17/2014. Modified and improved since.
% Uri Maoz, UCLA, Caltech: Significant modification, 11/16/2015: made the
% function general, accepts any divergence function.

[divergFunc,divergTh,neighborMat,neighborMatTh,nMinClustLen,nMaxGapUnify,...
    pMax4signif,bUnifPosNeg,nIter]=GetParams(params,cond1);

n1=size(cond1,1); n2=size(cond2,1);
conds=[cond1;cond2]; lbls=[ones(n1,1);-ones(n2,1)];
bsMaxT=zeros(1,nIter);
% Run actual data and save divergence-statistic sums over all clusters.
% Then run bootstrapped data and calculate largest divergence-statistic sum
% over all clusters for each bootsrap iteration.
rmi(-1);
[~,tSums,clustMask,pVals]=OneIter(conds,lbls,divergFunc,divergTh,...
    neighborMat,neighborMatTh,nMinClustLen,nMaxGapUnify,bUnifPosNeg,false);
if ((~isempty(pVals))&&(~any(isfinite(pVals)))), return; end
% iterLbls=NaN(length(lbls),nIter);
% parfor k=1:nIter
for k=1:nIter
%     [bsMaxT(k),~,~,~,iterLbls(:,k)]=...
    bsMaxT(k)=OneIter(conds,lbls,divergFunc,divergTh,neighborMat,...
        neighborMatTh,nMinClustLen,nMaxGapUnify,bUnifPosNeg,true);
    rmi(sprintf('; Done iteration %i of %i',k,nIter));
end

u=unique(clustMask(:)); uNeg=flipud(u(u<0)); uPos=u(u>0); 
for k=1:length(uNeg), clustMask(clustMask==uNeg(k))=-k; end
for k=1:length(uPos), clustMask(clustMask==uPos(k))=k; end
clustMask(~ismember(clustMask,-length(uNeg):length(uPos)))=0;
clustMask=permute(clustMask,[3 2 1]);

pVals=NaN(1,length(tSums));
for c=1:length(tSums)
    pVals(c)=1-sum(abs(tSums(c))>bsMaxT)/nIter;
end

p=[pVals(1:length(uNeg))';0;pVals(find(u>0)-1)']';
clustMaskSgnf=clustMask; 
clustMaskSgnf(ismember(clustMaskSgnf,u(p>pMax4signif)))=0; 
u=unique(clustMaskSgnf(:)); uNeg=flipud(u(u<0)); uPos=u(u>0); 
for k=1:length(uNeg), clustMaskSgnf(clustMaskSgnf==uNeg(k))=-k; end
for k=1:length(uPos), clustMaskSgnf(clustMask==uPos(k))=k; end
clustMaskSgnf(~ismember(clustMaskSgnf,-length(uNeg):length(uPos)))=0;
pValsSignif=pVals(pVals<pMax4signif);
end

% ============================== Subfunctions =============================
function [bsMaxT, tSums, clustMask, pVals, iterLbls] = ...
    OneIter (conds, lbls, divergFunc, divergTh, neighbor, neighborTh, ...
    nMinClustLen, nMaxGapUnify, bUnifyPosNeg, bBootstrap)
% Run one iteration of either the real data (BS=false) or bootstrapping
% (BS=true).
if (bBootstrap)
    % In bootsrapping mode, shuffle the conditions
    iterLbls=lbls(randperm(length(lbls)));
else
    % For actual data, compute the conditions as is
    iterLbls=lbls;
end
currCond1=conds(iterLbls==1,:,:,:);
currCond2=conds(iterLbls==-1,:,:,:);
divergVals=divergFunc(currCond1,currCond2);
divergVals=permute(divergVals,[4 3 2 1]);
% Find all clusters
% Search for positive and negative clusters separately and then
% unify them in the masking matrix
divergValsRows=divergVals(:);
if (bUnifyPosNeg) % Unify positive and negative clusters
    if (any(abs(divergValsRows)>=divergTh))
        bsClustMask=...
            FindClusts3D(abs(divergVals)>=divergTh,neighbor<=neighborTh,...
            nMinClustLen,nMaxGapUnify);
    else
        bsClustMask=zeros(size(divergVals));
    end
else
    if (any(divergValsRows>=divergTh))
        bsPosClustMask=...
            FindClusts3D(divergVals>=divergTh,neighbor<=neighborTh,...
            nMinClustLen,nMaxGapUnify);
    else
        bsPosClustMask=zeros(size(divergVals));
    end
    if (any(divergValsRows<=-divergTh))
        bsNegClustMask=...
            FindClusts3D(divergVals<=-divergTh,neighbor<=neighborTh,...
            nMinClustLen,nMaxGapUnify);
    else
        bsNegClustMask=zeros(size(divergVals));
    end
    bsClustMask=bsPosClustMask-bsNegClustMask;
end

u=unique(bsClustMask(:)); u=u(u~=0);
tBSSums=zeros(1,length(u));
for c=1:length(u)
    tBSSums(c)=sum(divergVals(bsClustMask==u(c)));
end
if (bBootstrap)
    % Find the largest patch according to sum of t values
    if (isempty(tBSSums))
        bsMaxT=0;
    else
        bsMaxT=max(abs(tBSSums));
    end
    tSums=[]; clustMask=[]; pVals=[];
else
    tSums=tBSSums;
    clustMask=bsClustMask;
    bsMaxT=[];
    if (isempty(tBSSums))
        pVals=NaN(1,length(tSums));
    else
        pVals=[];
    end
end
end

% -------------------------------------------------------------------------
function [divergFunc, divergTh, neighborMat, neighborMatTh, ...
    nMinClustLen, nMaxGapUnify, pMax4signif, bUnifyPosNeg, nIter] = ...
    GetParams(params, cond1)
% The function transforms the parameters structure into individual
% variables and assigns them default values if they do not exist
if (isfield(params,'divergFunc'))
    divergFunc=params.divergFunc; 
else
    divergFunc=[]; 
end
if (isfield(params,'divergTh'))
    divergTh=params.divergTh; 
else
    divergTh=2; 
end
if (isfield(params,'neighborMat')) && ~isempty(params.neighborMat)
    neighborMat=params.neighborMat; 
else
    neighborMat=ones(size(cond1,3),size(cond1,3)); 
end
if (isfield(params,'neighborMatTh'))
    neighborMatTh=params.neighborMatTh; 
else
    neighborMatTh=1; 
end
if (isfield(params,'nMinClustLen'))
    nMinClustLen=params.nMinClustLen; 
else
    nMinClustLen=10; 
end
if (isfield(params,'nMaxGapUnify'))
    nMaxGapUnify=params.nMaxGapUnify; 
else
    nMaxGapUnify=round(nMinClustLen/3); 
end
if (isfield(params,'pMax4signif'))
    pMax4signif=params.pMax4signif; 
else
    pMax4signif=.5; 
end
if (isfield(params,'bUnifyPosNeg'))
    bUnifyPosNeg=params.bUnifyPosNeg; 
else
    bUnifyPosNeg=false;
end
if (isfield(params,'nIter')), nIter=params.nIter; else nIter=1000; end
end