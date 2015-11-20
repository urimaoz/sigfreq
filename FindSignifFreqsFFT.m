function [freqRanges, pVals, tValsSum] = ...
    FindSignifFreqsFFT (trials, time, freqLim, smoothSpan, tTH, figures)
% function [freqRanges, pVals, tValsSum] = ...
%     FindSignifFreqsFFT (trials, time, freqLim, smoothSpan, tTH, figures)
% 
% Given 'trials' data (#trials x #samples), the 'time' (1 x #samples), the
% frequency limits, 'freqLim' [minimum frequency, maximum frequency], a
% smoothing span on the data, 'smoothSpan', and a threshold on the t test,
% 'tTH', the function computes which bands in the frequency are higher or
% lower than would be expected by chance. The function can also plot this
% out, if 'figures' is true. If length(figures)>1, the function also plots
% out the log-log plots of power vs. frequency with the robust regression
% and the computed alpha value. 
% The function is a spinoff of the Maris & Oostenvald (2007) method.
if (~exist('figures','var')||isempty(figures)), figures=false; end
if (exist('freqLim','var')&&~isempty(freqLim)&&(length(freqLim)==1))
    freqLim=[eps,freqLim]; 
end

% Minimum # Hz t statistic chunks, with minimum unification chunks
minHzWidth=1; minHzUnif=minHzWidth/3;

% Find t statistics of actual data
if (length(figures)>1), figure(figures(2)); clf; end
    % The smoothing in the function below needs to be replaced by a
    % multi-taper approach @@@ 
[tStat,logF,logFreq,logBase,allAlpha]=...
    Data2tStat(trials,time,freqLim,smoothSpan,length(figures)>1);

% Run code on random 1/f^a data
nIters=100;
% rmi(-1);
allExtr=NaN(2,nIters);
% dPrimeRnd=NaN(nIters,size(logF,2));
allAlphaRnd=NaN(nIters,length(allAlpha));
allTstatRand=NaN(nIters,length(tStat));
% % %%% Remove
% % figure(30); clf; 
% % plot3(logFreq,tStat,ones(size(tStat)),...
% %     logFreq([1 end]),[-2 -2],[1 1],logFreq([1 end]),[2 2],[1 1]); 
% % view (2); hold on
% % 
% % figure(31); clf; 
% % plot3(logFreq,median(logF),ones(size(logF)),'b'); 
% % view(2); hold on;
% % %%% Remove - end

tic;
parfor k=1:nIters
    % Create 1/f^alpha noise
    trialsRnd=NaN(size(trials));
    for n=1:size(trialsRnd,1)
%         trialsRnd(n,:)=PowerNoiseFilt(-allAlpha(n),size(trials,2)); %%#ok<PFBNS>
        trialsRnd(n,:)=PowerNoise(-allAlpha(n),size(trials,2),...
            'normalize','randpower'); %#ok<PFBNS> 
    end
    % Compute t statistics and chunk them, and compute d'
    [allTstatRand(k,:),~,~,~,allAlphaRnd(k,:)]=...
        Data2tStat(trialsRnd,time,freqLim,smoothSpan,false);
    valsRnd=TStat2Chunks(allTstatRand(k,:),logFreq,tTH,...
        minHzWidth,minHzUnif,logBase);
    if (~isempty(valsRnd))
        allExtr(:,k)=[min(valsRnd);max(valsRnd)];
    end
%     rmi(sprintf('; Done bootstrapping iternation %i of %i',k,nIters));
% %     figure(30); plot(logFreq,tStatRnd,'c')
% %     figure(31); plot(logFreq,median(logFRnd),'c');
end
fprintf('Overall bootstrapping ran for %0.2f s\n',toc);

% Only t statistics that are outside the 95% CI of the noise should be
% considered as part of significant chunks
rndTpctls=prctile(allTstatRand,[2.5 50 97.5]);
iTsignif=((tStat<rndTpctls(1,:))|(tStat>rndTpctls(3,:)));
t4chunks=tStat; t4chunks(~iTsignif)=0;
[tValsSum,inds]=TStat2Chunks(t4chunks,logFreq,tTH,minHzWidth,minHzUnif,...
    logBase);

pVals=getPvals(tValsSum,allExtr);
freqRanges=logBase.^logFreq(inds); 

%  Plotting
if (figures(1)>0)
    if (length(figures)>1), figure(figures(1)); clf; end
    cols='bcr';
    pTH=.05;
    % Plot power
    subplot(211); 
    m=median(logF,1); 
    s=iqr(logF,1)./1.349; 
    h=plot([3 4],[0 0],['-',cols(1)],[3 4],[0 0],['--',cols(1)],...
        [3 4],[0 0],['.-',cols(2)],[3 4],[0 0],['.-',cols(3)]);
    hold on;
    plot(logFreq,m,'-b',logFreq,m-s,'--b',logFreq,m+s,'--b',...
        logFreq([1 end]),[0 0],'-k'); 
    ylabel('log2(power), 1/f corrected');
    iSig=pVals<=pTH;
    for k=1:size(inds,2)
        currLogFreq=logFreq(inds(1,k):inds(2,k));
        plot(currLogFreq,m(inds(1,k):inds(2,k)),['.-',cols(2+iSig(k))]);
    end
    legend('Mean','mean+/-STD','t clusts','Signif t clusts',...
        'Location','Best');
    set(h,'Visible','Off');
    [rho,p]=corr([allAlpha',allAlphaRnd(end,:)']);
    title(sprintf('Corr data & pink noise alphas: %0.2f (p=%0.3f)',...
        rho(1,2),p(1,2)));
    % Plot t staitistics
    subplot(212);
    cols='bkg';
    h=plot([3 4],[0 0],['-',cols(1)],...
        [3 4],[0 0],['--',cols(2)],...
        [3 4],[0 0],['-',cols(3)],[3 4],[0 0],['--',cols(3)]);
    hold on;
    plot(logFreq,tStat,cols(1)); 
%     rndTpctls=prctile(allTstatRand,[2.5 50 97.5]);
    plot(logFreq([1 end]),[0 0],['-',cols(2)],... 
        logFreq([1 end]),[-tTH -tTH],['--',cols(2)],...
        logFreq([1 end]),[tTH tTH],['--',cols(2)],... 
        logFreq,rndTpctls(2,:),['-',cols(3)],...        
        logFreq,rndTpctls(1,:),['--',cols(3)],...
        logFreq,rndTpctls(3,:),['--',cols(3)]);
    ylabel('t stat'); xlabel('log2(frequency)');
    legend('Actual data',sprintf('|t threshold|=%g',tTH),...
        'pink noise','noise 95% CI');
    set(h,'Visible','Off');
    ylim([-tTH*5 tTH*5]);
end
end

% ============================ Subfunctions ===============================
function [tStat, logF, logFreq, logBase, allAlpha] = ...
    Data2tStat (trials, time, freqLim, smoothSpan, draw)
% Compute the t statistics over trials data
if (~exist('draw','var')||isempty(draw)), draw=false; end
nTrials=size(trials,1);
if (draw)
    % The screens are typically 9x16 now, so make the subplots accordingly
%     I=floor(sqrt(nTrials)); J=ceil(sqrt(nTrials)); 
    I=ceil(sqrt(nTrials/(9*16))*9); J=floor(sqrt(nTrials/(9*16))*16); 
    while (I*J<nTrials), J=J+1; end
    subplot(I,J,1);
end
[logF1,logFreq,alpha,logBase]=...
    FreqPeaks(trials(1,:),time,freqLim,smoothSpan, draw);
if (draw), title(sprintf('Tr %d, a=%0.2f',1,alpha)); axis('tight'); end
logF=NaN(size(trials,1),length(logF1)); logF(1,:)=logF1;
allAlpha=NaN(1,size(trials,1)); allAlpha(1)=alpha;
for k=2:nTrials
    if (draw), subplot(I,J,k); end
    [logF(k,:),~,allAlpha(k)]=...
        FreqPeaks(trials(k,:),time,freqLim,smoothSpan,draw);
    if (draw)
        title(sprintf('Tr %d, a=%0.2f',k,allAlpha(k))); 
        axis('tight'); 
        if (k==nTrials)
            xlabel('log freq'); xlabel('log power'); 
        end
    end
end
% figure(20); hold on; plot(2.^logFreq,2.^mean(logF))
tStat=GetTStatRobust(logF);
end

% -------------------------------------------------------------------------
function [logF, logFreq, alpha, logBase] = ...
    FreqPeaks (f, t, fLim, span, draw)
% function [logF, logFreq, alpha, logBase] = ...
%     FreqPeaks (f, t, fLim, span, draw)
% 
% Neural data tends to decrease with frequency as 1/(f^alpha). The
% function therefore: 
% 1. Takes the signal 'f' over time 't', finds its frequency power (with a
% Fourier transform). 
% 2. Smoothes the frequency power with the 'loess' method and a span of n
% samples, if 'span'>1 or a% of the data if 'span'<1. Default is no
% smoothing. 
% 3. Converts the data into log-space
% 4. Runs a robust line fit through its log-space representaiton, 
% 5. Subtract this line from the power to get output power 'F' over
% frequencies 'freq'. 
% 
% The function returns:
% logF      the detrended log-space frequency spectrum
% logFreqs  the log-space frequencies
% alpha     the 1/f^(alpha) parameter
% logBase   the base of the logarithm

[freq,F]=MyFFT(f,t,[],[]);
if (exist('fLim','var')&&~isempty(fLim))
    % If 'fLim' is a scalar, assume it is the high cutoff of frequency
    if (length(fLim)==1), fLim=[eps,fLim]; end
    iIncl=((freq>=fLim(1))&(freq<=fLim(2)));
else
    iIncl=true(size(freq));
end
if (~exist('span','var')||isempty(span))
    span=1; % No smoothing
end
[logF,logFreq,b,logBase]=FindPeaks(F.^2,freq,iIncl,span,draw);
alpha=b(2);
end

% -------------------------------------------------------------------------
function [Fsm, logFreq, b, logBase] = ...
    FindPeaks (F, freq, iIncl, span, draw)
%  An internal function that finds the frequency peaks for the actual data
%  and for the bootstrapping
freq=freq(iIncl);
if (~exist('span','var')||isempty(span))
    span=1; % No smoothing
end
if (span~=1)
    % Smooth 10% past the end of the iIncl and then run the iIncl to avoid
    % edge artifacts
            % This approach should be replaced by multi-tapering @@@
    i1=find(iIncl,1,'first'); iE=find(iIncl,1,'last'); 
    iSmoothEnd=round(i1+1.1*(iE-i1+1));
    try
        Fsm=smooth(F(1:iSmoothEnd),span,'loess')';
    catch
        Fsm=Smooth(F(1:iSmoothEnd),span)';
    end
    Fsm=Fsm(iIncl(1:iSmoothEnd));
else
    Fsm=F(iIncl);
end
% A hack to make sure we do not get <0 Fsm, resulting in complex logs
Fsm(Fsm<0)=min(Fsm(Fsm>0)); 
logBase=2;
logFreq=log(freq)/log(logBase);
logF=log(Fsm)/log(logBase);
warnState=warning('off','stats:statrobustfit:IterationLimit');
    % This approach should be replaced by non-linear power fitting using
    % the fit() function @@@
b=robustfit(logFreq(:),logF(:)); 
warning(warnState);
Fsm=logF-(b(1)+b(2)*logFreq);
if (draw)
    plot(logFreq,logF,logFreq([1 end]),logFreq([1 end])*b(2)+b(1));
end
end

% -------------------------------------------------------------------------
function [vals, inds] = TStat2Chunks (tStat, logFreq, tTH, ...
    minHzWidth, minHzUnif, logBase)
% Given the continuous t-statistic values, gather those > tTH or <-tTH into
% chunks. Discard chunks < min length and unify chunks < min unify apart
maxChunksNum=1000;
vals=NaN(1,maxChunksNum); inds=NaN(2,maxChunksNum); iVals=0;
i=[(tStat<-tTH);(tStat>tTH)]; 
for k=1:size(i,1)
    % Unify clusters that are close enough
    i(k,:)=UnifyClusts1D(i(k,:),logFreq,minHzUnif,logBase);
    if (any(i(k,:)))
        % Omit clusters that are too small
        [i1,iE]=Bool2Chunks(i(k,:));
        iRm=(logBase.^logFreq(iE)-logBase.^logFreq(i1)<minHzWidth);
        i1(iRm)=[]; iE(iRm)=[];
        % For each t statistic cluster remaining, save the area under the
        % curve, calculating it as a series of trapezoids
        for n=1:length(i1)
            vals(iVals+n)=sum(...
                (tStat(i1(n):(iE(n)-1))+tStat((i1(n)+1):iE(n))).*...
                diff(logFreq(i1(n):iE(n))))/2;
            %         vals(iVals+n)=sum(tStat(i1(n):iE(n)));
        end
        inds(:,iVals+(1:length(i1)))=[i1;iE];
        iVals=iVals+length(i1);
    end
end
vals=vals(1:iVals); inds=inds(:,1:iVals);
end

% -------------------------------------------------------------------------
function c = UnifyClusts1D (m, logFreq, unifyGap, logBase)
% Unify clusters along one dimensional array m and return the result in an
% array c with i for cluster i and 0's elsewhere
i1=find(diff(m)>0)+1; if (m(1)), i1=[1,i1]; end
i2=find(diff(m)<0); if (m(end)), i2=[i2,length(m)]; end
iUnif=(logBase.^logFreq(i1(2:end))-logBase.^logFreq(i2(1:end-1))<=...
    unifyGap);
c=m; 
if (any(iUnif)), c(MultiColon(i2(iUnif),i1([false,iUnif])))=true; end
end

% -------------------------------------------------------------------------
function [i1, iE] = Bool2Chunks (b)
% Given a boolean series (i.e., of true or false), divide it into chunks,
% with i1 holding the beginnig of the chunks and iE their ends

if (~any(b)), error('Should not be called with all-false input'); end
i1=find(diff(b)>0)+1; iE=find(diff(b)<0);
% If reached here, i1 and i2 cannot both be empty
if (isempty(i1)), i1=1; end
if (isempty(iE)), iE=length(b); end
if (i1(1)>iE(1)), i1=[1,i1]; end
if (iE(end)<i1(end)), iE=[iE,length(b)]; end
end

% -------------------------------------------------------------------------
function pVals = getPvals (vals, allExtr)
% Compute where the values 'vals' are in the distribution of maxima and
% minima 'allExtr'
pVals=NaN(size(vals));
for k=1:length(vals)
    if (vals(k)>0)
        pVals(k)=nansum(allExtr(2,:)>vals(k))/sum(isfinite(allExtr(2,:)));
    else
        pVals(k)=nansum(allExtr(1,:)<vals(k))/sum(isfinite(allExtr(1,:)));
    end
end
end