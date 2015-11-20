% %% Test MyFFT(). It is old code. For a series of sinusoids with powers
% % decreasing as 1/f^a, what do we see with cftool?
nLoops=50;
% t=0:.001:1; 
% a=1.25;
% res=zeros(1,nLoops);
% rmi(-1);
% for i=1:nLoops
%     x=0;
%     for k=1:.01:100
%         x=x+sin(k*pi*t+rand(1)*2*pi-pi)/k^a;
%     end;
%     [f,m]=MyFFT(x,t,[],0); m=m(f>0)'; f=f(f>0)';
%     model=fit(f,m,'power1');
%     res(i)=model.b;
%     rmi(sprintf('; Done loop %d of %d, ratio %g',i,nLoops,i/nLoops));
% end
% 

% cftool(f(f>0),m(f>0)); % Note: Choose POWER in top right dropdown menu
% cftool(log(f(f>0)),log(m(f>0)))

%% The above seems to work. So I assume MyFFT() is OK. Now create colored 
% noise using Matlab's DSP toolbox and test it
a=.5:.1:1.5;
res2=zeros(length(a),nLoops);
n=2048; 
t=linspace(0,1,n);
rmi(-1);
for j=1:length(a)
    for i=1:nLoops
        hColoredNoise=dsp.ColoredNoise(a(j),n,1);
        % rng default; % Create repeatable randomness for debugging purposes only
        x=step(hColoredNoise);
        [f,m]=MyFFT(x,t,[],[]);  
        fInd=((f>1)&(f<500)); f=f(fInd); m=m(fInd); 
        model=fit(f,m,'power1');
        res2(j,i)=model.b;
        rmi(sprintf('; a=%g -- Done loop %d of %d, ratio %g',...
            a(j),i,nLoops,i/nLoops));
    end
end

% cftool(f(f>0),m(f>0)); % Note: Choose POWER in top right dropdown menu
% cftool(log(f(f>0)),log(m(f>0)))
%%
% squeeze(mean(res2,2)), squeeze(std(res2,0,2))
figure(10); clf; hold on;
for j = 1:length(a)
    plot(a(j),-2*res2(j,:),'k.');
    plot(a(j),-2*mean(res2(j,:)),'sg');
    [~,p]=ttest(-2*res2(j,:),a(j));
    text(a(j),2,sprintf('p=%0.3f',p),'FontSize',14);
end
plot(a,a,'-r');
