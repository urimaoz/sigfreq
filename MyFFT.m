function [freq, mag, ang, F, a, b, len] = ...
    MyFFT (f, t, limitFFT, figNum, titleStr)
% function [freq, mag, ang, F, a, b, len] = ...
%     MyFFT (f, t, limitFFT, fig_num, title_str)
% 
% This function computes the FFT for signal 'f' over time 't' and
% optionally draws it. If 'limitFFT' (default []) is empty, half of
% length(t) is ploted. If it is a frequency>1, the data is only computed up 
% to that frequency. If it is between 0 and 1, the frequencies axis is cut
% off where all the higher frequencies 'limitFFT' ratio are below limitFFT
% of the maximum magnitude. If 'fig_num' (default 10) is given it is used
% as the figure on which to draw. If it is 0, the current axes are used; if
% [] nothing is drawn. 'title_str' (default '') is added to the titles of
% the axes.  
% It returns the returns the frequencies in 'freq', the magnitude (sqrt of
% power) 'mag' and angle 'ang' of the Fourier transform. It also returns
% the Fourier transform 'F' and the  corresponding frequencies (the x-axis
% of the generated plots, if generated). It further outputs 'a' and 'b'
% hold the coefficients of the cosine and sine parts of the Fourier
% transform respectively. 'len' holds the effective length of the FFT (as
% determined by 'limitFFT').  

removeZerosTail = 0;
if (nargin < 5)
    titleStr = '';
	if (nargin < 4)
        figNum = 10;
		if (nargin < 3)
            limitFFT = [];%round (length (t) ./ 2);
		end
	end
end

% dt=mean(diff(t)); 
% Fs=1./dt;
T=t(end)-t(1);
N0=length(t);
row2col=(size(f,1)<size(f,2));
t=t(:);
f=f(:);
% nFFT=2^nextpow2(N0);
% F=fft(f,nFFT)/N0;
F=fft(f)/N0;
if (row2col), F=F'; t=t'; end

[Fp,Fm]=cart2pol(real(F),imag(F)); Fm=Fm*2;
% Fundamental frequency(FF)=sampling frequency(SF)/#data point-1(N).
% SF=1/dt => FF=1/(dt*N). dt*N=T => FF=1/T.
if (row2col)
    freq=(0:(N0-1))/T;
else
    freq=(0:(N0-1))'/T;
end

if isempty (limitFFT)
    limitFFT=freq(round(length(t)./2));
elseif ((limitFFT>0) && (limitFFT<1))
    removeZerosTail=limitFFT;
elseif (limitFFT==0)
    limitFFT=freq(length(t));
end

if (removeZerosTail>0) % Remove 
    absFm=abs(Fm(1:round(length(t)./2)));
    relatCumulPower=cumsum(absFm)/sum(absFm);
    rmInds=find(relatCumulPower>limitFFT);
    limitFFT=rmInds(1);
    if (limitFFT<2), error('"limitFFT" is too small'); end
end
    
i=(freq<=limitFFT);
if (~isempty(figNum))
    if (figNum>=1), figure (figNum); clf; end
    plot(freq(i),Fm(i),'.-'); 
	xlabel ('Frequency'); ylabel ('FFT magnitude'); axis tight; 
    if (~isempty(titleStr)), title (titleStr); end
end

a=2/T.*real(F); 
b=-2/T.*imag(F); 

freq=freq(i);
mag=Fm(i);
ang=Fp(i);
len=limitFFT;
end