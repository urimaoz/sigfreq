function tstat = GetTStatRobust(x,m,varargin)
%GETSTAT  Robust version of one-sample and paired-sample t-test.
% 
%   The function acts the same way as GetTStat. But it is a robust version.
%   Thus it attempts to overcome noisy samples using the median as an
%   approximation for the mean and IQR as an approximation for the STD.
%   See also GETTSTAT2ROBUST, GETTSTAT, GETTSTAT2, 


if nargin < 2 || isempty(m)
    m = 0;
elseif ~isscalar(m) % paired t-test
    if ~isequal(size(m),size(x))
        error(message('stats:ttest:InputSizeMismatch'));
    end
    x = x - m;
    m = 0;
end

% Process remaining arguments
dim = '';

if nargin>=3
    if isnumeric(varargin{1})
        % Old syntax
        %    TTEST(X,M,ALPHA,TAIL,DIM)
        dim = varargin{1};
                
    elseif nargin==3 
            error(message('stats:ttest:BadAlpha'));
   
    else
        % Calling sequence with named arguments
        okargs =   {'dim'};
        defaults = {''};
        dim = internal.stats.parseArgs(okargs,defaults,varargin{:});
    end
end

if isempty(dim)
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

nans = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end
xmean = nanmedian(x,dim);
% sdpop = MyIQR(x,dim)/1.349; % iqr() works with NaN's
sdpop = iqr(x,dim)/1.349;
ser = sdpop ./ sqrt(samplesize);
tstat = (xmean - m) ./ ser;

