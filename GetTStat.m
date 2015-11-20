function tstat = GetTStat(x,m,varargin)
%GETSTAT  One-sample and paired-sample t-test.
%   TSTAT = GETSTAT(X) performs a t-test of the hypothesis that the data in
%   the vector X come from a distribution with mean zero, and returns the
%   resulting t statistic of the test in TSTAT. The data are assumed to
%   come from a normal distribution with unknown variance.
%
%   X can also be a matrix or an N-D array. For matrices, GETSTAT performs
%   separate t-tests along each column of X, and returns a vector of
%   results.  For N-D arrays, GETSTAT works along the first non-singleton
%   dimension of X.
%
%   GETSTAT treats NaNs as missing values, and ignores them.
%
%   TSTAT = GETSTAT(X,M) performs a t-test of the hypothesis that the data
%   in X come from a distribution with mean M.  M must be a scalar.
%
%   TSTAT = GETSTAT(X,Y) performs a paired t-test of the hypothesis that
%   two matched samples, in the vectors X and Y, come from distributions
%   with equal means. The difference X-Y is assumed to come from a normal
%   distribution with unknown variance.  X and Y must have the same length.
%   X and Y can also be matrices or N-D arrays of the same size.
%
%   [...] = GETSTAT(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following name/value pairs:
%
%       Parameter       Value
%       'dim'           Dimension DIM to work along. For example, specifying
%                       'dim' as 1 tests the column means. Default is the
%                       first non-singleton dimension.
%
%   See also TTEST, TTEST2, ZTEST, SIGNTEST, SIGNRANK, VARTEST.

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, page 206.

%   Copyright 1993-2012 The MathWorks, Inc.


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
% df = max(samplesize - 1,0);
xmean = nanmean(x,dim);
sdpop = nanstd(x,[],dim);
ser = sdpop ./ sqrt(samplesize);
tstat = (xmean - m) ./ ser;

