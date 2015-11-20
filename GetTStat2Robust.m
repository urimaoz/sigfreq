function [tstat] = GetTStat2Robust (x, y, varargin)
% function [tstat] = GetTStat2Robust (x, y, varargin)
% 
% This is a special version of 'ttest2' that is quicker in getting back
% only the t statistic (usually a part of the 'stats' output variable). It
% is also built to be more robust to outliers, relying on median and IQR
% instead of mean and STD.
%GETTSTAT2ROBUST Two-sample t-test with pooled or unpooled variance estimate.
%   T = GETTSTAT2ROBUST(X,Y) performs a t-test of the hypothesis that two
%   independent samples, in the vectors X and Y, come from distributions
%   with equal means, and returns the resulting test statistic in T. The
%   data are assumed to come from normal distributions with unknown, but
%   equal, variances. X and Y can have different lengths.
%
%   This function performs an unpaired two-sample t-test. For a paired
%   test, use the GETTSTATROBUST function.
%
%   X and Y can also be matrices or N-D arrays.  For matrices,
%   GETTSTAT2ROBUST performs separate t-tests along each column, and
%   returns a vector of results.  X and Y must have the same number of
%   columns.  For N-D arrays, GETTSTAT2 works along the first non-singleton
%   dimension.  X and Y must have the same size along all the remaining
%   dimensions. 
%
%   GETTSTAT2ROBUST treats NaNs as missing values, and ignores them.
%
%
%   [...] = GETTSTAT2ROBUST(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
%   one or more of the following name/value pairs:
%
%       Parameter       Value
%       'vartype'       'equal' to perform the default test assuming equal
%                       variances, or 'unequal', to perform the test
%                       assuming that the two samples come from normal
%                       distributions with unknown and unequal variances.
%                       This is known as the Behrens-Fisher problem. TTEST2
%                       uses Satterthwaite's approximation for the
%                       effective degrees of freedom.
%       'dim'           Dimension DIM to work along. For example, specifying
%                       'dim' as 1 tests the column means. Default is the
%                       first non-singleton dimension.
%
%   See also GETTSTATROBUST, GETTSTAT2, TTEST2, TTEST, RANKSUM, VARTEST2, ANSARIBRADLEY.

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 13.4. (Table 13.4.1 on page 210)

%   Copyright 1993-2012 The MathWorks, Inc.


if nargin < 2
    error(message('stats:ttest2:TooFewInputs'));
end

% Process remaining arguments
vartype = '';
dim = '';

if nargin>=3
    if isnumeric(varargin{1})
        % Old syntax
        %    TTEST2(X,Y,ALPHA,TAIL,VARTYPE,DIM)
        vartype = varargin{1};
        if nargin>=4
            dim = varargin{2};
        end
        
    elseif nargin==3
            error(message('stats:ttest2:BadAlpha'));
    
    else
        % Calling sequence with named arguments
        okargs =   {'vartype' 'dim'};
        defaults = {''      ''};
        [vartype,dim]=internal.stats.parseArgs(okargs,defaults,varargin{:});
    end
end

if isempty(vartype)
    vartype = 1;
elseif isnumeric(vartype) && isscalar(vartype) && ismember(vartype,[1 2])
    % OK, grandfathered
else
    [~,vartype] = internal.stats.getParamVal(vartype,...
        {'equal','unequal'},'''vartype''');
end

if isempty(dim)
    % Figure out which dimension mean will work along by looking at x. y
    % will have be compatible. If x is a scalar, look at y.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = find(size(y) ~= 1, 1); end
    if isempty(dim), dim = 1; end
    
    % If we haven't been given an explicit dimension, and we have two
    % vectors, then make y the same orientation as x.
    if isvector(x) && isvector(y)
        if dim == 2
            y = y(:)';
        else % dim == 1
            y = y(:);
        end
    end
end

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
    error(message('stats:ttest2:InputSizeMismatch'));
end

xnans = isnan(x);
if any(xnans(:))
    nx = sum(~xnans,dim);
else
    nx = size(x,dim); % a scalar, => a scalar call to tinv
end
ynans = isnan(y);
if any(ynans(:))
    ny = sum(~ynans,dim);
else
    ny = size(y,dim); % a scalar, => a scalar call to tinv
end

s2x = (iqr(x,dim)./1.349).^2;
s2y = (iqr(y,dim)./1.349).^2;
difference = nanmedian(x,dim) - nanmedian(y,dim);
if vartype == 1 % equal variances
    dfe = nx + ny - 2;
    sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
    se = sPooled .* sqrt(1./nx + 1./ny);
    tstat = difference ./ se;

elseif vartype == 2 % unequal variances
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
%     dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
    se = sqrt(s2xbar + s2ybar);
    tstat  = difference ./ se;
end

