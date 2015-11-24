function res = dPrime (x, y, dim, robust)
% function res = dPrime (x, y, dim, robust)
% 
% Compute d' of 'x' and 'y' along dimention 'dim':
% res=(mean(x)-mean(y))/sqrt(0.5*(var(x)+var(y));
% If 'robust' is true (default is false), compute a robust version of d':
% res=(median(x)-median(y))/sqrt(0.5*((iqr(x)./1.349)^2+(iqr(y)./1.349)^2);
% 
% Written by Uri Maoz, urim@caltech.edu 18 Jan 2015
if (~exist('y','var')||isempty(y)), y=zeros(size(x)); end
if (~exist('dim','var')||isempty(dim)), dim=1; end
if (~exist('robust','var')||isempty(robust)), robust=false; end
if (robust)
    mX=nanmedian(x,dim); mY=nanmedian(y,dim); 
    sX=iqr(x,dim)/1.349; sY=iqr(y,dim)/1.349; 
else
    mX=nanmean(x,dim); mY=nanmean(y,dim); 
    sX=nanstd(x,0,dim); sY=nanstd(y,0,dim); 
end
res=(mX-mY)./sqrt((sX.^2+sY.^2)./2);
end
