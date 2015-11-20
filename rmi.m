function hTic = rmi (dur, msg, hTic)
% function hTic = rmi (dur, msg, nTic)
% 
% If more than 'dur' seconds (default 10) have passed since the previous
% call to rmi, the function outputs where the program was when it called
% 'rmi', and how many seconds elapsed. It appends any message in 'msg'. The
% function must be initialized by the call rmi(-1).
% rmi (msg) uses the default sec=10.
% rmi is short for "where am i?"
persistent first
persistent last
persistent tstart

dflt_dur=10;
if (~exist('dur','var')||isempty(dur)), dur=dflt_dur; end
if (ischar(dur))
    if (nargin>1), hTic=msg; end
    msg=dur; dur=dflt_dur; 
end
if (~exist('msg','var')), msg=''; end
if (~exist('hTic','var')), hTic=[]; end
if (dur==-1);
    tstart=tic;
    hTic=tstart;
    first=toc(tstart); last=toc(tstart);
else
    if (toc(tstart)-last>dur)        
        if (isempty(hTic))
            last=toc(tstart);
        else
            last=toc(hTic);
        end
        st=dbstack;
        fprintf('%0.1fs elapsed. In "%s" ln %i%s\n',...
            last-first,st(end).file,st(end).line,msg);
    end
end
end