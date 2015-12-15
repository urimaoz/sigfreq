function beginAndEnd = continuousRunsOfTrue(vectorOfLogicals)
% SYNTAX: beginAndEnd = continuousRunsOfTrue(vectorOfLogicals)
%
% A handy function that finds the start and end indices of all of the
% continous runs of true within a vector. beginAndEnd is an n-by-2 array,
% where n is the number of continuous runs of true in vectorOfLogicals.
% Each row represents one run. The first column is the starting index, the
% second column is the ending index. Note that singleton trues are
% included. So if vectorOfLogicals is [0 0 1 1 0 1 0], beginAndEnd will
% return [3 4; 6 6];
%
% Written by Emily Mankin, (c) 2012


if size(vectorOfLogicals,1)>1
    vectorOfLogicals = vectorOfLogicals';
end
if numel(vectorOfLogicals)~=length(vectorOfLogicals)
    warning(['continuousRunsOfTrue only supports single-dimension vectors. ',...
        'Your data is being converted'])
    vectorOfLogicals = vectorOfLogicals(:);
end
if sum(vectorOfLogicals)==0
    beginAndEnd = zeros(0,2);
elseif sum(vectorOfLogicals)==length(vectorOfLogicals)
    beginAndEnd = [1,length(vectorOfLogicals)];
else
    v = vectorOfLogicals;
    dv = diff(v);
    endPoints = find(dv==-1);
    startPoints = find(dv==1)+1;
    if vectorOfLogicals(1)
        startPoints = [1,startPoints];
    end
    if vectorOfLogicals(end)
        endPoints = [endPoints,length(v)];
    end
    beginAndEnd = [startPoints',endPoints'];
end