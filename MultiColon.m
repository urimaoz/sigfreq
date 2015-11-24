function i = MultiColon (i1, i2)
% function i = MultiColon (i1, i2)
% 
% The function takes two arrays, 'i1' and 'i2', both (1 x n) and returns
% [i1(1):i2(1) i1(2):i2(2) ... i1(n):i2(n)]
% if (length(i1)~=length(i2))
%     error('The function''s two input arrays must be of the same length');
% end
i=cell2mat(cellfun(@(x,y) x:y,num2cell(i1),num2cell(i2),...
    'UniformOutput',false));
end
