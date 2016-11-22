function [str]=writecommaformat(a,f)
% translate a 1 or 2 dimensional array into a string

if(~isvector(a))
    error('this function only works for vectors')
end
    
str = '';

if(~exist('f','var'))
    f = '%f';
end

% unpack cell singleton
if(iscell(a) && length(a)==1)
    a = a{1};
end

if(size(a,2)>1)
    a = a';
end 

clear strj
for j=1:size(a,1)
    % colon separate a(i,j,:)
    strj{j} = sprintf(f,a(j));
%     for k =2:size(a,2)
%         strj{j} = [strj{j} ':' sprintf(f,a(j,k))];
%     end
end
str = strj{1};
for j=2:size(a,1)
    str = [str ',' strj{j}];
end

