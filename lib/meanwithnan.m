function [A]=meanwithnan(B,dim)

if(isempty(B))
    A=[];
    return
end

if(nargin==1)
    dim=1;
end
    
if(dim==2)
    B=B';
end

for i=1:size(B,2)
    if(all(isnan(B(:,i))))
        A(i) = NaN;
        continue
    end
    A(i) = mean(B(~isnan(B(:,i)),i));
end

if(dim==2)
    A=A';
end