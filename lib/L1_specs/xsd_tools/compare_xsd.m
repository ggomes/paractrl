function [result,why] = compare_xsd(fileA,fileB)
% compare two xsd files.
% return result=true if they are logically identical, result=false otherwise
% why: list of differences
% author: Gabriel Gomes

[result, why] = eq_xsd(xml_read(fileA),xml_read(fileB));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, why] = eq_xsd(struct1, struct2)

why = {};

% have we reached a leaf?
if(~isstruct(struct1) || ~isstruct(struct2))
    result = isequal(struct1,struct2);
    if(~result)
        why = {'Values are different'};
    end
    return
end

% order fields
struct1 = orderfields(struct1);
struct2 = orderfields(struct2);

% extract field names
fields1 = fieldnames(struct1);
fields2 = fieldnames(struct2);

% same size arrays?
if any(size(struct1) ~= size(struct2))
    result = false;
    why = {'Sizes are different'};
    return
end

% Same number of fields?
if length(fields1) ~= length(fields2)
    result = false;
    why = {'Number of fields are different'};
    return
end

% Same field names?
result = celleq(fields1,fields2);
result = all(result);
if ~result
    why = {'Field names are different'};
    return
end

%  order struct array by name (for xsd comparisons)
if(isfieldRecursive(struct1,'ATTRIBUTE','name'))
    for i=1:length(struct1)
        name1{i} = struct1(i).ATTRIBUTE.name;
        name2{i} = struct2(i).ATTRIBUTE.name;
    end
    in1not2 = setdiff(name1,name2);
    for i=1:length(in1not2)
        why = [why {['First contains "' in1not2{i} '", second does not.']}];
    end
    in2not1 = setdiff(name2,name1);
    for i=1:length(in2not1)
        why = [why {['Second contains "' in2not1{i} '", first does not.']}];
    end
    if(~isempty(why))
        result = false;
        return;
    end
    [ismem,ind]=ismember(name1,name2);
    struct2(ismem) = struct2(ind(ismem));
    clear name1 name2    
end

% traverse array, check each branch
for i = 1:numel(struct1)
    for j=1:length(fields1)        
        [subresult, subwhy] = structeq1(struct1(i).(fields1{j}),struct2(i).(fields2{j}));
        if(~subresult)
            for k=1:length(subwhy)
                why = [why;{sprintf('(%d).%s: %s',i,fields1{j},subwhy{k})}];
                result = false;
            end
        end
    end
end