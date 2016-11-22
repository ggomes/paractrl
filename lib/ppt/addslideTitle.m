function addslideTitle(op,titletext)
% addslideText(op,titletext)

% Get current number of slides:
slide_count = get(op.get('Slides'),'Count');

% Add a new slide (with title object):
slide_count = int32(double(slide_count)+1);
slide = invoke(op.get('Slides'),'Add',slide_count,1);

% Insert text into the title object:
if(~isempty(titletext))
    str{1} = titletext{1};
    c = 2;
    for i=2:length(titletext)
        str{c} = char(13);
        str{c+1} = titletext{i};
        c = c+2;
    end
    set(slide.Shapes.Title.TextFrame2.TextRange,'Text',char(str{:})');
end

return