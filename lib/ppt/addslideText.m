function addslideText(op,titletext,bodytext)
% addslideText(op,titletext,bodytext)

% Get current number of slides:
slide_count = get(op.Slides,'Count');

% Add a new slide (with title object):
slide_count = int32(double(slide_count)+1);
slide = invoke(op.Slides,'Add',slide_count,11);

% Insert text into the title object:
if(~isempty(titletext))
    set(slide.Shapes.Title.TextFrame.TextRange,'Text',titletext);
end

if(~isempty(bodytext))
    fromtop = 60;
    frombottom = 20;
    fromside = 20;
    slide_H = op.PageSetup.SlideHeight;
    slide_W = op.PageSetup.SlideWidth;
    Title1 = slide.Shapes.AddTextbox('msoTextOrientationHorizontal',fromside,fromtop,slide_W-2*fromside,slide_H-fromtop-frombottom);
    Title1.TextFrame.TextRange.Text = bodytext;
end

return