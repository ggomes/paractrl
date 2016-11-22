function addslide_with_text(op,titletext,bodytext,neworlast,pos,prnopt)
% addslide(op,titletext,neworlast,pos,prnopt)

% Capture current figure/model into clipboard:
if nargin<6
    print -dmeta
else
    print('-dmeta',prnopt)
end

if nargin<5
    pos = [0.5 0.5];
end

if nargin<4
    neworlast = 'new';
end

% Get current number of slides:
slide_count = get(op.Slides,'Count');

% Add a new slide (with title object):
if strcmp(neworlast,'new')
    slide_count = int32(double(slide_count)+1);
    slide = invoke(op.Slides,'Add',slide_count,11);
else
    slide_count = int32(double(slide_count));
    slide=invoke(invoke(op.Slides,'Range'),'Item',slide_count);
end

% Insert text into the title object:
if(~isempty(titletext) && strcmp(neworlast,'new') )
    set(slide.Shapes.Title.TextFrame.TextRange,'Text',titletext);
end

% Get height and width of slide:
slide_H = op.PageSetup.SlideHeight;
slide_W = op.PageSetup.SlideWidth;

% Paste the contents of the Clipboard:
pic1 = invoke(slide.Shapes,'Paste');

% Get height and width of picture:
pic_H = get(pic1,'Height');
pic_W = get(pic1,'Width');

% Center picture on page (below title area):
set(pic1,'Left',single(max([pos(1)*double(slide_W)-0.5*double(pic_W) 0])));
set(pic1,'Top' ,single(max([(1-pos(2))*double(slide_H)-0.5*double(pic_H) 0])));

% Insert body text
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