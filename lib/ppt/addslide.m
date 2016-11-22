function addslide(op,titletext,neworlast,pos,prnopt,scale)
% addslide(op,titletext,neworlast,pos,prnopt,scale)

% Capture current figure/model into clipboard:
if nargin<6
    scale = 1;
end
if nargin<5 || isempty(prnopt)
    print -dbitmap
else
    print('-dbitmap',prnopt)
end

if nargin<4 || isempty(pos)
    pos = [0.5 0.5];
end

if nargin<3 || isempty(neworlast)
    neworlast = 'new';
end

% Get current number of slides:
slide_count = get(op.get('Slides'),'Count');

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

set(pic1,'Height',scale*pic_H);
set(pic1,'Width',scale*pic_W);
pic_H = get(pic1,'Height');
pic_W = get(pic1,'Width');

% Center picture on page (below title area):
set(pic1,'Left',single(max([pos(1)*double(slide_W)-0.5*double(pic_W) 0])));
set(pic1,'Top' ,single(max([(1-pos(2))*double(slide_H)-0.5*double(pic_H) 0])));

return