function compile_xsd(infile,outfile)


infile = 'C:\Users\gomes\code\L0\L0-utilities\L1_specs\mo\xsd\scenario.xsd';
infolder = fileparts(infile);
outfile = 'C:\Users\gomes\code\L0\L0-utilities\L1_specs\mo\xsd\bla.xml';

outfid = fopen(outfile,'w+');
include_xsd(infolder,infile,outfid,true)
fclose(outfid);

function []=include_xsd(infolder,filename,outfid,includeheader)

fid=fopen(filename);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if(~isempty(strfind(tline,'<?xml')) && ~includeheader)
        continue
    end
    
    if(~isempty(strfind(tline,'<xs:schema')) && ~includeheader)
        continue
    end
    
    if(~isempty(strfind(tline,'</xs:schema')) && ~includeheader)
        continue
    end
    
    if(~isempty(strfind(tline,'<xs:include')))
        z=regexp(tline,'"');
        includefile = fullfile(infolder,tline(z(1)+1:z(2)-1));
        include_xsd(infolder,includefile,outfid,false)
    else
        fwrite(outfid,sprintf('%s\n',tline));
    end
end
fclose(fid);
