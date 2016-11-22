classdef CCTypes
    
    methods(Static)
        
        function A = get_type_table(tablename)
            A = CCTypes.read_table_file(tablename);
        end
         
        function p = get_type_with_name(tablename,typename)
            A = CCTypes.read_table_file(tablename);
            ind = strcmp(typename,{A.name});
            if(any(ind))
                p=A(ind);
            else
                p=[];
            end
        end
        
        function p = get_type_with_id(tablename,typeid)
            A = CCTypes.read_table_file(tablename);
            ind = typeid==[A.id];
            if(any(ind))
                p=A(ind);
            else
                p=[];
            end
        end
        
        function [A]=read_table_file(tablename)                  
            filename = fullfile(fileparts(mfilename('fullpath')),'db_printout',[tablename '_types.txt']);
            fid = fopen(filename);
            if(fid<0)
                error('bad type name')
            end
            data = textscan(fid,'%d%s%s','delimiter','\t');
            fclose(fid);
            A=cell2struct([num2cell(data{1}) data{2} data{3}]',{'id','name','description'});
        end
        
    end
   
end