classdef ScenarioPtr < handle
    
    properties (Access = public)
        scenario
        link_id_begin_end
    end
    
    methods (Access = public)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% load/save/clone
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [obj,success]=load(obj,configfile,validate,silent)
            % load scenario from an xml file.
            
            if(nargin<3)
                validate = true;
            end
            
            if(nargin<4)
                silent = false;
            end
            
            if(~exist(configfile,'file'))
                error('file not found')
            end
                        
            success = false;
            
            % validate
            if(validate)
                utils_path = fileparts(fileparts(mfilename('fullpath')));
                xml_validator = fullfile(utils_path,'lib','xml-validation','xsd11-validator.jar');
                schemafile = fullfile(utils_path,'simulation','beats.xsd');
                % schemafile = fullfile(utils_path,'L1_specs','mo','xsd','scenario_full.xsd');
                [status,result] = system(['java -jar "' xml_validator '" -sf "' schemafile '" -if "' configfile '"']);
                
                if(~isempty(result) || status~=0 )
                    disp('Input file not valid')
                    disp(result)
                    return
                end
            end
            
            if ~silent
                disp(['Loading ' configfile])
            end
            obj.scenario = xml_read(configfile);
            
            % check single network
            if obj.has_network && length(obj.scenario.NetworkSet.network)>1
                error('ScenarioPtr does not support scenarios with multiple networks')
            end
            
            % change xEnd to end
            if(obj.has_links)
                numlinks = length(obj.scenario.NetworkSet.network.LinkList.link);
                newLink = safe_rmfield(generate_mo('link',true),{'roads','dynamics'});
                newLink.ATTRIBUTE = safe_rmfield(newLink.ATTRIBUTE,{'mod_stamp','link_name','speed_limit','crudFlag','lane_offset','in_sync'});
                newLinks = repmat(newLink,1,numlinks);
                for j=1:numlinks
                    L = obj.scenario.NetworkSet.network.LinkList.link(j);
                    if(isfield(L,'position'))
                        newLinks(j).position = L.position;
                    end
                    if(isfield(L,'begin'))
                        newLinks(j).begin = L.begin;
                    end
                    if(isfield(L,'ATTRIBUTE'))
                        newLinks(j).ATTRIBUTE = L.ATTRIBUTE;
                    end
                    if(isfield(L,'xEnd'))
                        newLinks(j).end = L.xEnd;
                    end
                    if(isfield(L,'shape'))
                        newLinks(j).shape = L.shape;
                    end
                    if(isfield(L,'link_type'))
                        newLinks(j).link_type = L.link_type;
                    end
                end
                obj.scenario.NetworkSet.network.LinkList.link = newLinks;
            end
            
            % populate link_inputs and link_outputs
            obj.generate_link_id_begin_end();
            
            success = true;
            
        end
        
        function [obj]=load_fwy_from_excel(obj,excelfile)
            obj.scenario = FwyScenarioLoader.load_scenario_from_excel(excelfile);
            obj.generate_link_id_begin_end;
        end
        
        function [obj]=save(obj,outfile)

            
            disp(['Saving ' outfile])
            
            % copy scenario
            scenario = obj.scenario;
            uts = lower(obj.get_units());
            switch uts
                case 'us'
                    precision.splits = '%.5f';
                    precision.demands = '%.3f';
                case 'si'
                    precision.splits = '%.5f';
                    precision.demands = '%.7f';
                case 'metric'
                    precision.splits = '%.5f';
                    precision.demands = '%.3f';
                otherwise
                    precision.splits = '%.5f';
                    precision.demands = '%.3f';
            end
            
            % remove nans and empties from...
            
            % ... point.elevation, point.lat, point.lng
            if(obj.has_nodes && obj.has_links)
                scenario.NetworkSet.network = obj.replace_nan_with(scenario.NetworkSet.network,{'position','point','ATTRIBUTE','elevation'},0);
                for j=1:length(scenario.NetworkSet.network.NodeList.node)
                    scenario.NetworkSet.network.NodeList.node(j) = ...
                        obj.replace_nan_with(scenario.NetworkSet.network.NodeList.node(j),...
                        {'position','point','ATTRIBUTE','elevation'},0);
                    scenario.NetworkSet.network.NodeList.node(j) = ...
                        obj.replace_empty_with(scenario.NetworkSet.network.NodeList.node(j),...
                        {'position','point','ATTRIBUTE','elevation'},0);
                end
                for j=1:length(scenario.NetworkSet.network.LinkList.link)
                    scenario.NetworkSet.network.LinkList.link(j) = ...
                        obj.replace_nan_with(scenario.NetworkSet.network.LinkList.link(j),...
                        {'position','point','ATTRIBUTE','elevation'},0);
                    scenario.NetworkSet.network.LinkList.link(j) = ...
                        obj.replace_empty_with(scenario.NetworkSet.network.LinkList.link(j),...
                        {'position','point','ATTRIBUTE','elevation'},0);
                    scenario.NetworkSet.network.LinkList.link(j) = ...
                        obj.replace_nan_with(scenario.NetworkSet.network.LinkList.link(j),...
                        {'position','point','ATTRIBUTE','lat'},0);
                    scenario.NetworkSet.network.LinkList.link(j) = ...
                        obj.replace_empty_with(scenario.NetworkSet.network.LinkList.link(j),...
                        {'position','point','ATTRIBUTE','lat'},0);
                    scenario.NetworkSet.network.LinkList.link(j) = ...
                        obj.replace_nan_with(scenario.NetworkSet.network.LinkList.link(j),...
                        {'position','point','ATTRIBUTE','lng'},0);
                    scenario.NetworkSet.network.LinkList.link(j) = ...
                        obj.replace_empty_with(scenario.NetworkSet.network.LinkList.link(j),...
                        {'position','point','ATTRIBUTE','lng'},0);
                end
                
                
                % set network position
                if(isfieldRecursive(scenario.NetworkSet.network,'position','point'))
                    scenario.NetworkSet.network.position.point(2:end)=[];
                else
                    points = repmat(struct('lat',nan,'lng',nan),1,length(scenario.NetworkSet.network.NodeList.node));
                    for j=1:length(scenario.NetworkSet.network.NodeList.node)
                        if(isfield(scenario.NetworkSet.network.NodeList.node(j),'position'))
                            p = scenario.NetworkSet.network.NodeList.node(j).position.point.ATTRIBUTE;
                            points(j).lat = p.lat;
                            points(j).lng = p.lng;
                        end
                    end
                    scenario.NetworkSet.network.position.point.ATTRIBUTE.elevation = 0;
                    scenario.NetworkSet.network.position.point.ATTRIBUTE.lat = mean([points.lat]);
                    scenario.NetworkSet.network.position.point.ATTRIBUTE.lng = mean([points.lng]);
                end
            end
            
            % ... sensors
            if(obj.has_sensors)
                for i=1:length(scenario.SensorSet.sensor)
                    S = scenario.SensorSet.sensor(i);
                    S.ATTRIBUTE = obj.remove_if_nan(S.ATTRIBUTE,{'health_status','lane_number'});
                    S = obj.replace_nan_with(S,{'display_position','point','ATTRIBUTE','elevation'},0);
                    scenario.SensorSet.sensor(i)=S;
                end
            end
            
            % ... jam densities
            if(obj.has_fds)
                for i=1:length(scenario.FundamentalDiagramSet.fundamentalDiagramProfile)
                    fd = scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i);
                    for j=1:length(fd.fundamentalDiagram)
                        fd.fundamentalDiagram(j).ATTRIBUTE=obj.remove_if_nan(fd.fundamentalDiagram(j).ATTRIBUTE,{'jam_density'});
                    end
                    scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i) = fd;
                end
            end
            
            % ... format demands
            if(obj.has_demands)
                for i=1:length(scenario.DemandSet.demandProfile)
                    dP = scenario.DemandSet.demandProfile(i);
                    if(isfield(dP,'demand'))
                        for j=1:length(dP.demand)
                            if(ischar(dP.demand(j).CONTENT))
                                x = str2double(dP.demand(j).CONTENT);
                                if(isnan(x))
                                    x=eval(['[' dP.demand(j).CONTENT ']']);
                                end
                            else
                                x = dP.demand(j).CONTENT;
                            end
                            scenario.DemandSet.demandProfile(i).demand(j).CONTENT = writecommaformat(x,precision.demands);
                        end
                    end
                end
            end
            
            % ... format split ratios
            if(obj.has_splits)
                for i=1:length(scenario.SplitRatioSet.splitRatioProfile)
                    sr = scenario.SplitRatioSet.splitRatioProfile(i);
                    for j=1:length(sr.splitratio)
                        if(isfield(sr.splitratio(j),'CONTENT'))
                            if(ischar(sr.splitratio(j).CONTENT))
                                x = str2double(sr.splitratio(j).CONTENT);
                                if(isnan(x))
                                    x=eval(['[' sr.splitratio(j).CONTENT ']']);
                                end
                            else
                                x = sr.splitratio(j).CONTENT;
                            end
                            scenario.SplitRatioSet.splitRatioProfile(i).splitratio(j).CONTENT = writecommaformat(x,precision.splits);
                        end
                        if(isfield(sr,'concentrationParameters')) && (~isempty(sr.concentrationParameters))
                            if(isfield(sr.concentrationParameters(j),'CONTENT'))
                                if(ischar(sr.splitratio(j).CONTENT))
                                    x = str2double(sr.concentrationParameters(j).CONTENT);
                                else
                                    x = sr.concentrationParameters(j).CONTENT;
                                end
                                scenario.SplitRatioSet.splitRatioProfile(i).concentrationParameters(j).CONTENT = writecommaformat(x,precision.demands);
                            end
                        end
                    end
                end
            end
            
            % add project_id and id to sets
            sets = {'NetworkSet','EventSet','ControllerSet'};
            for i=1:length(sets)
                if(isfield(scenario,sets{i}))
                    if(~isfield(scenario.(sets{i}),'ATTRIBUTE'))
                        X = generate_mo(sets{i},true);
                        scenario.(sets{i}).ATTRIBUTE = X.ATTRIBUTE;
                        clear X
                    end
                    scenario.(sets{i}).ATTRIBUTE = obj.default_to(scenario.(sets{i}).ATTRIBUTE,'project_id',0);
                    scenario.(sets{i}).ATTRIBUTE = obj.default_to(scenario.(sets{i}).ATTRIBUTE,'id',0);
                end
            end
            
            xml_write(outfile,scenario);
            clear scenario
            
            %%% POST PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fid_in = fopen(outfile,'r');
            fname = [num2str(round(rand*10000000000)) '.xml'];
            fid_out = fopen(fname,'w+');
            while 1
                tline = fgetl(fid_in);
                if ~ischar(tline), break, end
                
                % remove empty link references from sensors
                if(~isempty(strfind(tline,'sensor')) && ~isempty(strfind(tline,'link_id=""')))
                    tline = strrep(tline,'link_id=""', '');
                end
                
                % remove empty sensor reference from fd
                if(~isempty(strfind(tline,'fundamentalDiagramProfile')) && ~isempty(strfind(tline,'sensor_id=""')))
                    tline = strrep(tline,'sensor_id=""', '');
                end
                
                % remove empty positions
                if(~isempty(strfind(tline,'<position/>')))
                    continue
                end
                
                % remove empty shapes
                if(~isempty(strfind(tline,'<shape/>')))
                    continue
                end
                
                fprintf(fid_out,'%s\n',tline);
            end
            fclose(fid_in);
            fclose(fid_out);
            if ispc
                system(['move ' fname ' "' outfile '"']);
            elseif isunix || ismac
                system(['mv ' fname ' "' outfile '"']);
            end
            
        end
        
        function [out] = clone(obj)
            out = ScenarioPtr;
            out.scenario = obj.scenario;
            out.link_id_begin_end = obj.link_id_begin_end;
        end
        
        % generate link/node map ..........................................
        function [] = generate_link_id_begin_end(obj)
            % populate link_inputs and link_outputs
            numlinks = obj.get_num_links;
            obj.link_id_begin_end = nan(numlinks,3);
            for i=1:numlinks
                L=obj.scenario.NetworkSet.network.LinkList.link(i);
                obj.link_id_begin_end(i,:) = [  L.ATTRIBUTE.id ...
                    L.begin.ATTRIBUTE.node_id ...
                    L.end.ATTRIBUTE.node_id];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [errors]=check_network(obj)
            
            errors = [];
            
            % nodelist node and link ids
            x = [obj.scenario.NetworkSet.network.NodeList.node];
            x = [x.outputs];
            x = [x.output];
            x = [x.ATTRIBUTE];
            nodelist.link_ids = [x.link_id];
            x = [obj.scenario.NetworkSet.network.NodeList.node];
            x = [x.ATTRIBUTE];
            nodelist.node_ids = [x.id];
            
            % linklist node and link ids
            linklist.link_ids = obj.link_id_begin_end(:,1)';
            linklist.node_ids = unique([obj.link_id_begin_end(:,2);obj.link_id_begin_end(:,3)])';
            
            % ids are unique
            if(length(nodelist.node_ids)~=length(unique(nodelist.node_ids)))
                errors = [errors {'node ids are not unique'}];
            end
            if(length(linklist.link_ids)~=length(unique(linklist.link_ids)))
                errors = [errors {'link ids are not unique'}];
            end
            
            % compare node ids
            node_minus_link = setdiff(nodelist.node_ids,linklist.node_ids);
            if(~isempty(node_minus_link))
                errors = [errors {['node ids in nodelist and not linklist: ' num2str(node_minus_link)]}];
            end
            link_minus_node = setdiff(linklist.node_ids,nodelist.node_ids);
            if(~isempty(link_minus_node))
                errors = [errors {['node ids in linklist and not nodelist: ' num2str(link_minus_node)]}];
            end
            
            % compare link ids
            node_minus_link = setdiff(nodelist.link_ids,linklist.link_ids);
            if(~isempty(node_minus_link))
                errors = [errors {['link ids in nodelist and not linklist: ' num2str(node_minus_link)]}];
            end
            link_minus_node = setdiff(linklist.link_ids,nodelist.link_ids);
            if(~isempty(link_minus_node))
                errors = [errors {['link ids in linklist and not nodelist: ' num2str(link_minus_node)]}];
            end
            
            % check that node io according to nodelist agrees with linklist
            bad_nodes = [];
            for i=1:length(obj.scenario.NetworkSet.network.NodeList.node)
                node = obj.scenario.NetworkSet.network.NodeList.node(i);
                
                node_id = node.ATTRIBUTE.id;
                
                % get output ids
                if(isfieldRecursive(node,'outputs','output'))
                    x=[node.outputs.output.ATTRIBUTE];
                    out_ids = [x.link_id];
                else
                    out_ids = [];
                end
                
                % get input ids
                if(isfieldRecursive(node,'inputs','input'))
                    x=[node.inputs.input.ATTRIBUTE];
                    in_ids = [x.link_id];
                else
                    in_ids = [];
                end
                
                % according to link list:
                out_links = obj.link_id_begin_end(obj.link_id_begin_end(:,2)==node_id,1);
                in_links = obj.link_id_begin_end(obj.link_id_begin_end(:,3)==node_id,1);
                
                % are they the same?
                if(~ScenarioPtr.arrays_are_same(out_ids,out_links))
                    bad_nodes = [bad_nodes struct( 'node_id' , node_id , ...
                        'type','out',...
                        'node_list',out_ids,...
                        'link_list',out_links)];
                    
                end
                
                if(~ScenarioPtr.arrays_are_same(in_ids,in_links))
                    bad_nodes = [bad_nodes struct( 'node_id' , node_id , ...
                        'type','in',...
                        'node_list',in_ids,...
                        'link_list',in_links)];
                end
                
            end
            
            if(~isempty(bad_nodes))
                errors = [errors {bad_nodes}]';
            end
            
            if(isempty(errors))
                disp('Network passed')
            end
            
        end
        
        function []=fix_nodelist(obj)
            
            
            % nodelist node and link ids
            x = [obj.scenario.NetworkSet.network.NodeList.node];
            x = [x.outputs];
            x = [x.output];
            x = [x.ATTRIBUTE];
            nodelist.link_ids = [x.link_id];
            x = [obj.scenario.NetworkSet.network.NodeList.node];
            x = [x.ATTRIBUTE];
            nodelist.node_ids = [x.id];
            
            % linklist node and link ids
            linklist.link_ids = obj.link_id_begin_end(:,1)';
            linklist.node_ids = unique([obj.link_id_begin_end(:,2);obj.link_id_begin_end(:,3)])';
            
            % remove nodes that are in nodelist but not linklist
            remove_nodes = setdiff(nodelist.node_ids,linklist.node_ids);
            obj.scenario.NetworkSet.network.NodeList.node(ismember(remove_nodes,obj.get_node_ids)) = [];
            
            % add nodes that are in linklist but node nodelist
            add_nodes = setdiff(linklist.node_ids,nodelist.node_ids);
            if(~isempty(add_nodes))
                error('not implemented!')
            end
            
            % regenerate node inputs and outputs
            for i=1:length(obj.scenario.NetworkSet.network.NodeList.node)
                node = obj.scenario.NetworkSet.network.NodeList.node(i);
                node_id = node.ATTRIBUTE.id;
                in = obj.link_id_begin_end(obj.link_id_begin_end(:,3)==node_id,1);
                if(~isempty(in))
                    x = repmat(generate_mo('input'),1,length(in));
                    for j=1:length(in)
                        x(j).ATTRIBUTE.link_id = in(j);
                    end
                    node.inputs.input = x;
                else
                    node.inputs = [];
                end
                out = obj.link_id_begin_end(obj.link_id_begin_end(:,2)==node_id,1);
                if(~isempty(out))
                    x = repmat(generate_mo('output'),1,length(out));
                    for j=1:length(out)
                        x(j).ATTRIBUTE.link_id = out(j);
                    end
                    node.outputs.output = x;
                else
                    node.outputs = [];
                end
                obj.scenario.NetworkSet.network.NodeList.node(i) = node;
            end
            
            
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% get/set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % settings and vehicle types ......................................
        
        function [u] = get_units(obj)
            u='';
            if(~isfield(obj.scenario,'settings'))
                return
            end
            if(~isfield(obj.scenario.settings,'units'))
                return
            end
            u = lower(obj.scenario.settings.units);
        end
        
        function [vt] = get_vehicle_types(obj)
            vt = [];
            if(~isfield(obj.scenario,'VehicleTypeSet'))
                return
            end
            if(~isfield(obj.scenario.VehicleTypeSet,'vehicleType'))
                return
            end
            numvt = length(obj.scenario.VehicleTypeSet.vehicleType);
            vt = repmat(struct('id',[],'name',''),1,numvt);
            for i=1:numvt
                x = obj.scenario.VehicleTypeSet.vehicleType(i);
                vt(i).id = x.ATTRIBUTE.id;
                vt(i).name = x.ATTRIBUTE.name;
            end
        end
        
        function [obj] = change_units_to(obj,newunits)
            
            if(strcmpi(obj.get_units,newunits))
                return
            end
            
            cnv = unit_conversion_factor(obj.get_units(),newunits);
            
            % initial density
            if(obj.has_initial_densities)
                error('NEED TO IMPLEMENT THIS!')
            end
            
            % downstream capacity profiles
            if(obj.has_capacity_profiles)
                error('NEED TO IMPLEMENT THIS!')
            end
            
            % links
            for i=1:obj.get_num_links
                L = obj.scenario.NetworkSet.network.LinkList.link(i);
                L = obj.safe_conv(L,cnv.length,'ATTRIBUTE','length');
                obj.scenario.NetworkSet.network.LinkList.link(i) = L;
            end
            
            % FDs
            if(obj.has_fds)
                for i=1:numel(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile)
                    for j=1:numel(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i).fundamentalDiagram)
                        FD = obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i).fundamentalDiagram(j);
                        FD = obj.safe_conv(FD,cnv.flw,'ATTRIBUTE','capacity');
                        FD = obj.safe_conv(FD,cnv.flw,'ATTRIBUTE','capacity_drop');
                        FD = obj.safe_conv(FD,cnv.flw,'ATTRIBUTE','std_dev_capacity');
                        FD = obj.safe_conv(FD,cnv.spd,'ATTRIBUTE','free_flow_speed');
                        FD = obj.safe_conv(FD,cnv.spd,'ATTRIBUTE','congestion_speed');
                        FD = obj.safe_conv(FD,cnv.spd,'ATTRIBUTE','critical_speed');
                        FD = obj.safe_conv(FD,cnv.spd,'ATTRIBUTE','std_dev_free_flow_speed');
                        FD = obj.safe_conv(FD,cnv.spd,'ATTRIBUTE','std_dev_congestion_speed');
                        FD = obj.safe_conv(FD,cnv.dty,'ATTRIBUTE','jam_density');
                        FD = obj.safe_conv(FD,cnv.dty,'ATTRIBUTE','transition_density');
                        obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i).fundamentalDiagram(j) = FD;
                    end
                end
            end
            
            % demands
            if(obj.has_demands)
                for i=1:numel(obj.scenario.DemandSet.demandProfile)
                    dP = obj.scenario.DemandSet.demandProfile(i);
                    for j=1:length(dP.demand)
                        dP.demand(j) = obj.safe_conv(dP.demand(j),cnv.flw,'CONTENT');
                    end
                    obj.scenario.DemandSet.demandProfile(i) = dP;
                end
            end
            
            obj.scenario.settings.units = newunits;
            
        end
        
        % make single vehicle type
        function [] = change_to_single_vtype(obj, vtype_id)
            
            if numel(obj.scenario.VehicleTypeSet.vehicleType) == 1
                return
            end
            
            if nargin < 2 || isempty(vtype_id)
                vtype_id = obj.scenario.VehicleTypeSet.vehicleType(1).ATTRIBUTE.id;
            end
            
            % remove excess vehicle types
            temp = [obj.scenario.VehicleTypeSet.vehicleType(:).ATTRIBUTE];
            vtype_idx = [temp.id] == vtype_id;
            obj.scenario.VehicleTypeSet.vehicleType( ~vtype_idx) = [];
            
            clear temp
            
            % remove demands for extra vehicle types
            % assumes all demands in one demandProfile have same dt and
            % start/end time
            for i_dp = 1 : numel(obj.scenario.DemandSet.demandProfile)
                unused_demand = zeros(size(obj.scenario.DemandSet.demandProfile( i_dp ).demand(1).CONTENT));
                % iterate through the demands and add up the ones for extra
                % vehicles
                unsetter_array = false(1, numel( obj.scenario.DemandSet.demandProfile( i_dp ).demand));
                for i_dem = 1 : numel(obj.scenario.DemandSet.demandProfile( i_dp ).demand)
                    if obj.scenario.DemandSet.demandProfile( i_dp ).demand( i_dem ).ATTRIBUTE.vehicle_type_id ~= vtype_id
                        unused_demand = unused_demand + obj.scenario.DemandSet.demandProfile( i_dp ).demand(i_dem).CONTENT;
                        unsetter_array( i_dem ) = true;
                    else
                        this_dp_vtype_idx = i_dem;
                    end
                end
                obj.scenario.DemandSet.demandProfile(i_dp).demand( this_dp_vtype_idx ).CONTENT = ...
                    obj.scenario.DemandSet.demandProfile(i_dp).demand( this_dp_vtype_idx ).CONTENT + unused_demand;
                
                % unset the demands for this dp that are of the extra
                % vehicle types
                obj.scenario.DemandSet.demandProfile(i_dp).demand( unsetter_array ) = [];
            end
            
            % remove splits for extra vehicle types
            for i_srp = 1:numel(obj.scenario.SplitRatioSet.splitRatioProfile)
                unsetter_array = false(1, numel( obj.scenario.SplitRatioSet.splitRatioProfile( i_dp ).splitratio));
                for i_splitr = 1:numel(obj.scenario.SplitRatioSet.splitRatioProfile(i_srp).splitratio)
                    if obj.scenario.SplitRatioSet.splitRatioProfile(i_srp).splitratio(i_splitr).ATTRIBUTE.vehicle_type_id ~= vtype_id
                        unsetter_array( i_splitr ) = true;
                    end
                end
                obj.scenario.SplitRatioSet.splitRatioProfile(i_srp).splitratio(unsetter_array) = [];
            end
            
        end
        
        function [] = remove_hov_links(obj)
            
            if numel(obj.scenario.VehicleTypeSet.vehicleType) ~= 1
                fprintf('Method remove_hov_links only supported for single vtype scenarios\n');
                return
            end
            
            links_per_network = numel(obj.scenario.NetworkSet.network.LinkList.link);
            nodes_per_network = numel(obj.scenario.NetworkSet.network.NodeList.node);
            
            links_offset = cumsum(links_per_network);
            nodes_offset = cumsum(nodes_per_network);
            
            links_removed = [];
            nodes_removed = [];
            
            offset_link_inds = 1:links_offset(1);
            offset_node_inds = 1:nodes_offset(1);
            
            linkids = obj.get_link_ids;
            linktypes = obj.get_link_types;
            linktypes = linktypes(offset_link_inds);
            links_to_remove = strcmpi(linktypes,'hov')';
            link_ids_to_remove = linkids(links_to_remove);
            
            nodeids = obj.get_node_ids;
            nodeids_this_network = nodeids(offset_node_inds);
            
            % collect the list of nodes to remove (ie, nodes that only
            % connect hov links)
            nodes_to_remove = false(numel(obj.scenario.NetworkSet.network.NodeList.node),1);
            for i_node = 1:numel(nodes_to_remove)
                node = obj.scenario.NetworkSet.network.NodeList.node(i_node);
                links_with_this_node=[];
                outputs_to_remove = [];
                if ~isempty(node.outputs)
                    for i_output = 1:numel(node.outputs.output)
                        links_with_this_node(end+1) = node.outputs.output(i_output).ATTRIBUTE.link_id;
                    end
                    [~, outputs_to_remove, ~] = intersect(links_with_this_node,link_ids_to_remove);
                end
                numoutputs = numel(links_with_this_node);
                inputs_to_remove = [];
                if ~isempty(node.inputs)
                    for i_input = 1:numel(node.inputs.input)
                        links_with_this_node(end+1) = node.inputs.input(i_input).ATTRIBUTE.link_id;
                    end
                    [~, inputs_to_remove, ~] = intersect(links_with_this_node(numoutputs+1:end),link_ids_to_remove);
                end
                if numel([outputs_to_remove, inputs_to_remove]) == numel(links_with_this_node) % if every link attached to this node is hov-type and going away
                    nodes_to_remove(i_node) = true;
                end
                if ~isempty(outputs_to_remove)
                    obj.scenario.NetworkSet.network.NodeList.node(i_node).outputs.output(outputs_to_remove) = [];
                end
                if ~isempty(inputs_to_remove)
                    obj.scenario.NetworkSet.network.NodeList.node(i_node).inputs.input(inputs_to_remove) = [];
                end
            end % i_node
            
            node_ids_to_remove = nodeids_this_network(nodes_to_remove);
            
            % removing stuff here
            obj.scenario.NetworkSet.network.LinkList.link(links_to_remove) = [];
            obj.scenario.NetworkSet.network.NodeList.node(nodes_to_remove) = [];
            
            links_removed = vertcat(links_removed, link_ids_to_remove);
            nodes_removed = vertcat(nodes_removed, node_ids_to_remove);
            
            % splits
            numsplitprofiles = numel(obj.scenario.SplitRatioSet.splitRatioProfile);
            splitprofiles_to_remove = false(numsplitprofiles,1);
            for i_split = 1:numsplitprofiles
                splitprofile = obj.scenario.SplitRatioSet.splitRatioProfile(i_split);
                node_for_this_split = splitprofile.ATTRIBUTE.node_id;
                if ismember(node_for_this_split, nodes_removed)
                    splitprofiles_to_remove(i_split) = true;
                    continue;
                end
                splitdata = [ splitprofile.splitratio(:).ATTRIBUTE];
                in_links = [splitdata.link_in];
                hov_in_links = ismember(in_links, links_removed);
                obj.scenario.SplitRatioSet.splitRatioProfile(i_split).splitratio(hov_in_links) = [];
                
                splitprofile = obj.scenario.SplitRatioSet.splitRatioProfile(i_split); % refresh it after deleting the hov-in SRs
                splitdata = [ splitprofile.splitratio(:).ATTRIBUTE];
                
                in_links = [splitdata.link_in];
                out_links = [splitdata.link_out];
                hov_out_links = ismember(out_links, links_removed);
                
                unique_hov_out_link_ids = unique(out_links(hov_out_links));
                num_unique_hov_out_links = numel(unique_hov_out_link_ids);
                if num_unique_hov_out_links < 0
                    continue;
                end
                split_to_remove_cum = false(size(out_links)); % initialization of hov-involved splits that will be pruned after reassigning their split numbers
                for i_linkout = 1:num_unique_hov_out_links
                    this_out_link_id = unique_hov_out_link_ids(i_linkout);
                    inlinks_with_this_out = unique( in_links( out_links == this_out_link_id ) );
                    for i_linkin = 1:numel(inlinks_with_this_out)
                        this_in_link_id = inlinks_with_this_out(i_linkin);
                        split_to_remove = and( in_links == this_in_link_id, out_links == this_out_link_id );
                        if nnz(split_to_remove) ~= 1
                            error('Something is odd. Splitratioprofile number %g may have duplicate splits', i_split);
                        end
                        splits_to_divide_among = and( in_links == this_in_link_id, out_links ~= this_out_link_id );
                        constant = 1/nnz(splits_to_divide_among);
                        split_inds_to_divide_among = find( splits_to_divide_among );
                        for i_reallocate = 1:numel(split_inds_to_divide_among)
                            obj.scenario.SplitRatioSet.splitRatioProfile(i_split).splitratio(split_inds_to_divide_among(i_reallocate)).CONTENT = ...
                                obj.scenario.SplitRatioSet.splitRatioProfile(i_split).splitratio(split_inds_to_divide_among(i_reallocate)).CONTENT + ...
                                constant * obj.scenario.SplitRatioSet.splitRatioProfile(i_split).splitratio(split_to_remove).CONTENT;
                        end
                        split_to_remove_cum = or(split_to_remove_cum, split_to_remove);
                    end
                end
                obj.scenario.SplitRatioSet.splitRatioProfile(i_split).splitratio(split_to_remove_cum) = [];
            end %i_split
            
            % FDs
            num_fds = numel(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile);
            fds_to_remove = false(1,num_fds);
            for i_fd = 1:num_fds
                if ismember(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i_fd).ATTRIBUTE.link_id, links_removed)
                    fds_to_remove(i_fd) = true;
                end
            end
            obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(fds_to_remove) = [];
            
        end % function end
        
        % networks ........................................................
        
        function OneOrZero = get_num_networks(obj)
            if obj.has_network
                OneOrZero = 1;
            else
                OneOrZero = 0;
            end
        end
        
        function [b] = has_network(obj)
            b = isfieldRecursive(obj.scenario,'NetworkSet','network');
        end
        
        function [] = add_link(obj,newLink)
            if(length(newLink)>1)
                error('add links one at a time')
            end
            
            % check that the start and end nodes are present
            start_node_id = newLink.begin.ATTRIBUTE.node_id;
            end_node_id = newLink.end.ATTRIBUTE.node_id;
            node_ids = obj.get_node_ids;
            if(~ismember(start_node_id,node_ids) || ~ismember(end_node_id,node_ids))
                error('unknown start or end node ')
            end
            if(~obj.has_links)
                obj.scenario.NetworkSet.network.LinkList.link = [];
            end
            
            try
                obj.scenario.NetworkSet.network.LinkList.link = [...
                    obj.scenario.NetworkSet.network.LinkList.link ...
                    newLink];
            catch
                error('new links could not be added to the network')
            end
            
            obj.link_id_begin_end = [obj.link_id_begin_end ; ...
                [newLink.ATTRIBUTE.id start_node_id end_node_id] ];
        end
        
        function [] = add_links(obj,links)
            for i=1:length(links)
                obj.add_link(links(i))
            end
        end
        
        function [] = add_node(obj,newNode)
            if(length(newNode)>1)
                error('add nodes one at a time')
            end
            if(~obj.has_nodes)
                obj.scenario.NetworkSet.network.NodeList.node =[];
            end
            try
                obj.scenario.NetworkSet.network.NodeList.node = [...
                    obj.scenario.NetworkSet.network.NodeList.node ...
                    newNode];
            catch
                error('new nodes could not be added to the network')
            end
        end
        
        function [] = add_nodes(obj,newNodes)
            for i=1:length(newNodes)
                obj.add_node(newNodes(i))
            end
        end
        
        function obj = extract_subnetwork(obj,linksToKeep)
            linkIds = obj.get_link_ids;
            linkIndecesToPrune = ~ismember(linkIds,linksToKeep);
            linksToPrune = linkIds(linkIndecesToPrune);
            
            obj.scenario.NetworkSet.network.LinkList.link(linkIndecesToPrune) = [];
            
            nodesRemoved = [];
            nodeIndecesToRemove = [];
            for i = 1:numel(obj.scenario.NetworkSet.network.NodeList.node)
                node = obj.scenario.NetworkSet.network.NodeList.node(i);
                if ~isempty(node.outputs) 
                    outputs = [node.outputs.output.ATTRIBUTE];
                    if ~isempty(outputs)
                        outputLinks = [outputs.link_id];
                        outputsRemoved = ismember(outputLinks,linksToPrune);
                        obj.scenario.NetworkSet.network.NodeList.node(i).outputs.output(outputsRemoved) = [];
                    end
                else
                    outputsRemoved = [];
                end
                if ~isempty(node.inputs)
                    inputs = [node.inputs.input.ATTRIBUTE];
                    if ~isempty(inputs)
                        inputLinks = [inputs.link_id];
                        inputsRemoved = ismember(inputLinks,linksToPrune);
                        obj.scenario.NetworkSet.network.NodeList.node(i).inputs.input(inputsRemoved) = [];
                    end
                else
                    inputsRemoved = [];
                end
                
                if all(outputsRemoved) && all(inputsRemoved)
                    nodesRemoved = vertcat(nodesRemoved, obj.scenario.NetworkSet.network.NodeList.node(i).ATTRIBUTE.id);
                    nodeIndecesToRemove = vertcat(nodeIndecesToRemove, i );
                end
            end
            obj.scenario.NetworkSet.network.NodeList.node(nodeIndecesToRemove) = [];
            
            if(isfieldRecursive(obj.scenario,'DemandSet','demandProfile'))
                demandProfilesRemoved = [];
                for i = 1:numel(obj.scenario.DemandSet.demandProfile)
                    if ismember(obj.scenario.DemandSet.demandProfile(i).ATTRIBUTE.link_id_org, linksToPrune)
                        demandProfilesRemoved = vertcat(demandProfilesRemoved, i);
                    end
                end
                obj.scenario.DemandSet.demandProfile(demandProfilesRemoved) = [];
            end
            
            if(isfieldRecursive(obj.scenario,'FundamentalDiagramSet','fundamentalDiagramProfile'))
                FDsRemoved = [];
                for i = 1:numel(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile)
                    if ismember(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i).ATTRIBUTE.link_id, linksToPrune)
                        FDsRemoved = vertcat(FDsRemoved, i);
                    end
                end
                obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(FDsRemoved) = [];
            end
            
            if(isfieldRecursive(obj.scenario,'SplitRatioSet','splitRatioProfile'))
                SRPsRemoved = [];
                for i = 1:numel(obj.scenario.SplitRatioSet.splitRatioProfile)
                    if ismember(obj.scenario.SplitRatioSet.splitRatioProfile(i).ATTRIBUTE.node_id, nodesRemoved) || ...
                            ( numel(obj.get_node_byID(obj.scenario.SplitRatioSet.splitRatioProfile(i).ATTRIBUTE.node_id).inputs.input) <= 1 && ...
                            numel(obj.get_node_byID(obj.scenario.SplitRatioSet.splitRatioProfile(i).ATTRIBUTE.node_id).outputs.output) <= 1)
                        SRPsRemoved = vertcat(SRPsRemoved, i);
                    end
                end
                obj.scenario.SplitRatioSet.splitRatioProfile(SRPsRemoved) = [];
            end
            
            if(isfieldRecursive(obj.scenario,'SensorSet','sensor'))
                sensorsRemoved = [];
                for i = 1:numel(obj.scenario.SensorSet.sensor)
                    if ~isfield(obj.scenario.SensorSet.sensor(i).ATTRIBUTE,'link_id') || ...
                            ismember(obj.scenario.SensorSet.sensor(i).ATTRIBUTE.link_id, linksToPrune)
                        sensorsRemoved = vertcat(sensorsRemoved, i);
                    end
                end
                obj.scenario.SensorSet.sensor(sensorsRemoved) = [];
            end
            
            % regenerate the link node map
            obj.generate_link_id_begin_end();

        end
        
        % links ...........................................................
        
        function [b] = has_links(obj)
            b = isfieldRecursive(obj.scenario,'NetworkSet','network','LinkList','link');
        end
        
        function [n] = get_num_links(obj)
            if ~obj.has_links
                n = 0;
            else
                n = length(obj.scenario.NetworkSet.network.LinkList.link);
            end
        end
        
        function [linklanes] = get_link_lanes(obj)
            linklanes=[];
            if ~obj.has_links
                return
            end
            numlinks = obj.get_num_links;
            linklanes = nan(1,numlinks);
            for j=1:numlinks
                linklanes(j) = obj.scenario.NetworkSet.network.LinkList.link(j).ATTRIBUTE.lanes;
            end
        end
        
        function [linklengths] = get_link_lengths(obj,desired_units)
            
            if(nargin<2)
                desired_units = 'si';
            end
            linklengths=[];
            if ~obj.has_links
                return
            end
            numlinks = obj.get_num_links;
            linklengths = nan(1,numlinks);
            for j=1:numlinks
                linklengths(j) = obj.scenario.NetworkSet.network.LinkList.link(j).ATTRIBUTE.length;
            end
            
            % unit conversion
            cnv = unit_conversion_factor(obj.get_units(),desired_units);
            linklengths = linklengths*cnv.length;
            
        end
        
        function [linkids] = get_link_ids(obj,ind)
            linkids=[];
            if ~obj.has_links
                return
            end
            numlinks = obj.get_num_links;
            linkids = nan(1,numlinks);
            for j=1:numlinks
                linkids(j) = obj.scenario.NetworkSet.network.LinkList.link(j).ATTRIBUTE.id;
            end
            
            % option to mask with ind
            if(nargin>1)
                linkids=linkids(ind);
            end
        end
        
        function [linktypes] = get_link_types(obj,ind)
            linktypes={};
            if ~obj.has_links
                return
            end
            numlinks = obj.get_num_links;
            linktypes = repmat({''},1,numlinks);
            for j=1:numlinks
                L = obj.scenario.NetworkSet.network.LinkList.link(j);
                if(isfield(L,'link_type'))
                    linktypes{j} = L.link_type.ATTRIBUTE.name;
                end
            end
            
            % option to mask with ind
            if(nargin>1)
                linktypes=linktypes(ind);
            end
        end
        
        function [link_names] = get_link_names(obj)
            link_names={};
            if ~obj.has_links
                return
            end
            numlinks = obj.get_num_links();
            link_names = repmat({''},1,numlinks);
            for j=1:numlinks
                L = obj.scenario.NetworkSet.network.LinkList.link(j);
                if(isfield(L.ATTRIBUTE,'link_name'))
                    link_names{j} = L.ATTRIBUTE.link_name;
                end
            end
        end
        
        function [L] = get_link(obj,link_index)
            L = [];
            if( ~obj.has_network )
                return
            end
            if( any(link_index>length(obj.scenario.NetworkSet.network.LinkList.link)) )
                return
            end
            L = obj.scenario.NetworkSet.network.LinkList.link(link_index);
        end
        
        function [L] = get_link_byID(obj,ids)
            L = [];
            if(~obj.has_links)
                return
            end
            L = obj.scenario.NetworkSet.network.LinkList.link(ismember(obj.get_link_ids,ids));
        end
        
        function [is_src] = is_source_link(obj)
            numlinks = size(obj.link_id_begin_end,1);
            is_src = false(1,numlinks);
            for i=1:numlinks
                begin_node = obj.link_id_begin_end(i,2);
                is_src(i) = obj.is_source_node(begin_node);
            end
        end
        
        function [is_snk] = is_sink_link(obj)
            numlinks = size(obj.link_id_begin_end,1);
            is_snk = false(1,numlinks);
            for i=1:numlinks
                end_node = obj.link_id_begin_end(i,3);
                is_snk(i) = obj.is_sink_node(end_node);
            end
        end
       
        function [is_type] = get_link_is_type(obj,str)
            is_type = strcmpi(obj.get_link_types,str);
        end
        
        function [link_id_begin_end] = get_link_id_begin_end(ptr)
            link_id_begin_end = ptr.link_id_begin_end;
        end
        
        function [fd] = get_link_fd(obj)
            fd2link = obj.get_fd_link_map();
            
            link_ids = obj.get_link_ids;
            num_links = obj.get_num_links;
            fdstruct = struct(  'capacity',nan,...
                'congestion_speed',nan,...
                'free_flow_speed',nan,...
                'jam_density',nan,...
                'transition_density',nan);
            fd = repmat(fdstruct,1,num_links);
            
            for i=1:num_links
                ind = fd2link(:,2) == link_ids(i);
                if(~any(ind))
                    continue
                end
                myfd = obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(ind).fundamentalDiagram(1);
                fd(i).capacity = myfd.ATTRIBUTE.capacity;
                fd(i).congestion_speed = myfd.ATTRIBUTE.congestion_speed;
                fd(i).free_flow_speed = myfd.ATTRIBUTE.free_flow_speed;
                if(isfield(myfd.ATTRIBUTE,'jam_density'))
                    fd(i).jam_density = myfd.ATTRIBUTE.jam_density;
                end
                if(isfield(myfd.ATTRIBUTE,'transition_density'))
                    fd(i).transition_density = myfd.ATTRIBUTE.transition_density;
                end
            end
            
        end
        
        function [pos] = get_link_pos(obj,link_ids)
            pos = nan(length(link_ids),2);
            for i=1:length(link_ids)
                link = obj.get_link_byID(link_ids(i));
                if(isempty(link))
                    continue
                end
                node = obj.get_node_byID(link.begin.ATTRIBUTE.node_id);
                if(isempty(node))
                    continue
                end
                pos(i,1) = node.position.point.ATTRIBUTE.lat;
                pos(i,2) = node.position.point.ATTRIBUTE.lng;
            end
        end
        
        function [ids] = get_link_ids_by_type(obj, type)
            %returns an (n x 1) array listing the ids of the links of the
            %entered type. Case insensitive.
            if (isa(type,'char'))
                link_ids =obj.get_link_ids;
                link_types=obj.get_link_types;
                mask=ismember(lower(link_types), lower(cellstr(type)));
                if (any(mask))
                    ids=reshape(link_ids(1,mask),[],1);
                else
                    error('Link type does not exist in the scenario.');
                end
            else
                error('Type must be a string.');
            end
        end

        % nodes ...........................................................
        
        function [b] = has_nodes(obj)
            b = isfieldRecursive(obj.scenario,'NetworkSet','network','NodeList','node');
        end
        
        function [n] = get_num_nodes(obj)
            if ~obj.has_nodes
                n = 0;
            else
                n = length(obj.scenario.NetworkSet.network.NodeList.node);
            end
        end
        
        function [nodeids] = get_node_ids(obj,ind)
            nodeids=[];
            if ~obj.has_nodes
                return
            end
            numnodes = obj.get_num_nodes;
            nodeids = nan(1,numnodes);
            for j=1:length(obj.scenario.NetworkSet.network.NodeList.node)
                nodeids(j) = obj.scenario.NetworkSet.network.NodeList.node(j).ATTRIBUTE.id;
            end
            
            % option to mask with ind
            if(nargin>1)
                nodeids=nodeids(ind);
            end
        end
        
        function [N] = get_node(obj,node_index)
            N = [];
            if( any(node_index>length(obj.scenario.NetworkSet.network.NodeList.node)) )
                return
            end
            N = obj.scenario.NetworkSet.network.NodeList.node(node_index);
        end
        
        function [is_source] = is_source_node(obj,node_id)
            is_valid_id = ismember(node_id,obj.get_node_ids);
            isnt_end_node = ~any(obj.link_id_begin_end(:,3)==node_id);
            is_source = is_valid_id && isnt_end_node;
        end
        
        function [is_sink] = is_sink_node(obj,node_id)
            is_valid_id = ismember(node_id,obj.get_node_ids);
            isnt_begin_node = ~any(obj.link_id_begin_end(:,2)==node_id);
            is_sink = is_valid_id && isnt_begin_node;
        end
        
        function [isTerminal] = is_terminal_node(obj,node_id)
            isTerminal = obj.is_source_node(node_id) || obj.is_sink_node(node_id);
        end
        
        function [x]=get_node_io(obj,node_ids)
            x = repmat(struct('id',[],'link_in',[],'link_out',[]),1,length(node_ids));
            for i=1:length(x)
                x(i).id = node_ids(i);
                x(i).link_in  = obj.link_id_begin_end(obj.link_id_begin_end(:,3)==node_ids(i),1);
                x(i).link_out = obj.link_id_begin_end(obj.link_id_begin_end(:,2)==node_ids(i),1);
            end
        end
        
        % find terminal nodes based on LinkList
        % return the list of node_id for origin nodes and destination nodes
        % and the list of begin node ids, end node ids, and link types
        function [origin_node_ids,dest_node_ids,begin_node_ids,end_node_ids,link_types]=find_terminal_nodes(obj)
            links = obj.scenario.NetworkSet.network.LinkList.link;
            link_id = zeros(numel(links),1);
            begin_node_ids = zeros(numel(links),1);
            end_node_ids = zeros(numel(links),1);
            link_types = cell(numel(links),1);
            for l = 1 : numel(links)
                link_id(l) = links(l).ATTRIBUTE.id;
                begin_node_ids(l) = links(l).begin.ATTRIBUTE.node_id;
                end_node_ids(l) = links(l).end.ATTRIBUTE.node_id;
                link_types{l} = links(l).link_type.ATTRIBUTE.name;
            end
            
            internal_node_ids = intersect(begin_node_ids,end_node_ids);
            origin_node_ids = setdiff(begin_node_ids,internal_node_ids);
            dest_node_ids = setdiff(end_node_ids,internal_node_ids);
        end
        
        function [N] = get_node_byID(obj,id)
            N = [];
            if ~obj.has_nodes
                return
            end
            N = obj.scenario.NetworkSet.network.NodeList.node(obj.get_node_ids==id);
        end
        
        % sensors .........................................................
        
        function [b] = has_sensors(obj)
            b = isfieldRecursive(obj.scenario,'SensorSet','sensor');
        end
        
        function [n] = get_num_sensors(obj)
            n=0;
            if(~obj.has_sensors)
                return
            end
            n = length(obj.scenario.SensorSet.sensor);
        end
        
        function [sensor_id] = get_sensor_ids(obj)
            sensor_id = [];
            if(obj.has_sensors)
                numsensor = length(obj.scenario.SensorSet.sensor);
                sensor_id = nan(1,numsensor);
                for i=1:numsensor
                    sensor_id(i)=obj.scenario.SensorSet.sensor(i).ATTRIBUTE.id;
                end
            end
        end
        
        function [S] = get_sensor(obj,index)
            S = [];
            if( any(index>obj.get_num_sensors()) )
                return
            end
            S = obj.scenario.SensorSet.sensor(index);
        end
        
        function [vds2id] = get_sensor_vds2id_map(obj)
            if(obj.has_sensors)
                numsensor=length(obj.scenario.SensorSet.sensor);
            else
                numsensor = 0;
            end
            vds2id = nan(numsensor,2);
            for i=1:numsensor
                S=obj.scenario.SensorSet.sensor(i);
                vds2id(i,:) = [S.ATTRIBUTE.sensor_id_original S.ATTRIBUTE.id];
            end
        end
        
        function [sensor2link] = get_sensor_link_map(obj)
            if(obj.has_sensors)
                numsensor=length(obj.scenario.SensorSet.sensor);
            else
                numsensor = 0;
            end
            sensor2link = nan(numsensor,2);
            for i=1:numsensor
                S=obj.scenario.SensorSet.sensor(i);
                sensor2link(i,1) = S.ATTRIBUTE.id;
                if(isfield(S.ATTRIBUTE,'link_id'))
                    sensor2link(i,2) = S.ATTRIBUTE.link_id;
                end
            end
        end
        
        function [sensorlanes] = get_sensor_lanes(obj)
            if(obj.has_sensors)
                numsensor=length(obj.scenario.SensorSet.sensor);
            else
                numsensor = 0;
            end
            sensorlanes = nan(numsensor,2);
            for i=1:numsensor
                S=obj.scenario.SensorSet.sensor(i);
                sensorlanes(i,1) = S.ATTRIBUTE.id;
                if(isfield(S.ATTRIBUTE,'lane_number'))
                    sensorlanes(i,2) = S.ATTRIBUTE.lane_number;
                end
            end
        end
        
        function [sensortype] = get_sensor_link_type(obj)
            if(obj.has_sensors)
                numsensor=length(obj.scenario.SensorSet.sensor);
            else
                numsensor = 0;
            end
            sensortype = cell(numsensor,2);
            for i=1:numsensor
                S=obj.scenario.SensorSet.sensor(i);
                sensortype{i,1} = S.ATTRIBUTE.id;
                if(isfield(S.ATTRIBUTE,'link_type_original'))
                    sensortype{i,2} = S.ATTRIBUTE.link_type_original;
                end
            end
        end
        
        function [latlng] = get_sensor_latlng(obj)
            if(obj.has_sensors)
                numsensor=length(obj.scenario.SensorSet.sensor);
            else
                numsensor = 0;
            end
            latlng = nan(numsensor,2);
            for i=1:numsensor
                S=obj.scenario.SensorSet.sensor(i);
                if(isfieldRecursive(S,'display_position','point','ATTRIBUTE'))
                    latlng(i,1) = S.display_position.point.ATTRIBUTE.lat;
                    latlng(i,2) = S.display_position.point.ATTRIBUTE.lng;
                end
            end
        end
        
        function [vds] = get_vds_for_link_id(obj,link_id)
            
            % load maps
            vds2id = obj.get_sensor_vds2id_map();
            sensor_link_map = obj.get_sensor_link_map();
            
            vds = nan(1,length(link_id));
            for i=1:length(link_id)
                
                % get sensor id
                ind = link_id(i)==sensor_link_map(:,2);
                if(~any(ind))
                    continue
                end
                sensor_id = sensor_link_map(ind,1);
                
                % get corresponding vds
                [~,ind]=ismember(sensor_id,vds2id(:,2));
                if(~any(ind))
                    continue
                end
                myvds = vds2id(ind,1);
                if(length(myvds)>1)
                    warning('Multiple vds stations on this link. Keeping only the first encountered.');
                    myvds = myvds(1);
                end
                vds(i) = myvds;
            end
        end
        
        function success = create_sensors_for_mainline_links(obj, ~, links_to_sense)
            if obj.has_sensors
                error('ScenarioPtr already has a SensorSet. Please remove it before calling create_sensors_for_links.');
            end
            
            numlinks = numel(obj.scenario.NetworkSet.network.LinkList.link);
            linktypes=cell(1,numlinks);
            for i_link = 1:numlinks
                linktypes{i_link} = obj.scenario.NetworkSet.network.LinkList.link(i_link).link_type.ATTRIBUTE.name;
            end
            
            freeway_links = find(strcmpi(linktypes,'Freeway'));
            
            if nargin < 3 || isempty(links_to_sense)
                links_to_sense = freeway_links;
            elseif any(~ismember(links_to_sense,freeway_links))
                error('One or more of the link indices provided are not freeway links');
            end
            obj.scenario.SensorSet.sensor = struct();
            obj.scenario.SensorSet.ATTRIBUTE = struct();
            obj.scenario.SensorSet.ATTRIBUTE.id = 1;
            obj.scenario.SensorSet.ATTRIBUTE.lockedForEdit = 'false';
            obj.scenario.SensorSet.ATTRIBUTE.lockedForHistory = 'false';
            obj.scenario.SensorSet.ATTRIBUTE.name = 'Synthetic-SensorSet';
            obj.scenario.SensorSet.ATTRIBUTE.project_id = 0;
            
            for i_sensor = 1:numel(links_to_sense)
                link_for_this_sensor = obj.scenario.NetworkSet.network.LinkList.link(links_to_sense(i_sensor));
                obj.scenario.SensorSet.sensor(i_sensor).display_position.point = struct();
                obj.scenario.SensorSet.sensor(i_sensor).parameters = struct();
                obj.scenario.SensorSet.sensor(i_sensor).sensor_type = struct();
                obj.scenario.SensorSet.sensor(i_sensor).sensor_type.CONTENT = [];
                obj.scenario.SensorSet.sensor(i_sensor).sensor_type.ATTRIBUTE.description = '';
                obj.scenario.SensorSet.sensor(i_sensor).sensor_type.ATTRIBUTE.id = 1;
                obj.scenario.SensorSet.sensor(i_sensor).sensor_type.ATTRIBUTE.name = 'LOOP';
                obj.scenario.SensorSet.sensor(i_sensor).ATTRIBUTE.id = i_sensor;
                obj.scenario.SensorSet.sensor(i_sensor).ATTRIBUTE.lane_number = link_for_this_sensor.ATTRIBUTE.lanes;
                obj.scenario.SensorSet.sensor(i_sensor).ATTRIBUTE.link_id = link_for_this_sensor.ATTRIBUTE.id;
                obj.scenario.SensorSet.sensor(i_sensor).ATTRIBUTE.sensor_id_original = i_sensor;
            end
            success = true;
        end
        
        function [vds_ordered] = get_ordered_vds(obj)
            
            vds_ordered = [];
            numsensor = obj.get_num_sensors();
            
            if(numsensor==0)
                return
            end
            
            % extract ordered list of freeway indices
            ordered_ind = extract_linear_fwy_indices(obj);
            
            % link ids and types
            linkid = obj.get_link_ids();
            linktype = obj.get_link_types();
            
            % ordered fwy ids
            ordered_id = linkid(ordered_ind);
            ordered_type = linktype(ordered_ind);
            
            % construct map of sensor id, vds, link ref
            sinfo = cell(numsensor,5);
            for i=1:numsensor
                S = obj.get_sensor(i);
                sinfo{i,1} = S.ATTRIBUTE.id;
                sinfo{i,2} = S.ATTRIBUTE.sensor_id_original;
                if(isfield(S.ATTRIBUTE,'link_id') && ~isempty(S.ATTRIBUTE.link_id) )
                    sinfo{i,3} = safe_str2double(S.ATTRIBUTE.link_id);
                    ind = sinfo{i,3}==ordered_id;
                    if(~any(ind))
                        error('link not found')
                    end
                    sinfo{i,4} = find(ind);
                    sinfo{i,5} = ordered_type{ind};
                else
                    sinfo{i,3} = nan;
                    sinfo{i,4} = nan;
                    sinfo{i,5} = '';
                end
            end
            
            % sort
            sinfo=sortrows(sinfo,4);
            
            % return a structure
            vds_ordered = cell2struct(sinfo,{'sensor_id','sensor_vds','link_id','link_ind','link_type'},2);
            
        end
        
        function [is_attached] = get_sensor_isattached(obj)
            is_attached = [];
            if ~obj.has_sensors()
                return
            end
            is_attached = false(1,length(obj.scenario.SensorSet.sensor));
            for i=1:length(is_attached)
                S = obj.scenario.SensorSet.sensor(i);
                if ~isfieldRecursive(S,'ATTRIBUTE','link_id')
                    continue
                end
                is_attached(i) = ~isempty(S.ATTRIBUTE.link_id) & ~isnan(S.ATTRIBUTE.link_id);
            end
        end
        
        % demands .........................................................
        
        function [b] = has_demands(obj)
            b = isfieldRecursive(obj.scenario,'DemandSet','demandProfile');
        end
        
        function [demand2link] = get_demandprofile_link_map(obj)
            if(obj.has_demands)
                numdemand=length(obj.scenario.DemandSet.demandProfile);
            else
                numdemand = 0;
            end
            demand2link = nan(numdemand,2);
            for i=1:numdemand
                dp=obj.scenario.DemandSet.demandProfile(i);
                demand2link(i,1) = dp.ATTRIBUTE.id;
                if(isfield(dp.ATTRIBUTE,'link_id_org'))
                    demand2link(i,2) = dp.ATTRIBUTE.link_id_org;
                end
            end
        end
        
        function [l2d] = get_link_demand_map(obj)
        % returns a hasmap from link id to demands
            l2d = [];
            if(~obj.has_demands)
                return
            end
            numdemand=length(obj.scenario.DemandSet.demandProfile);
            link_id = nan(1,numdemand);
            for i=1:numdemand
                dp=obj.scenario.DemandSet.demandProfile(i);
                if(isfield(dp.ATTRIBUTE,'link_id_org'))
                    link_id(i) = dp.ATTRIBUTE.link_id_org;
                end
                   
                if(isfield(dp.ATTRIBUTE,'dt'))
                    dt = dp.ATTRIBUTE.dt;
                else
                    dt = 300;
                end
                
                if(isfield(dp.ATTRIBUTE,'start_time'))
                    start_time = dp.ATTRIBUTE.start_time;
                else
                    start_time = 0;
                end
                
                if(isfieldRecursive(dp.demand(1),'CONTENT'))
                    nsteps = length(dp.demand(1).CONTENT);
                else 
                    nsteps = 0;                    
                end
                
                if(isfield(dp.ATTRIBUTE,'knob'))
                    knob = dp.ATTRIBUTE.knob;
                else
                    knob = 1;
                end
                
                demand{i} = struct( ...
                    'dp_id',dp.ATTRIBUTE.id,...
                    'time',start_time + (0:dt:dt*(nsteps-1)), ...
                    'knob',knob, ...
                    'value',zeros(1,nsteps) ...
                );
                    
                if(isfield(dp,'demand'))
                    for j=1:length(dp.demand)
                        demand{i}.value = demand{i}.value + dp.demand(j).CONTENT;
                    end
                end
            end            
            l2d = containers.Map(link_id,demand);
        end
        
        function [t] = get_demand_start_time(obj)
            t = [];
            if ~obj.has_demands
                return
            end
            try
                temp = [obj.scenario.DemandSet.demandProfile.ATTRIBUTE];
                t = min([ temp.start_time ]);
            catch
                return
            end
        end
        
        function [t] = get_demand_end_time(obj)
            t = [];
            if ~obj.has_demands
                return
            end
            try
                temp = [obj.scenario.DemandSet.demandProfile.ATTRIBUTE];
                temp2 = [obj.scenario.DemandSet.demandProfile.demand];
                temp2 = { temp2.CONTENT };
                t = max( cellfun( @numel,temp2 ) .* [temp.dt] );
            catch
                return
            end
        end
        
        function [Dem]=get_total_source_flw(obj)
                        
            link_ids = obj.get_link_ids;
            is_source_link = obj.is_source_link;
            source_link_ids = link_ids(is_source_link);
            l2d = obj.get_link_demand_map;
                        
            for i=1:length(source_link_ids)
                source_link_id = source_link_ids(i);
                
                D = l2d(source_link_id);
                
                if(i==1)
                    demand_time = D.time;
                    demand = D.value*D.knob;
                else
                    if(length(demand_time)~=length(D.time) | ~all(demand_time==D.time))
                        error('demand profiles have different lengths')
                    else
                        demand = demand + D.value*D.knob;
                    end
                end
                
            end
            
            Dem.flw = demand;
            Dem.time = demand_time;
            Dem.links = source_link_ids;

        end
        
        function success = set_demandprofile(obj, new_demandSet)
            % new_demandSet: A t x n matrix where n is the number of
            % demands
            % Not robust at all, use with caution
            success = false;
            if ~obj.has_demands
                return
            end
            try
                for index = 1 : numel(obj.scenario.DemandSet.demandProfile)
                    obj.scenario.DemandSet.demandProfile(index).demand.CONTENT = new_demandSet(index, :);
                end
            catch
                return
            end
            success = true;
        end
        
        function success = add_demand_noise_multiplicative(obj, std_dev)
            if ~obj.has_demands
                return
            end
            try
                for index = 1 : numel(obj.scenario.DemandSet.demandProfile)
                    obj.scenario.DemandSet.demandProfile(index).ATTRIBUTE.std_dev_mult = std_dev;
                end
            catch
                return
            end
            success = true;
        end
        
        function [] = set_knob_values(obj, demand_profile_id, new_knob_values)
            %demand_profile_id: a 1 x n matrix where n is the number of knobs to
            %change containing the ids of their respective 'demandProfile'.
            %new_knobs_values : a 1 x n matrix where n is the number of
            %knobs to change containing the new value of each knob,
            %corresponding to the ids given in new_knobs_id.
            demand2link = obj.get_demandprofile_link_map;
            for i=1:size(demand_profile_id)
                ind = demand_profile_id(i)==demand2link(:,1);
                if(~any(ind))
                    error('non-existent demand profile id');
                end
                obj.scenario.DemandSet.demandProfile(ind).ATTRIBUTE.knob = new_knob_values(i);
            end
        end
        
        function [dps]=get_demandprofiles_with_linkIDs(obj,linkids)
            dps = repmat(ScenarioPtr.dp_struct(),1,length(linkids));
            dplink = obj.get_demandprofile_link_map;
            ind = index_into(linkids,dplink(:,2));
            dp_ids = nan(1,length(linkids));
            dp_ids(ind>0) = dplink(ind(ind>0),1);
            demandProfile = obj.get_demandprofiles_with_IDs(dp_ids(ind>0));
            dps(ind>0) = ScenarioPtr.unpack_dps(demandProfile);
        end
        
        function [dps]=get_demandprofiles_with_IDs(obj,ids)
            dps = nan(1,length(ids));
            if(obj.has_dps)
                numdps = length(obj.scenario.DemandSet.demandProfile);
                dpids = nan(1,numdps);
                for i=1:numdps
                    dpids(i)=obj.scenario.DemandSet.demandProfile(i).ATTRIBUTE.id;
                end
                dps = obj.scenario.DemandSet.demandProfile(index_into(ids,dpids));
            end
        end
           
        function [obj] = add_demandprofile(obj,link_id,demand,units)
            
            % convert to correct units
            cnv = unit_conversion_factor(units,obj.get_units);
            demand.val = demand.val * cnv.flw;
            
            % initialize 
            if ~obj.has_demands
                dset = generate_mo('DemandSet');
                dset.ATTRIBUTE.id = 0;
                dset.ATTRIBUTE.project_id = obj.scenario.ATTRIBUTE.project_id;
                obj.scenario.DemandSet = dset;
                obj.scenario.DemandSet.demandProfile = [];
                clear dset
            end
            
            dem_link_map = obj.get_demandprofile_link_map;
            ind = dem_link_map(:,2)==link_id;
            packDem = ScenarioPtr.pack_demandProfile(link_id,demand.val,demand.dt,demand.vehType);
            if(isempty(dem_link_map))
                packDem.ATTRIBUTE.id = 0;
            else
                packDem.ATTRIBUTE.id = max(dem_link_map(:,1))+1;
            end
            if(any(ind))
                obj.scenario.DemandSet.demandProfile(ind) = packDem;
            else
                obj.scenario.DemandSet.demandProfile = ...
                    [obj.scenario.DemandSet.demandProfile packDem];
            end
        end
        
        % split ratios ....................................................
        
        function [b] = has_splits(obj)
            b = isfieldRecursive(obj.scenario,'SplitRatioSet','splitRatioProfile');
        end
        
        function [split2node] = get_split_node_map(obj)
            if(obj.has_splits)
                numsplits=length(obj.scenario.SplitRatioSet.splitRatioProfile);
            else
                numsplits = 0;
            end
            split2node = nan(numsplits,2);
            for i=1:numsplits
                srp=obj.scenario.SplitRatioSet.splitRatioProfile(i);
                split2node(i,1) = srp.ATTRIBUTE.id;
                if(isfield(srp.ATTRIBUTE,'node_id'))
                    split2node(i,2) = srp.ATTRIBUTE.node_id;
                end
            end
        end
        
        function success = set_splitRatioSet(obj, new_SRset)
            % new_SRset: A n-length cell array whre n is the number of
            % split ratio profiles. Each entry in the cell array should be
            % an m-length cell array, where m is the number of split ratios
            % in the split ratio profile. Each entry in this level of cell
            % array should be a t-length vector of doubles, where t is is
            % the number of timesteps to write in the new split ratio.
            success = false;
            if ~obj.has_splits
                return
            end
            if ~iscell(new_SRset) || numel(new_SRset) ~= numel(obj.scenario.SplitRatioSet.splitRatioProfile)
                error('The new split ratio set must be a cell array of the same length as the number of SR profiles in the scenario''s split ratio set.');
            end
            for i_SRP = 1:numel(new_SRset)
                if numel(new_SRset{i_SRP}) ~= numel(obj.scenario.SplitRatioSet.splitRatioProfile(i_SRP).splitratio)
                    warning('Split ratio profile number %g has a different number of split ratios than the new profile. Skipping this profile.', i_SRP);
                    continue
                end
                for i_SR = 1:numel(new_SRset{i_SRP})
                    obj.scenario.SplitRatioSet.splititRatioProfile(i_SRP).splitratio(i_SR).CONTENT = ...
                        new_SRset{i_SRP}{i_SR};
                end
            end
            success = true;
        end
        
        function [srp,ind] = get_split_ratios_for_node_id(obj,node_id)
            ind=0;
            srp = [];
            if ~obj.has_splits
                return;
            end
            for i=1:length(obj.scenario.SplitRatioSet.splitRatioProfile)
                if(node_id==obj.scenario.SplitRatioSet.splitRatioProfile(i).ATTRIBUTE.node_id)
                    srp = obj.scenario.SplitRatioSet.splitRatioProfile(i);
                    ind = i;
                    return
                end
            end
            return
        end
        
        function [obj] = add_splitprofile(obj,node_id,dt,splits)
            
            
            % initialize 
            if ~obj.has_splits
                srset = generate_mo('SplitRatioSet');
                srset.ATTRIBUTE.id = 0;
                srset.ATTRIBUTE.project_id = obj.scenario.ATTRIBUTE.project_id;
                obj.scenario.SplitRatioSet = srset;
                obj.scenario.SplitRatioSet.splitRatioProfile = [];
                clear srset
            end
                        
            split_node_map = obj.get_split_node_map;
            ind = split_node_map(:,2)==node_id;
            packSplit = ScenarioPtr.pack_splitProfile(node_id,dt,splits);
            if(isempty(split_node_map))
                packSplit.ATTRIBUTE.id = 0;
            else
                packSplit.ATTRIBUTE.id = max(split_node_map(:,1))+1;
            end
            if(any(ind))
                obj.scenario.SplitRatioSet.splitRatioProfile(ind) = packSplit;
            else
                obj.scenario.SplitRatioSet.splitRatioProfile = ...
                    [obj.scenario.SplitRatioSet.splitRatioProfile packSplit];
            end
            
        end
       
        % fundamental diagrams ............................................
        
        function [b] = has_fds(obj)
            b = isfieldRecursive(obj.scenario,'FundamentalDiagramSet','fundamentalDiagramProfile');
        end
        
        function [fd2link] = get_fd_link_map(obj)
            if obj.has_fds
                numfd=length(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile);
            else
                numfd = 0;
            end
            fd2link = nan(numfd,2);
            for i=1:numfd
                fdp=obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i);
                fd2link(i,1) = fdp.ATTRIBUTE.id;
                if(isfield(fdp.ATTRIBUTE,'link_id'))
                    fd2link(i,2) = fdp.ATTRIBUTE.link_id;
                end
            end
        end
        
        function success = add_fd_capacity_noise_multiplicative(obj, std_dev)
            if ~obj.has_fds
                return
            end
            try
                for index = 1 : numel(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile)
                    for index2 = 1:numel(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(index).fundamentalDiagram)
                        obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(index).fundamentalDiagram(index2).ATTRIBUTE.std_dev_capacity = ...
                            obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(index).fundamentalDiagram(index2).ATTRIBUTE.capacity * std_dev;
                    end
                end
            catch
                return
            end
            success = true;
        end
        
        function [fds]=get_fds_with_linkIDs(obj,linkids)
            fds = repmat(ScenarioPtr.fd_struct(),1,length(linkids));
            fdlink = obj.get_fd_link_map;
            ind = index_into(linkids,fdlink(:,2));
            fd_ids = nan(1,length(linkids));
            fd_ids(ind>0) = fdlink(ind(ind>0),1);
            FundDiag = obj.get_fds_with_IDs(fd_ids(ind>0));
            fds(ind>0) = ScenarioPtr.unpack_fds(FundDiag);
        end
        
        function [fds]=get_fds_with_IDs(obj,ids)
            fds = nan(1,length(ids));
            if(obj.has_fds)
                numfds = length(obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile);
                fdids = nan(1,numfds);
                for i=1:numfds
                    fdids(i)=obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(i).ATTRIBUTE.id;
                end
                fds = obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(index_into(ids,fdids));
            end
        end
           
        function [obj] = add_fd(obj,link_id,fd,units)
            
            % convert to correct units
            cnv = unit_conversion_factor(units,obj.get_units);
            fd.capacity = fd.capacity * cnv.flw;
            fd.free_flow_speed = fd.free_flow_speed * cnv.spd;
            fd.jam_density = fd.jam_density * cnv.dty;
            ncrit = fd.capacity / fd.free_flow_speed;
            fd.congestion_speed = fd.capacity/(fd.jam_density-ncrit);
            
            % initialize fds
            if ~obj.has_fds
                fdset = generate_mo('FundamentalDiagramSet');
                fdset.ATTRIBUTE.id = 0;
                fdset.ATTRIBUTE.project_id = obj.scenario.ATTRIBUTE.project_id;
                obj.scenario.FundamentalDiagramSet = fdset;
                obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile = [];
                clear fdset
            end
            
            fd_link_map = obj.get_fd_link_map;
            ind = fd_link_map(:,2)==link_id;
            packFd = ScenarioPtr.pack_fds(link_id,fd);
            if(isempty(fd_link_map))
                packFd.ATTRIBUTE.id = 0;
            else
                packFd.ATTRIBUTE.id = max(fd_link_map(:,1))+1;
            end
            if(any(ind))
                obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile(ind) = packFd;
            else
                obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile = ...
                    [obj.scenario.FundamentalDiagramSet.fundamentalDiagramProfile packFd];
            end
             
        end
        
        % intial densities ................................................
        
        function [b] = has_initial_densities(obj)
            b = isfieldRecursive(obj.scenario,'InitialDensitySet','density');
        end
        
        % capacity profiles ...............................................
        
        function [b] = has_capacity_profiles(obj)
            b = isfieldRecursive(obj.scenario,'DownstreamBoundaryCapacitySet','downstreamBoundaryCapacityProfile');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% more complicated getters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [T] = get_sensor_table(obj)
                    
            sensor_ids = obj.get_sensor_ids;
            is_attached = obj.get_sensor_isattached;
            latlng = obj.get_sensor_latlng;
            sensor_link = obj.get_sensor_link_map;
            sensor_link_type = obj.get_sensor_link_type;
            vds2id = obj.get_sensor_vds2id_map;
            sensor_lanes = obj.get_sensor_lanes;
            link_ids = obj.get_link_ids;
            link_types = obj.get_link_types;

            for i=1:length(sensor_ids)
                C{i,1} = sensor_ids(i);
                C{i,2} = double(is_attached(i));
                C{i,3} = latlng(i,1);
                C{i,4} = latlng(i,2);
                C{i,5} = sensor_link(i,2);
                C{i,6} = sensor_link_type(i,2);
                C{i,7} = vds2id(vds2id(:,2)==sensor_ids(i),1);
                C{i,8} = sensor_lanes(i,2);
                if is_attached(i)
                    C{i,9} = link_types(link_ids==sensor_link(i,2));
                else
                    C{i,9} = '';
                end
            end
            
            if ~isempty(sensor_ids)
                T = cell2table(C,'VariableNames',{ ...
                    'id' ...
                    'is_attached' ...
                    'lat' ...
                    'lng' ...
                    'link_id' ...
                    'sensor_link_type' ...
                    'vds' ...
                    'lanes' ...
                    'link_type' });
            else
                T = [];
            end
        end
        
        function [T]=get_link_table(obj)
            
            link_types = obj.get_link_types;
            link_lanes = obj.get_link_lanes;
            link_lengths = obj.get_link_lengths;
            sensor_link = obj.get_sensor_link_map;
            sensor_id_vds = obj.get_sensor_vds2id_map;
            fds = obj.get_fds_with_linkIDs(obj.get_link_ids);
            for i=1:size(obj.link_id_begin_end,1)
                link_id = obj.link_id_begin_end(i,1);
                begin_node = obj.get_node_byID(obj.link_id_begin_end(i,2));
                C{i,1} = link_id;
                C{i,2} = link_types{i};
                C{i,3} = link_lanes(i);
                C{i,4} = link_lengths(i);
                C{i,5} = obj.link_id_begin_end(i,2);
                C{i,6} = obj.link_id_begin_end(i,3);
                C{i,7} = begin_node.position.point.ATTRIBUTE.lat;
                C{i,8} = begin_node.position.point.ATTRIBUTE.lng;
                sensor_ind = link_id==sensor_link(:,2);
                if(any(sensor_ind))
                    sensor_ids = sensor_link(sensor_ind,1);
                    sensor_vds = sensor_id_vds(index_into(sensor_ids,sensor_id_vds(:,2)),1);
                else
                    sensor_ids = [];
                    sensor_vds = [];
                end
                C{i,9} = sensor_ids;
                C{i,10} = sensor_vds;
                C{i,11} = fds(i).free_flow_speed;
                C{i,12} = fds(i).capacity;
                C{i,13} = fds(i).congestion_speed;
            end
            
            T = cell2table(C,'VariableNames',{ ...
                    'id' ...
                    'type' ...
                    'lanes' ...
                    'length' ...
                    'begin_node' ...
                    'end_node' ...
                    'lat' ...
                    'lng' ...
                    'sensors' ...
                    'vds',...
                    'free_flow_speed',...
                    'capacity',...
                    'congestion_speed'});
            
        end
        
        function [data]=get_5min_date_for_day(obj,day,data_folder)
            
            vds2id = obj.get_sensor_vds2id_map;
            vds = vds2id(:,1);
            
            matfile = fullfile(data_folder,['pems5min_' datestr(day,'yyyy') '_' datestr(day,'mm') '_' datestr(day,'dd')]);
            
            if(~exist([matfile '.mat'],'file'))
                error(['Data file not found (' matfile ')'])
            end
            
            load(matfile)
            
            % construct time vector
            numTime = 288;
            dayvec = datevec(day);
            dayvec = repmat(dayvec,numTime,1);
            
            dayvec(:,4) = reshape(repmat((0:23),12,1),288,1);
            dayvec(:,5) = repmat((0:5:55)',24,1);
            
            time_vec = datenum(dayvec);
            time_vec = round(time_vec*1440)/1440;    % round to the nearest minute
            
            [hasvds,indexvds]=ismember(vds,pems.vds);
            flw = nan(numTime,length(vds));
            occ = nan(numTime,length(vds));
            spd = nan(numTime,length(vds));
            for v=1:length(vds)
                if(~hasvds(v))
                    continue
                end
                pemsflw = pems.data(indexvds(v)).flw;
                pemsocc = pems.data(indexvds(v)).occ;
                pemsspd = pems.data(indexvds(v)).spd;
                pemstime = pems.data(indexvds(v)).time;
                pemstime = round(pemstime*1440)/1440;    % round to the nearest minute
                [hastime,indextime]=ismember(pemstime,time_vec);
                if(any(hastime))
                    flw(indextime(hastime),v) = pemsflw(hastime);
                    occ(indextime(hastime),v) = pemsocc(hastime);
                    spd(indextime(hastime),v) = pemsspd(hastime);
                end
            end
            
            data.flw = flw;
            data.occ = occ;
            data.spd = spd;
            data.timeAbs = time_vec;
            data.timeSec = (0:300:86100)';
        end
        
        function [health] = get_loop_health_for_day(obj,day,data_folder)
            matfile = fullfile(data_folder,'loop_health');
            if(~exist([matfile '.mat'],'file'))
                error(['Data file not found (' matfile ')'])
            end
            load(matfile)
            vds2id = obj.get_sensor_vds2id_map;
            vds = vds2id(:,1);
            health = nan(1,length(vds));
            day_ind = loop_health.days==day;
            for i=1:length(vds)
                vds_ind=vds(i)==loop_health.vds;
                health(i) = loop_health.percent_observed(vds_ind,day_ind)==100;
            end
        end
        
        %% Freeway specific
        
        function [T] = get_linear_fwy_sensors(obj)
            
            ordered_ind = obj.extract_linear_fwy_indices';
            link_ids = obj.get_link_ids';
            link_ids = link_ids(ordered_ind);
            
            link_types = obj.get_link_types';
            link_types = link_types(ordered_ind);
            
            sensor_link = obj.get_sensor_link_map;
            
            sensor_link(isnan(sensor_link(:,2)),:) = [];
            
            sensor_ids = nan(length(link_ids),1);
            [ismem,ind] = ismember(link_ids,sensor_link(:,2));
            sensor_ids(ismem) = sensor_link(ind(ismem),1);
            
            vds2id = obj.get_sensor_vds2id_map;
            
            [ismem,ind] = ismember(sensor_ids,vds2id(:,2));
            sensor_vds = nan(length(link_ids),1);
            sensor_vds(ismem) = vds2id(ind(ismem),1);
            
            T = table(link_ids,sensor_ids,sensor_vds,link_types);
            
        end
             
        function [ordered_ind] = extract_linear_fwy_indices(obj)
            % traverses a freeway network, where every node has at most
            % two entry links (one freeway/one onramp) and two exit nodes
            % (one freeway, one offramp). Return an array of ordered indices to
            % LinkList. For junctions with an onramp and/or offramp, the order is:
            % [entering freeway , offramp , onramp , exiting freeway]/
            
            %% gather
            
            % from nodes, allnode_id
            node_ids = obj.get_node_ids();
            
            % from links
            numlinks = obj.get_num_links();
            link_ids = obj.get_link_ids();
            link_types = obj.get_link_types();
            is_fwy_link = strcmpi(link_types,'freeway');
            is_source_link = obj.is_source_link;
            link_begin_ind = nan(1,numlinks);
            link_end_ind = nan(1,numlinks);
            next_is_fwy_link = false(1,numlinks);
            end_node_is_simple = false(1,numlinks);
            for i=1:numlinks
                L = obj.get_link(i);
                link_begin_ind(i)=find(L.begin.ATTRIBUTE.node_id==node_ids);
                link_end_ind(i)=find(L.end.ATTRIBUTE.node_id==node_ids);
                node_io = obj.get_node_io(L.end.ATTRIBUTE.node_id);
                if( length(node_io.link_out)==1 )
                    nextLink = obj.get_link(node_io.link_out==link_ids);
                    next_is_fwy_link(i) = strcmpi(nextLink.link_type.ATTRIBUTE.name,'freeway');
                end
                end_node_is_simple(i) = length(node_io.link_in)<=1 && length(node_io.link_out)<=1;
            end
            clear L N
            
            %% find first freeway link
            % fwy_ind = find(is_source_link & (end_node_is_simple & next_is_fwy_link));
            fwy_ind = find( (is_fwy_link & is_source_link) | ...
                (is_source_link  & end_node_is_simple & next_is_fwy_link) );
            if(isempty(fwy_ind))
                error('first fwy link not found')
            end
            if(length(fwy_ind)>1)
                error('more than one freeway source')
            end
            
            %% traverse
            ordered_ind = fwy_ind;
            remaining_link_id = link_ids;
            remaining_link_id(fwy_ind) = nan;
            % remaining_link_id(~is_fwy_link & ~is_or_link & ~is_fr_link) = nan;
            while(true)
                
                end_node = obj.get_node(link_end_ind(fwy_ind));
                
                % stop if we've reached a terminal node
                if(is_end_node(end_node))
                    break
                end
                
                % no more freeway links to include
                if(all(isnan(remaining_link_id)))
                    error('loop found')
                end
                
                % get onramp, offramp, and fwy
                fwy_ind = get_incident(end_node.outputs.output,remaining_link_id,is_fwy_link);
                out_not_fwy = get_incident(end_node.outputs.output,remaining_link_id,~is_fwy_link);
                in_not_fwy = get_incident(end_node.inputs.input,remaining_link_id,~is_fwy_link);
                
                %     fr_ind = get_incident(end_node.outputs.output,remaining_link_id,is_fr_link);
                %     fwy_ind = get_incident(end_node.outputs.output,remaining_link_id,is_fwy_link);
                %     or_ind = get_incident(end_node.inputs.input,remaining_link_id,is_or_link);
                
                if(length(fwy_ind)~=1)
                    error('node has 0 or >1 freeway outputs')
                end
                
                % record
                ordered_ind = [ordered_ind out_not_fwy in_not_fwy fwy_ind];
                
                remaining_link_id(fwy_ind) = nan;
                remaining_link_id(out_not_fwy) = nan;
                remaining_link_id(in_not_fwy) = nan;
                
            end
            
            function [x] = get_incident(Nlinks,alllink_id,alllink_istype)
                x = [];
                for ii=1:length(Nlinks)
                    ind = Nlinks(ii).ATTRIBUTE.link_id==alllink_id;
                    if(~any(ind))
                        continue;
                    end
                    if(alllink_istype(ind))
                        x = [x find(ind)];
                    end
                end
            end
            
        end
           
        function [F] = get_freeway_structure(obj)
            
            F.linear_fwy_sensors = obj.get_linear_fwy_sensors;
            F.linear_fwy_ind = obj.extract_linear_fwy_indices;
            F.linear_fwy_link_ids = obj.get_link_ids(F.linear_fwy_ind);
            
            % get mainline links
            ml_ind = find(obj.get_link_is_type('Freeway'));
            F.linear_ml_link_ind = F.linear_fwy_ind(ismember(F.linear_fwy_ind,ml_ind));
            F.linear_ml_link_ids = obj.get_link_ids(F.linear_ml_link_ind);
            
            % get ml nodes
            x = index_into(F.linear_ml_link_ids,obj.link_id_begin_end(:,1));
            begin_nodes = obj.link_id_begin_end(x,2)';
            end_nodes = obj.link_id_begin_end(x,3)';
            
            if(~all(begin_nodes(2:end)==end_nodes(1:end-1)))
                error('the nodes dont line up');
            end
            
            F.linear_ml_node_ids = [begin_nodes(1) end_nodes];
            
            % make tall matrices
            F.linear_fwy_ind = F.linear_fwy_ind';
            F.linear_fwy_link_ids = F.linear_fwy_link_ids';
            F.linear_ml_link_ind = F.linear_ml_link_ind';
            F.linear_ml_link_ids = F.linear_ml_link_ids';
            F.linear_ml_node_ids = F.linear_ml_node_ids';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% plots/reporting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [] = report(obj,pptfile)
            report_scenario(obj,pptfile)
        end
        
        function [T]=report_CFL_table(obj)

            cnv = unit_conversion_factor(obj.get_units,'SI');
            link_ids = obj.get_link_ids;
            fds = obj.get_fds_with_linkIDs(link_ids);

            link_lengths = obj.get_link_lengths;
            length_in_meters = link_lengths*cnv.length;
               
            maxdt_per_link = nan(1,length(fds));
            for i=1:length(fds)
               if(isfield(fds(i),'free_flow_speed') & ~isempty(fds(i).free_flow_speed))
                   vf = fds(i).free_flow_speed * cnv.spd;
               else
                   vf = 60 * 0.4470;
               end
               
               if(isfield(fds(i),'congestion_speed') & ~isempty(fds(i).congestion_speed))
                   w = fds(i).congestion_speed * cnv.spd;
               else
                   w = 20 * 0.4470;
               end
            
               maxdt_per_link(i) = length_in_meters(i) / max([vf w]);
               
            end

            T=table(link_ids',obj.get_link_types',maxdt_per_link',...
                'VariableNames',{'link_id' 'link_type' 'max_dt'});
            T = sortrows(T,'max_dt');
        end
                
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% private
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = private)
        
        function X=remove_if_nan(~,X,f)
            for i=1:length(f)
                if(isfield(X,f{i}) && isnan(X.(f{i})) )
                    X = safe_rmfield(X,f(i));
                end
            end
        end
        
        function X=replace_nan_with(obj,X,f,val)
            if(~isfield(X,f{1}))
                return
            end
            if(length(f)==1)
                if(isnan(X.(f{1})))
                    X.(f{1})=val;
                end
            else
                for i=1:length(X.(f{1}))
                    X.(f{1})(i)=obj.replace_nan_with(X.(f{1})(i),f(2:end),val);
                end
            end
        end
        
        function X=replace_empty_with(obj,X,f,val)
            if(~isfield(X,f{1}))
                return
            end
            if(length(f)==1)
                if(isempty(X.(f{1})))
                    X.(f{1})=val;
                end
            else
                for i=1:length(X.(f{1}))
                    X.(f{1})(i)=obj.replace_empty_with(X.(f{1})(i),f(2:end),val);
                end
            end
        end
        
        function in=default_to(~,in,field,val)
            if(~isfield(in,field) || isempty(in.(field)))
                in.(field) = val;
            end
        end
        
        function [X] = safe_conv(obj,X,c,varargin)
            if(~isfield(X,varargin{1}))
                return
            end
            if(length(varargin)>1)
                X.(varargin{1}) = obj.safe_conv(X.(varargin{1}),c,varargin{2:end});
            else
                X.(varargin{1}) = X.(varargin{1})*c;
            end
        end
        
        function [b] = has_dps(obj)
            b = isfieldRecursive(obj.scenario,'DemandSet','demandProfile');
        end

       
    end

    methods (Static)

        function [b] = arrays_are_same(A,B)
            if(isempty(A) && isempty(B))
                b = true;
                return;
            end
            A = reshape(A,1,length(A));
            B = reshape(B,1,length(B));
            b = isequal(sort(A),sort(B));
            return;
        end
        
        function [x]=fd_struct()
            x=struct('free_flow_speed',[],...
                                'capacity',[],...
                                'critical_density',[],...
                                'congestion_speed',[],...
                                'jam_density',[],...
                                'capacity_drop',[]);
        end
                
        function [x]=dp_struct()
            x=struct('dt',[],...
                                'id',[],...
                                'knob',[],...
                                'link_id_org',[],...
                                'demand',[]);
        end
        
        function [fds]=unpack_fds(FDs)
            fds = repmat(ScenarioPtr.fd_struct,1,length(FDs));
            for i=1:length(FDs)
                A = FDs(i).fundamentalDiagram.ATTRIBUTE;
                fds(i).free_flow_speed = A.free_flow_speed;
                fds(i).congestion_speed = A.congestion_speed;
                fds(i).capacity = A.capacity;
                if(isfield(A,'critical_density'))
                    fds(i).critical_density = A.critical_density;
                end
                if(isfield(A,'jam_density'))
                    fds(i).jam_density = A.jam_density;
                end
                if(isfield(A,'capacity_drop'))
                    fds(i).capacity_drop = A.capacity_drop;
                end           
            end
        end
        
        function [FDs]=pack_fds(link_id,fds)
            x = generate_mo('fundamentalDiagramProfile');
            x.fundamentalDiagram = generate_mo('fundamentalDiagram');
            x.ATTRIBUTE.id = 0;
            x.ATTRIBUTE.link_id = link_id;
            FDs = repmat(x,1,length(fds));
            for i=1:length(fds)
               FDs(i).fundamentalDiagram.ATTRIBUTE.id = 0;
               FDs(i).fundamentalDiagram.ATTRIBUTE.capacity = fds(i).capacity;
               FDs(i).fundamentalDiagram.ATTRIBUTE.free_flow_speed = fds(i).free_flow_speed;
               FDs(i).fundamentalDiagram.ATTRIBUTE.congestion_speed =fds(i).congestion_speed;
            end
        end        
        
        function [x]=pack_demandProfile(link_id,dem,dt,vtype)
            x = generate_mo('demandProfile');
            x.ATTRIBUTE.id = 0;
            x.ATTRIBUTE.start_time = 0;
            x.ATTRIBUTE.dt = dt;
            x.ATTRIBUTE.knob = 1;
            x.ATTRIBUTE.link_id_org = link_id;
            x.demand = generate_mo('demand',true);
            x.demand.ATTRIBUTE.vehicle_type_id = vtype;
            x.demand.CONTENT = writecommaformat(dem,'%f');
        end        
        
        function [x]=pack_splitProfile(node_id,dt,splits)
            x = generate_mo('splitRatioProfile');
            x.ATTRIBUTE.id = 0;
            x.ATTRIBUTE.start_time = 0;
            x.ATTRIBUTE.dt = dt;
            x.ATTRIBUTE.node_id = node_id;
            x.splitratio = repmat(generate_mo('splitratio'),1,length(splits));
            for i=1:length(splits)
                x.splitratio(i).ATTRIBUTE.link_in = splits(i).in_link;
                x.splitratio(i).ATTRIBUTE.link_out = splits(i).out_link;
                x.splitratio(i).ATTRIBUTE.vehicle_type_id = splits(i).vehType;
                x.splitratio(i).CONTENT = writecommaformat(splits(i).val,'%f');
            end
        end        
        
        
        function [dps]=unpack_dps(DPs)
            dps = repmat(ScenarioPtr.dp_struct,1,length(DPs));
            for i=1:length(DPs)
                A = DPs(i).ATTRIBUTE;
                dps(i).dt = A.dt;
                dps(i).id = A.id;
                dps(i).knob = A.knob;
                dps(i).link_id_org = A.link_id_org;
                dps(i).demand = DPs(i).demand.CONTENT;
            end
        end    
                
    end
    
end
