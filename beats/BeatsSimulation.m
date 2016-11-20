classdef BeatsSimulation < handle
    
    properties
        
        scenario_ptr        % ScenarioPtr
        configfile          % string name of the config file
        
        % used by beats persistent mode
        beats    
        initialized
        temp_output
        
        % state of the object ..............
        config_loaded       % boolean
        simulation_done     % boolean
        sim_output_loaded   % boolean
        
        % simulation information
        sim_dt              % in seconds
        out_dt              % in seconds
        
        % simulation output ................
        time
        density_veh
        outflow_veh
        inflow_veh
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Construction, loading scenarios                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        % Construction ......................................
        function [obj] = BeatsSimulation()
            obj.reset_object();
        end
                    
        function [out] = clone(obj)
            out = BeatsSimulation;
            out.scenario_ptr = obj.scenario_ptr.clone;
            out.configfile = obj.configfile;  
            
            % in lieu of Object.clone() in matlab, we need to create
            % a new beats object
            if ~isempty(obj.beats)
                out.beats = obj.create_beats_object( struct('SIM_DT',obj.sim_dt) );
                out.initialized = false;
                out.temp_output = '';
            end
            
            out.config_loaded = obj.config_loaded;   
            out.simulation_done = obj.simulation_done;   
            out.sim_output_loaded = obj.sim_output_loaded; 
            
            out.sim_dt = obj.sim_dt;
            out.out_dt = obj.out_dt;
            
            out.time = obj.time;   
            out.density_veh = obj.density_veh;   
            out.outflow_veh = obj.outflow_veh;   
            out.inflow_veh = obj.inflow_veh;  
        end
            
        % load scenario from xml ............................
        function [obj]=load_scenario(obj,cfgfile)
            obj.reset_object();
            [~,success] = obj.scenario_ptr.load(cfgfile);
            if(success)
                obj.config_loaded = true;
                obj.configfile = cfgfile;
            else
                error('load failed.')
            end
        end
        
        % loadfwy scenario from Excel ............................
        function [obj]=load_fwy_from_excel(obj,excelfile)
            obj.reset_object();
            obj.scenario_ptr.load_fwy_from_excel(excelfile);
        end
        
        % load scenario from properties ............................
        function [obj]=load_properties(obj,prop_file)
            
            % load properties
            props = BeatsSimulation.read_properties_file(prop_file);
            
            % load the configuration file
            obj=obj.load_scenario(props(strcmpi({props.name},'scenario')).value);
            
            % check whether the output is there. if so, load it
            out_prefix = props(strcmpi({props.name},'OUTPUT_PREFIX')).value;
            have_run = ~isempty(dir([out_prefix '_*.txt']));
            if(have_run)
                obj = obj.load_simulation_output(out_prefix);
            end
        end
        
        function [obj]=set_scenario(obj,scenario_ptr)
            obj.reset_simulation();
            obj.scenario_ptr = scenario_ptr;
            obj.configfile = '';
            obj.config_loaded = true;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  running BeATS                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        % new runner, creates a properties file from the param struct 
        % and runs beats with -s flag
        function [obj]=run_beats(obj,param)

            warning('on','all')
            
            % deal with inputs
            if(nargin<2)
                param = struct;
            end
                        
            % create temporary scenario
            temp_scenario = [tempname '.xml'];
            param.SCENARIO = temp_scenario;
            obj.scenario_ptr.save(temp_scenario);
            
            % build and write prop file
            [param,prop_file] = obj.create_properties_file(param);

            % choose beats jar
            if(isfield(param,'JAR'))
                beats_jar = param.JAR;
            else
                root = fileparts(mfilename('fullpath'));
                beats_jar = fullfile(root,'beats-0.1-SNAPSHOT-jar-with-dependencies.jar');
            end
            
            % run beats
            status = system(['java -jar "' beats_jar '" -s ' prop_file]);
            
            % delete prop file
            system(sprintf('del %s',prop_file));            
            
            % normal termination
            if(status==0)
                obj.simulation_done = true;
                obj.load_simulation_output(obj.temp_output);
                if ~isempty(obj.temp_output)
                    system(sprintf('del "%s_*.txt"',obj.temp_output));
                end
            end

            % abnormal termination
            if(status~=0)
                if ~isempty(obj.temp_output)
                    system(sprintf('del "%s_*.txt"',obj.temp_output));
                end
                error('failed to complete the simulation')
            end

        end
        
        function [obj]=reset_simulation(obj)
            obj.reset_simdata();
            if ~isempty(obj.beats)
                obj.beats.reset();
            end
        end
        
        function [obj]=load_simulation_output(obj,outprefix)
                        
            if(obj.sim_output_loaded)
                warning('overwriting existing simulation data.')
                obj.reset_simdata();
            end
            
            disp('Loading simulation result')
            
            % ... vehicle type names
            vtypes = obj.scenario_ptr.get_vehicle_types();
            
            % load data
            obj.time = load(sprintf('%s_%s_0.txt',outprefix,'time'));
            numlinks = length(obj.scenario_ptr.scenario.NetworkSet.network(1).LinkList.link);
            numvt = length(vtypes);
            
            obj.density_veh = cell(1,numvt);
            obj.outflow_veh = cell(1,numvt);
            obj.inflow_veh = cell(1,numvt);
            
            [obj.density_veh{1:numvt}] = deal(nan(length(obj.time),numlinks));
            [obj.outflow_veh{1:numvt}] = deal(nan(length(obj.time),numlinks));
            [obj.inflow_veh{1:numvt}] = deal(nan(length(obj.time),numlinks));
            for v=1:numvt
                obj.density_veh{v} = load(sprintf('%s_%s_%s_0.txt',outprefix,'density',vtypes(v).name));
                obj.outflow_veh{v} = load(sprintf('%s_%s_%s_0.txt',outprefix,'outflow',vtypes(v).name));
                obj.inflow_veh{v} = load(sprintf('%s_%s_%s_0.txt',outprefix,'inflow',vtypes(v).name));
            end
            
            obj.sim_output_loaded = true;
            
            % fix for working with cc-sim
            if(length(obj.time)<size(obj.density_veh{1},1))
                n = size(obj.density_veh{1},1)-length(obj.time);
                obj.time(end:end+n)=obj.time(end);
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  running BeATS with persistent java object                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
       
        function obj = create_beats_object(obj, param)
            %obj.scenario_ptr must be loaded first.
            
            warning('on','all')
            
            obj.sim_dt = param.SIM_DT;
            
            % create temporary scenario
            temp_scenario = [tempname '.xml'];
            param.SCENARIO = temp_scenario;
            obj.scenario_ptr.save(temp_scenario);

            % load the configuration into ScenarioPtr
            obj.beats = edu.berkeley.path.beats.Jaxb.create_scenario_from_xml(param.SCENARIO);
            obj.initialized = false;
            
            % print beats version
%             fprintf('BeATS version: %s\n',char(edu.berkeley.path.beats.Version.getGitHash));
                        
            % delete temporary scenario
            system(sprintf('del "%s"',temp_scenario));
        end
        
        function [obj,param]=intialize_beats_persistent( obj , param )
            
            % build and write propparam
            [param,prop_file] = obj.create_properties_file(param);
            
            % intialize beats
            BeatsProperties=edu.berkeley.path.beats.simulator.BeatsProperties(prop_file);            
            try
                obj.beats.initialize_with_properties(BeatsProperties);
            catch JavaException
                error(['Error in initializing the BeATS scenario: ', JavaException.message]);
            end
            
            % delete properties file
            system(sprintf('del %s',prop_file));

            obj.initialized = true;
        end
        
        function [obj]=run_beats_persistent(obj,duration,param)  

            % initialize
            if ~obj.initialized
                obj.intialize_beats_persistent( param );
            end
            
            % run beats
            obj.beats.set.end_time(duration-obj.sim_dt);   % NOTE: this should be "duration", but for a bug in BeATS
            obj.beats.run();
            obj.simulation_done = true;
            
            % load the result
            obj.load_simulation_output(obj.temp_output);
            
            % delete output files
            if ~isempty(obj.temp_output)
                system(sprintf('del "%s_*.txt"',obj.temp_output));
            end
        end
        
%         function [obj]=advance_beats_persistent( obj , advance_sec )
%             if ~obj.initialized
%                 warning('initialize first')
%                 return
%             end 
%             obj.beats.advanceNSeconds(advance_sec);
%             obj.simulation_done = true;
%         end
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  query the simulation                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
       
        function [X]=get_output_for_link_id(obj,link_id)
            
            if(~obj.sim_output_loaded)
                error('run simulation first')
            end
            
            numtime = length(obj.time);
            X = struct('link_id', link_id,...
                       'time_in_sec',obj.time, ...
                       'flw_in_veh',nan(numtime-1,1), ...
                       'flw_in_vph',nan(numtime-1,1), ...
                       'flw_out_veh',nan(numtime-1,1), ...
                       'flw_out_vph',nan(numtime-1,1), ...
                       'dty_veh',nan(numtime,1) , ...
                       'dty_vpm',nan(numtime,1) );
            
            link_ids = obj.scenario_ptr.get_link_ids;
            ind = link_id==link_ids;
            if(~any(ind))
                return;
            end
            
            % flow
            X.flw_in_veh = zeros(numtime-1,1);        
            X.flw_out_veh = zeros(numtime-1,1);      
            for i=1:length(obj.outflow_veh)
                X.flw_in_veh = X.flw_in_veh + obj.inflow_veh{i}(:,ind);
                X.flw_out_veh = X.flw_out_veh + obj.outflow_veh{i}(:,ind);
            end
            X.flw_in_vph = X.flw_in_veh/((obj.time(2)-obj.time(1))/3600);
            X.flw_out_vph = X.flw_out_veh/((obj.time(2)-obj.time(1))/3600);
            
            % density
            X.dty_veh = zeros(numtime,1);
            link_lengths = obj.scenario_ptr.get_link_lengths('us');
            for i=1:length(obj.density_veh)
                X.dty_veh = X.dty_veh + obj.density_veh{i}(:,ind);
            end            
            X.dty_vpm = X.dty_veh/link_lengths(ind);

        end
       
        function [spd,link_ids]=compute_speed(obj)
            if(~obj.sim_output_loaded)
                error('no simulation data.')
            end
            
            vtypes = obj.scenario_ptr.get_vehicle_types();
            numvt = length(vtypes);
            
            % ordered freeway indices
            lintypes = obj.scenario_ptr.get_link_types();
            ordered_ind = extract_linear_fwy_indices(obj.scenario_ptr);
            isfwy = strcmpi(lintypes(ordered_ind),'freeway');
            ordered_fwy_ind = ordered_ind(isfwy);
            link_ids = obj.scenario_ptr.get_link_ids(ordered_fwy_ind);
            
            % initialize
            numlinks = length(ordered_fwy_ind);
            density = zeros(length(obj.time),numlinks);
            flow = zeros(length(obj.time)-1,numlinks);
            
            % totals
            for v=1:numvt
                density = density+obj.density_veh{v}(:,ordered_fwy_ind);
                flow = flow+obj.outflow_veh{v}(:,ordered_fwy_ind);
            end
            
            % units
            outdt = (obj.time(2)-obj.time(1))/3600;
            flow = flow/outdt;
            
            conv = unit_conversion_factor(obj.scenario_ptr.get_units(),'us');
            
            for i=1:numlinks
                link_ind = ordered_fwy_ind(i);
                lgth = obj.scenario_ptr.scenario.NetworkSet.network(1).LinkList.link(link_ind).ATTRIBUTE.length;
                lgth = lgth*conv.length;
                density(:,i) = density(:,i)/lgth;
            end
            
            % remove final density
            density = density(1:end-1,:);
            
            % compute speed
            spd = flow./density;
            ind = density<0.01;
            spd(ind) = 60;
        end
        
        function [D] = get_total_source_flw(obj)
            if(~obj.sim_output_loaded)
                error('run a simulation first.')
            end
            link_ids = obj.scenario_ptr.get_link_ids;
            is_source_link = obj.scenario_ptr.is_source_link;
            source_link_ids = link_ids(is_source_link);
            outdt = (obj.time(2)-obj.time(1))/3600;
            for i=1:length(source_link_ids)
                lind = source_link_ids(i)==link_ids;
                if(~any(lind))
                    continue
                end
                for j=1:length(obj.outflow_veh)
                    f = obj.outflow_veh{j}(:,lind)/outdt;
                    if(i==1 && j==1)
                        flw = f;
                    else
                        flw = flw+f;
                    end
                end
            end
            D.flw = flw;
            D.time = obj.time(1:end-1);
            D.links = source_link_ids;
        end
        
        % travel times
        function [TT] = get_travel_times(obj)
            TT = obj.compute_travel_times();
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  plots and reports                                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        % scenario report ...................................
        function [obj]=report_scenario(obj,pptfile)
            if(~obj.config_loaded)
                error('load a configuration first.')
            end
            obj.scenario_ptr.report(pptfile);
        end
        
        % simulation report ...................................   
        function [obj]=report_simulation(obj,pptfile)

            % open ppt
            [ppt,op]=openppt(pptfile,true);

            addslideTitle(op,{'Simulation report'});

%             addslideTitle(op,{'Contour plots'});
%             X = obj.plot_freeway_contour;
%             for i=1:length(X.figs)
%                 figure(X.figs(i))
%                 addslide(op,'');
%             end
            
            addslideTitle(op,{'Performance'});
            X = obj.plot_performance;
            for i=1:length(X.figs)
                figure(X.figs(i))
                addslide(op,'');
            end
            
            closeppt(ppt,op)
            close all

        end

        function [obj]=plot_source_queues_and_demands(obj,pptfile)
            
            to_ppt = nargin>=2;
                
            if(to_ppt)
                [ppt,op]=openppt(pptfile,true);
            end

            if(~obj.sim_output_loaded)
                error('run a simulation first.')
            end
            
            link_ids = obj.scenario_ptr.get_link_ids;
            is_source_link = obj.scenario_ptr.is_source_link;
            source_link_ids = link_ids(is_source_link);
            l2d = obj.scenario_ptr.get_link_demand_map;
            
            
            cnv = unit_conversion_factor(obj.scenario_ptr.get_units,'US');
            
            for i=1:length(source_link_ids)
                
                % get demand index
                source_link_id = source_link_ids(i);
                
                lind = link_ids==source_link_id;
                if(~any(lind))
                    warning('bad link id: %d',source_link_id)
                end
                
                try
                    demand = l2d(source_link_id);
                catch 
                    warning('no demand found for link %d',source_link_id)
                    continue
                end
                
                % simulated values    
                sim = obj.get_output_for_link_id(source_link_id);
                
                if(to_ppt)
                    figure('Visible','off')
                else
                    figure('Visible','on')
                end
                
                subplot(211)
                plot(demand.time/3600,demand.value*demand.knob*cnv.flw,'k','LineWidth',3)
                hold on
                plot(sim.time_in_sec(1:end-1)/3600,sim.flw_out_vph,'r')
                grid
                set(gca,'XLim',[0 sim.time_in_sec(end)/3600])
                xlabel('time [hr]')
                ylabel('[veh/hr]')
                legend('demand','out flow')
                
                tit=sprintf('Link %d, demand profile %d',source_link_id,demand.dp_id);
                    
                if(~to_ppt)
                    title(tit)
                end

                subplot(212)
                plot(sim.time_in_sec/3600,sim.dty_veh,'k','LineWidth',2)
                set(gca,'XLim',[0 sim.time_in_sec(end)/3600])
                grid
                xlabel('time [hr]')
                ylabel('[veh]')
                legend('queue')
                
                if(to_ppt)
                    addslide(op,tit)
                    close
                end
            end
            
            if(to_ppt)
                closeppt(ppt,op)
                disp(['Wrote to ' pptfile])
            end
            
        end
        
        function [obj]=plot_source_queues_contour(obj)
            if(~obj.sim_output_loaded)
                error('run a simulation first.')
            end
            
            link_ids = obj.scenario_ptr.get_link_ids;
            source_link_ids = link_ids(obj.scenario_ptr.is_source_link);
            
            for i=1:length(source_link_ids)
                sim = obj.get_output_for_link_id(source_link_ids(i));
                queue(i,:) = sim.dty_veh;
            end
            
            figure
            imagesc(queue)
            colorbar
            
        end
        
        function [obj]=plot_sink_flows(obj,pptfile)
            
            to_ppt = nargin>=2;
                
            if(to_ppt)
                [ppt,op]=openppt(pptfile,true);
            end

            if(~obj.sim_output_loaded)
                error('run a simulation first.')
            end
            
            link_ids = obj.scenario_ptr.get_link_ids;
            sink_link_ids = link_ids(obj.scenario_ptr.is_sink_link);
            
            for i=1:length(sink_link_ids)

                sim = obj.get_output_for_link_id(sink_link_ids(i));
                                
                if(to_ppt)
                    figure('Visible','off')
                else
                    figure('Visible','on')
                end
                
                plot(sim.time_in_sec(1:end-1)/3600,sim.flw_out_vph,'b','LineWidth',2)
                grid
                set(gca,'XLim',[0 sim.time_in_sec(end)/3600])
                xlabel('time [hr]')
                ylabel('[veh/hr]')
                
                tit=sprintf('Link %d',sink_link_ids(i));
                if(to_ppt)
                    addslide(op,tit)
                    close
                else
                    title(tit)
                end
            end
            
            if(to_ppt)
                closeppt(ppt,op)
                disp(['Wrote to ' pptfile])
            end
            
        end
        
        % contour plot ........................................
        function [X]=plot_freeway_contour(obj,varargin)
            % argument is a structure with fields
            % pptfile .. name of powerpoint for output, empty for screen only, default empty
            % orientation ... e.g. down_right
            % visible .. ['on']|'off'
            % vds_only ... true|[false]
            % numtimeticks ... int, default 13
                    
            if(~obj.sim_output_loaded)
                error('no simulation data.')
            end
            
            default_options = struct('pptfile','', ...
                                     'orientation','down_right', ...
                                     'visible','on', ...
                                     'vds_only',false, ...
                                     'numtimeticks',13 );
            
            if(length(varargin)>=1)
                in_options = varargin{1};
            else 
                in_options = struct();
            end
            
            if(~isstruct(in_options))
                error('bad input')
            end
               
            df=fieldnames(default_options);
            f =fieldnames(in_options); 
            for i=1:length(df)
                if ismember(df{i},f)
                    options.(df{i}) = in_options.(df{i});
                else
                    options.(df{i}) = default_options.(df{i});
                end
            end
            clear f df in_options

            % vehicle types
            vtypes = obj.scenario_ptr.get_vehicle_types();
            numvt = length(vtypes);
            
            % ordered freeway indices
            lintypes = obj.scenario_ptr.get_link_types();
            ordered_ind = extract_linear_fwy_indices(obj.scenario_ptr);
            ordered_fwy_ind = ordered_ind(  strcmpi(lintypes(ordered_ind),'freeway') );
            clear linktypes ordered_ind
            if options.vds_only
                ordered_fwy_linkids = obj.scenario_ptr.get_link_ids(ordered_fwy_ind);
                vds_info = obj.scenario_ptr.get_ordered_vds;
                vds_info(~strcmpi({vds_info.link_type},'freeway')) = [];
                vds_link_ids = [vds_info.link_id];
                ordered_fwy_ind = ordered_fwy_ind(index_into(vds_link_ids,ordered_fwy_linkids));
                ordered_vds = [vds_info.sensor_vds];

                clear vds_info vds_link_ids ordered_fwy_linkids
            else

            end
            
            % initialize
            numtime = length(obj.time)-1;   % -1 to accomodate CC-sim runs
            numlinks = length(ordered_fwy_ind);
            density = zeros(numtime,numlinks);
            flow = zeros(numtime,numlinks);
            
            % totals
            for v=1:numvt
                density = density+obj.density_veh{v}(1:numtime,ordered_fwy_ind);
                flow = flow+obj.outflow_veh{v}(1:numtime,ordered_fwy_ind);
            end
            
            % units
            outdt = (obj.time(2)-obj.time(1))/3600;
            
            if(outdt==0)
                warning('outdt=0, setting to 300/3600')
                outdt=300/3600;
            end
            
            flow = flow/outdt;
            conv = unit_conversion_factor(obj.scenario_ptr.get_units(),'us');
                   
            for i=1:numlinks
                link_ind = ordered_fwy_ind(i);
                lgth = obj.scenario_ptr.scenario.NetworkSet.network(1).LinkList.link(link_ind).ATTRIBUTE.length;
                lgth = lgth*conv.length; % in miles
                density(:,i) = density(:,i)/lgth;
            end

            to_ppt = ~isempty(options.pptfile);
            
            if(to_ppt)
                [ppt,op]=openppt(options.pptfile,true);
            end
            
            % plot density (dont plot density of first link)
            fig1 = figure('Position',[158 138 940 503],'Visible',options.visible);
            obj.drawcontour(fig1,density,options.orientation,options.numtimeticks);
            colorbar
            
            % label space axis with vds
            if options.vds_only & strcmp(options.orientation,'right_up')
                set(gca,'XTick',0.5+(1:length(ordered_vds)))
                set(gca,'XTickLabel',num2str(ordered_vds'))
            end
            
            if(to_ppt)
                addslide(op,'Densities [veh/mile]')
            else
                title('Densities [veh/mile]')
            end
            
            % plot flow
            fig2 = figure('Position',[158 138 940 503],'Visible',options.visible);
            obj.drawcontour(fig2,flow,options.orientation,options.numtimeticks);
            colorbar
            
            % label space axis with vds
            if options.vds_only & strcmp(options.orientation,'right_up')
                set(gca,'XTick',0.5+(1:length(ordered_vds)))
                set(gca,'XTickLabel',num2str(ordered_vds'))
            end
            
            if(to_ppt)
                addslide(op,'Flows [veh/hr]')
            else
                title('Flows [veh/hr]')
            end
            
            % compute speed
            speed = flow./density;
            ind = density<0.01;
            speed(ind) = 60;
            
            % BIG HACK!!!!!
            speed = min(speed,ones(size(speed))*65);
            
            % plot speed (except first link)
            fig3 = figure('Position',[158 138 940 503],'Visible',options.visible);
            obj.drawcontour(fig3,speed,options.orientation,options.numtimeticks);
            colorbar
            
            % label space axis with vds
            if options.vds_only & strcmp(options.orientation,'right_up')
                set(gca,'XTick',0.5+(1:length(ordered_vds)))
                set(gca,'XTickLabel',num2str(ordered_vds'))
            end
            
            if(to_ppt)
                addslide(op,'Speeds [mile/hr]')
            else
                title('Speeds [mile/hr]')
            end
            
            if(to_ppt)
                closeppt(ppt,op)
            end            
            
            if(nargout>0)
                X.flw = flow;
                X.dty = density;
                X.spd = speed;
                link_ids = obj.scenario_ptr.get_link_ids;
                X.link_ids = link_ids(ordered_fwy_ind);
                X.figs = [fig1 fig2 fig3];
            end
                
            
        end
        
        function [X]=compute_performance(obj,link_mask)
            % link_mask is a logical array with length = total number of
            % links. this defines the subnetwork on which to compute
            % performance

            X = [];
            
            if(~obj.sim_output_loaded)
                warning('no simulation data.')
                return
            end
            
            % vehicle types
            vtypes = obj.scenario_ptr.get_vehicle_types();
            numvt = length(vtypes);
            
            % initialize
            numlinks_total = obj.scenario_ptr.get_num_links();
            
            if(nargin<2)
                link_mask = true(1,numlinks_total);
            end
            
            if(length(link_mask)~=numlinks_total)
                error('the length of the link mask should equal the total number of links')
            end
            clear numlinks_total
            
            numlinks = sum(link_mask);
            
            dty_veh = zeros(length(obj.time),numlinks);
            flw_veh = zeros(length(obj.time)-1,numlinks);
            
            % totals
            for v=1:numvt
                dty_veh = dty_veh+obj.density_veh{v}(:,link_mask);
                flw_veh = flw_veh+obj.outflow_veh{v}(:,link_mask);
            end
            
            % link length 
            meters_per_mile = 1609.34;
            X.link_length_miles = obj.scenario_ptr.get_link_lengths('si')/meters_per_mile;
            X.link_length_miles = X.link_length_miles(link_mask);
            
            % vehicle hours and miles
            X.time = obj.time(1:end-1);
            numtime = length(X.time);
            X.dt_hr = (obj.time(2)-obj.time(1))/3600;

            % total vehicle hours
            X.tot_veh = sum(dty_veh,2);
            X.tot_veh = X.tot_veh(1:end-1);
            
            % total vehicle miles
            X.tot_flux = sum(flw_veh.*repmat(X.link_length_miles,numtime,1),2)/X.dt_hr;   % veh.mile/hr
        end
        
        % performance plot ...................................
        function [X]=plot_performance(obj,pptfile)

            to_ppt = nargin>=2;
                
            if(to_ppt)
                [ppt,op]=openppt(pptfile,true);
            end
            
            P=obj.compute_performance();
            
            % total vehicle hours          
            if(to_ppt)
                X.figs(1) = figure('Visible','off');
            else
                X.figs(1) = figure('Visible','on');
            end
                
            obj.timeplot(P.time,P.tot_veh,5,'Total vehicles [veh]',X.figs(1))
            tit = sprintf('TVH = %.0f veh.hr',sum(P.tot_veh)*P.dt_hr);
                        
            if(to_ppt)
                addslide(op,tit)
                close
            else
                title(tit)
            end
                
            % total vehicle miles  
            if(to_ppt)
                X.figs(2) = figure('Visible','off');
            else
                X.figs(2) = figure('Visible','on');
            end
            
            obj.timeplot(P.time,P.tot_flux,5,'Total flux [veh.miles/hr]',X.figs(2))
            tit = sprintf('TVM = %.0f veh.mile',sum(P.tot_flux)*P.dt_hr);       
            if(to_ppt)
                addslide(op,tit)
                close
            else
                title(tit)
            end
            
            % average speed
            if(to_ppt)
                X.figs(3) = figure('Visible','off');
            else
                X.figs(3) = figure('Visible','on');
            end
            
            obj.timeplot(P.time,P.tot_flux./P.tot_veh,5,'Weighted speed [mile/hr]',X.figs(3))
            avg_speed = P.tot_flux./P.tot_veh;
            tit = sprintf('Average experienced speed = %.0f mile/hr',mean(avg_speed(~isinf(avg_speed))));                   
            if(to_ppt)
                addslide(op,tit)
                close
            else
                title(tit)
            end
            
            if(to_ppt)
                closeppt(ppt,op)
                disp(['Wrote to ' pptfile])
            end
            
        end
        
        % link state plot ....................................
        function [t,d,f]=plot_link_state(obj,link_id)
           
            if(~obj.sim_output_loaded)
                error('no simulation data.')
            end
            
            link_ids = obj.scenario_ptr.get_link_ids();
            ind = link_id==link_ids;
            
            if(~any(ind))
                error('bad link id')
            end
           
            numtime = length(obj.time)-1; 
            
            d = zeros(numtime+1,1);
            f = zeros(numtime,1);
            for i=1:length(obj.density_veh)
                d = d+obj.density_veh{i}(:,ind);
                f = f+obj.outflow_veh{i}(:,ind);
            end
            d = d(1:end-1);
            t = obj.time(1:end-1);
                        
            %new_figure(560,610);
            figure
            
            subplot(311)
            plot(t,d,'LineWidth',2)
            xlabel('time')
            ylabel('density [veh]')
            grid
            
            title(['Link id ' num2str(link_id)])
            
            subplot(312)
            plot(t,f,'LineWidth',2)
            xlabel('time')
            ylabel('flow [veh]')
            grid
            
            subplot(313)
            plot(t,f./d,'LineWidth',2)
            xlabel('time')
            ylabel('speed [-]')
            grid
                        
            
        end
        
        % demands ............................................
        function [dem,sim] = plot_total_demands(obj)

            sim = obj.get_total_source_flw();
            dem = obj.scenario_ptr.get_total_source_flw();

            figure
            plot(dem.time,dem.flw,'k','LineWidth',2), hold on
            plot(sim.time,sim.flw,'r--','LineWidth',2)
            grid
            legend('demand','simulated flow')
            hold on
        end
        
        % compare with data ..................................
        function [test] = sim_vs_pems(obj,day,processed_folder,excel_filename,ppt_filename,exclude_from_fhwa,good_vds)
            
            if ~obj.sim_output_loaded || ~obj.simulation_done
                error('run the simulation first')
            end
            
            if length(obj.density_veh)>1
                warning('aggregating over vehicle types')
            end
            
            write_to_excel = nargin>=4;
            write_to_ppt = nargin>=5;
            given_good_sensors = nargin>=7;
            if nargin<6
                exclude_from_fhwa = [];
            end
            
            LOOP_HEALTH_THRESHOLD = 90;
                       
            % check that the simulation was run
            if ~obj.sim_output_loaded 
                error('load simulation data first')
            end
            
            % sensor table
            sensor_table = obj.scenario_ptr.get_sensor_table;
            
            % keep only those that are attached to the network
            sensor_table = sensor_table( logical(sensor_table.is_attached) , :);
            
            % keep only those that are mainline sensors
            ind1 = strcmpi(sensor_table.link_type,'Freeway');
            ind2 = strcmpi(sensor_table.sensor_link_type,'ML');
            if ~all(ind1 == ind2)
                error('disagreement between link_type and sensor_link_type')
            end
            sensor_table = sensor_table(ind1,:);
            clear ind1 ind2

            % keep only those that are good
            if given_good_sensors
                sensor_table = sensor_table( ismember(good_vds,sensor_table.vds),:);   
            else
                beats_vds = sensor_table.vds;
                
                % load health information
                load(fullfile(processed_folder,'loop_health'))
                if ~all(ismember(beats_vds,loop_health.vds))
                    error('do not have loop health information for all vds')
                end
                if ~ismember(day,loop_health.days)
                    error('do not have loop health information for that day')
                end
                
                is_good = loop_health.percent_observed( ...
                    index_into(beats_vds,loop_health.vds) , ...
                    loop_health.days==day) > LOOP_HEALTH_THRESHOLD;
                            
                sensor_table = sensor_table( is_good ,:);
                clear is_good loop_health beats_vds
            end

            % record and order the remaining vds
            ordered = make_ordered_vds(obj.scenario_ptr);
            ind = ismember([ordered.sensor_vds],sensor_table.vds);
            ordered_vds = [ordered(ind).sensor_vds];
            ordered_links = [ordered(ind).link_id];
            clear sensor_table ordered ind
            
            % load traffic data
            pems = PeMS5minData;
            pems.load(processed_folder,ordered_vds,day);
            if any( cellfun( @(x) isempty(x) , pems.data) )
                error('did not find data for some vds on this day')
            end
            pems_data = pems.get_data_batch_aggregate(ordered_vds,day,'var',{'flw','occ','spd'},'fill',true,'smooth',true);
            clear pems
            % HACK: correct for 5 minute displacement
            pems_data.flw = circshift(pems_data.flw,-2);
            pems_data.occ = circshift(pems_data.occ,-2);
            pems_data.spd = circshift(pems_data.spd,-2);
            
            % aggregate pems data
            a = repmat({ones(1,12)},1,24);
            fivemin_to_hr = blkdiag(a{:})/12;
            clear a
            pems_flw_hr = fivemin_to_hr*pems_data.flw;
            pems_spd_hr = fivemin_to_hr*pems_data.spd;
            
            % collect simulation outputs
            beats_flw_5min = zeros(288,length(ordered_links));
            beats_dty_5min = zeros(288,length(ordered_links));
            for i=1:length(ordered_links)
                X = obj.get_output_for_link_id(ordered_links(i));
                beats_flw_5min(:,i) = X.flw_out_vph;
                beats_dty_5min(:,i) = X.dty_vpm(1:end-1,:);
            end
            beats_spd_5min = beats_flw_5min./beats_dty_5min;
            beats_spd_5min(beats_spd_5min>70) = 70;
            beats_spd_5min(isnan(beats_spd_5min)) = 70;
            clear X i
            
            % aggregate beats data
            beats_flw_hr = fivemin_to_hr*beats_flw_5min;
            beats_spd_hr = fivemin_to_hr*beats_spd_5min;
            
            % individual link flows
            flw_range = zeros(size(pems_flw_hr));
            flw_range(pems_flw_hr>700 & pems_flw_hr<2700) = 1;
            flw_range(pems_flw_hr>2700) = 2;
            
            flw_error = abs(beats_flw_hr-pems_flw_hr);
            flw_error_pct = flw_error./pems_flw_hr;
           
%             % keep only reportable vds
%             report_vds =  true(1,32);
%             report_vds(ismember(ordered_vds,exclude_from_fhwa)) = false;
%             flw_range = flw_range(:,report_vds);
%             flw_error = flw_error(:,report_vds);
%             flw_error_pct = flw_error_pct(:,report_vds);
% 
%             % remove first hour 
%             flw_range = flw_range(2:end,:);
%             flw_error = flw_error(2:end,:);
%             flw_error_pct = flw_error_pct(2:end,:);
            % same for ...
            report_beats_flw = beats_flw_hr; %(2:end,report_vds);
            report_pems_flw = pems_flw_hr; %(2:end,report_vds);
            % evaluate link flows tests
            passes = false(size(flw_range));
            passes(flw_range==0) = flw_error(flw_range==0)<100;
            passes(flw_range==1) = flw_error_pct(flw_range==1)<0.15;
            passes(flw_range==2) = flw_error(flw_range==2)<400;
            
            
            test(1).name = 'individual link flows, low';
            test(1).result = sum(sum(passes(flw_range==0)))/sum(sum(flw_range==0));
            test(1).passes = test(1).result>0.85;
            
            test(2).name = 'individual link flows, medium';
            test(2).result = sum(sum(passes(flw_range==1)))/sum(sum(flw_range==1));
            test(2).passes = test(2).result>0.85;
            
            test(3).name = 'individual link flows, high';
            test(3).result = sum(sum(passes(flw_range==2)))/sum(sum(flw_range==2));
            test(3).passes = test(3).result>0.85;
            
            test(4).name = 'sum of all link flows';
            total_beats_flw = sum(sum(report_beats_flw));
            total_pems_flw = sum(sum(report_pems_flw));
            test(4).result = abs(total_pems_flw-total_beats_flw)/total_pems_flw;
            test(4).passes = test(4).result<0.05;
            
            test(5).name = 'GEH individual link flows';
            gehLessthanFive = BeatsSimulation.GEH(report_beats_flw,report_pems_flw) < 5;
            test(5).result = sum(sum(gehLessthanFive))/numel(gehLessthanFive);
            test(5).passes = test(5).result> 0.85;
            
            test(6).name = 'GEH sum link flows';
            test(6).result = BeatsSimulation.GEH(total_beats_flw,total_pems_flw);
            test(6).passes = test(6).result<4;
            
            %  travel times
            beats_travel_times = obj.compute_travel_times(beats_spd_hr,ordered_links);
            pems_travel_times = obj.compute_travel_times(pems_spd_hr,ordered_links);
            tt_error_min = abs(beats_travel_times-pems_travel_times)*60;
            tt_error_pct = abs(beats_travel_times-pems_travel_times)./pems_travel_times;
            short_trips = pems_travel_times*60 < 1;

            test(7).name = 'link travel times (should be journey times) (short)';
            num_short = sum(sum(short_trips));
            if num_short>0
                test(7).result = sum(tt_error_min(short_trips)<1) / num_short;
            test(7).passes = test(7).result>0.85;
            else
                test(7).result = 1;
                test(7).passes = true;
            end

            test(8).name = 'link travel times (should be journey times) (long)';
            num_long = sum(sum(~short_trips));
            if num_long>0
                test(8).result = sum(tt_error_pct(~short_trips)<0.15) / num_long;
            test(8).passes = test(8).result>0.85;
            else
                test(8).result = 1;
                test(8).passes = true;
            end
            if write_to_excel
                here = fileparts(mfilename('fullpath'));
                copyfile(fullfile(here,'fhwa_template.xls'),excel_filename,'f');
                xlswrite(excel_filename,[test(1:6).result]','E4:E9');
                xlswrite(excel_filename,[test(7:8).result]','E11:E12');
            end
            
            if write_to_ppt
                lower_bnd = 0*pems_flw_hr;
                lower_bnd(flw_range==0) = pems_flw_hr(flw_range==0)-100;
                lower_bnd(flw_range==1) = pems_flw_hr(flw_range==1)*0.85;
                lower_bnd(flw_range==2) = pems_flw_hr(flw_range==2)-400;
                upper_bnd = 0*pems_flw_hr;
                upper_bnd(flw_range==0) = pems_flw_hr(flw_range==0)+100;
                upper_bnd(flw_range==1) = pems_flw_hr(flw_range==1)*1.15;
                upper_bnd(flw_range==2) = pems_flw_hr(flw_range==2)+400;
                is_good = beats_flw_hr>=lower_bnd & beats_flw_hr<=upper_bnd;
                ylim = round(max(max([upper_bnd beats_flw_hr pems_flw_hr]))/100)*100;
                [ppt,op]=openppt(ppt_filename,true);
                addslideTitle(op,{'simulated vs measured flows'});
                for i=1:length(ordered_links)
                    figure('Visible','off')
                    x = 0:23;
                    b = beats_flw_hr(:,i);
                    p = pems_flw_hr(:,i);
                    l = lower_bnd(:,i)';
                    u = upper_bnd(:,i)';
                    g = is_good(:,i);
                    jbfill(gcf,(0:23),l,u,'c','c',true,0.4);
                    hold on
                    z=plot(x,[b p],'-','LineWidth',2);
                    plot(x(g),b(g),'bo','LineWidth',2);
                    plot(x(~g),b(~g),'rx','LineWidth',2);
                    set(gca,'XLim',[0 23]);
                    set(gca,'YLim',[0 ylim]);
                    hline(700,'--')
                    hline(2700,'--')
                    grid
                    xlabel('time [hr]')
                    ylabel('flow [vph]')
                    legend(z,'beats','pems','Location','best')
                    addslide(op,sprintf('link id = %d, vds = %d',ordered_links(i),ordered_vds(i)))
                    close
                end
                closeppt(ppt,op)
            end
        
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  private methods                                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = private)
        
        function [obj] = reset_object(obj)
            obj.scenario_ptr = ScenarioPtr;
            obj.configfile = '';
            obj.config_loaded = false;
            obj.sim_output_loaded = false;
            obj.reset_simdata();
        end
        
        function [obj] = reset_simdata(obj)
        	obj.time = [];
        	obj.density_veh = [];
        	obj.outflow_veh = [];
        	obj.inflow_veh = [];
            obj.sim_output_loaded = false;
            obj.simulation_done = false;
        end
        
        function [obj]=drawcontour(obj,fig,Z,orientation,numtimeticks)
            
            if nargin<5
                numtimeticks = 20;
            end
            
            figure(fig)
            
            % draw the contour
            Zx = [Z Z(:,1);Z(1,:) nan];
            Zx(isnan(Zx)) = -1;
            h=pcolor(Zx);
            set(h,'EdgeAlpha',0)
            
            % coloring
%             C = colormap();
%             C = [[0.5 0 0];C];
%             colormap(C);
            colorbar
            
            [timetick,timeticklabel]=obj.get_time_ticks(numtimeticks);
            
            set(gca,'YTick',timetick)
            set(gca,'YTickLabel',timeticklabel)
            
            % orientation
            switch(orientation)
                case {'right_up','PeMS','MM'}
                    view(90,270)
                case {'down_right','FREQ'}
                    view(0,270)
                case {'up_right','TOPL','Aimsun'}
                    view(0,90)
                case 'right_down'
                    view(90,90)
            end
            
        end
        
        function []=timeplot(obj,X,Y,numticks,ylab,fig)
            [timetick,timeticklabel,f]=obj.get_time_ticks(numticks);
            
            if(nargin<6)
                fig = figure;
            end
            figure(fig)
            fillplot(X,Y,{'b'},gcf);
            hold on
            plot(X,Y,'b','LineWidth',2)
            set(gca,'XTick',obj.time(timetick))
            set(gca,'XTickLabel',timeticklabel)
            set(gca,'XLim',[X(1) X(end)])
            grid
            xlabel(['Time ' f])
            ylabel(ylab)
        end
        
        function [timetick,timeticklabel,f]=get_time_ticks(obj,numtick)
            numtime = length(obj.time);
            timetick = unique(round(linspace(1,numtime,min(numtick,numtime))));
            duration = obj.time(end)-obj.time(1);
            if(duration>3600)
                f = 'HH:MM';
            else
                f = 'MM:SS';
            end
            timeticklabel = datestr(obj.time(timetick)/86400,f);
        end
        
        function [TT] = compute_travel_times(obj,given_speeds,give_links_ids)
            
            if nargin<2
                [spd,link_ids] = obj.compute_speed;
            else
                spd = given_speeds;
                link_ids = give_links_ids;
            end
            
            if size(spd,2) ~= length(link_ids)
                error('inconsistent dimensions')
            end
            
            link_ind = index_into(link_ids,obj.scenario_ptr.get_link_ids);
            if any(link_ind==0)
                error('bad link id')
            end
            
            link_lengths = obj.scenario_ptr.get_link_lengths('us');
            TT = ( ones(size(spd,1),1) * link_lengths(link_ind) ) ./ spd;
        end
        
        function [param,prop_file] = create_properties_file(obj,param)
                        
            if(~isfield(param,'SIM_DT'))
                warning('Missing SIM_DT, using to 5 seconds.')
                param.SIM_DT = 5;
            end
            obj.sim_dt = param.SIM_DT;
            
            if(~isfield(param,'OUTPUT_DT'))
                param.OUTPUT_DT = 300;
            end
            obj.out_dt = param.OUTPUT_DT;
            
            if(isfield(param,'OUTPUT_PREFIX'))
                obj.temp_output = '';
            else
                obj.temp_output = tempname;
                param.OUTPUT_PREFIX = obj.temp_output;
            end
            
            % write properties file
            prop_file = sprintf('%s.properties',tempname);
            fid = fopen(prop_file,'w+');
            S = fieldnames(param);
            for i = 1:length(S)
                if(ischar(param.(S{i})))
                    fprintf(fid,'%s = %s\n',S{i},regexprep(param.(S{i}),'\\','\\\\'));
                else
                    fprintf(fid,'%s = %f\n',S{i},param.(S{i}));
                end
            end
            fclose(fid);
        end

    end
    
    methods (Static, Access = public)
        
        function [obj] = import_beats_classes(obj,beats_path)
            
            if (nargin<2)
                beats_path=fileparts(mfilename('fullpath'));
            end    

            warning('off','MATLAB:Java:DuplicateClass')
            javaaddpath(fullfile(beats_path,'beats-0.1-SNAPSHOT-jar-with-dependencies.jar'));
            
%             % pull in the BeATS java classes
%             javaaddpath(fullfile( beats_path, 'target', 'classes') );
%             
%             % Pull in BeATS lib jars
%             % returns a structure with info about the files in the
%             % directory
%             lib_jars = dir(fullfile(beats_path,'target','lib','*.jar'));
%             
%             % need to create a cell array of strings of the filepath to
%             % concatenate with lib_jars' entries
%             [prefix{1:numel(lib_jars)}] = deal(fullfile(beats_path,'target','lib'));
%             
%             % cellfun call just applies fullfile to each similarly-indexed
%             % pair in the cell arrays prefix and {lib_jars.name}
%             javaaddpath( cellfun(@fullfile, prefix, {lib_jars.name},'UniformOutput',false) );
            
        end
        
        function  [props] = read_properties_file(prop_file)
            props=repmat(struct('name','','value',''),1,0);
            fid=fopen(prop_file);
            if(fid<0),return,end
            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                
                % remove comments
                ind=strfind(tline,'#');
                if(~isempty(ind))
                    tline=tline(1:min(ind)-1);
                end
                if(isempty(tline)),continue,end
                
                x=strsplit(tline,'=');
                if(length(x)<2),continue,end
                
                props(end+1)=struct('name',strtrim(x{1}),'value',strtrim(x{2}));                
            end
            fclose(fid);
        end
        
        function [x] = GEH(A,B)
           x = sqrt( 2*( (A-B).^2 )./(A+B) ); 
        end
        
    end
    
end
