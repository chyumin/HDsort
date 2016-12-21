classdef Sorting < handle
    properties
        DS
        dPath
        name;
        
        files
        filesLoaded
        
        groupFolder
        groupFilesList
        
        %gridType
        sortjob
        postprocessjob
    end
    
    methods
        % -----------------------------------------------------------------
        function self = Sorting(DS, dPath, name, varargin)
            
            P.maxElPerGroup = 9;
            P = util.parseInputs(P, varargin, 'error');
            %self.gridType = P.gridType;
            
            self.DS = DS;
            self.dPath = dPath;
            self.name = name;
            self.init();
        end
        
        % -----------------------------------------------------------------
        function init(self)
            is_letter = isstrprop(self.name, 'alpha');
            
            % Create group file and groups folder
            self.files.groupFile     = fullfile(self.dPath, ['sort_' self.name], 'groupFile.mat');
            self.groupFolder         = fullfile(self.dPath, ['sort_' self.name], 'groups');
            [dir_exists,mess,messid] = mkdir( self.groupFolder );
            assert(dir_exists, 'Output directory could not be created!');
            self.files.GroupStruct   = fullfile(self.groupFolder, 'G_struct.mat');
            
            self.filesLoaded.groupFile     = [];
            self.filesLoaded.GroupStruct   = [];
            self.filesLoaded.sortingResult = [];
            
            keySet = {'spikes_det',
                'spikes_det_merged',
                'spikes_cut',
                'spikes_cut_bin',
                'spikes_aligned',
                'cov',
                'spikes_prewhitened',
                'spikes_features',
                'clusters_meanshift',
                'botm_matching',
                'clusters_meanshift_merged',
                'P',
                'templates'};
            valueSet =  {'.020spikes_det.mat',
                '.030spikes_det_merged.mat',
                '.040spikes_cut.mat',
                '.040spikes_cut.mat_cutSpikes.bin',
                '.050spikes_aligned.mat',
                '.060cov.mat',
                '.070spikes_prewhitened.mat',
                '.080spikes_features.mat',
                '.090clusters_meanshift.mat',
                '.100botm_matching.mat',
                '.110clusters_meanshift_merged.mat',
                '.P.mat',
                '_templates.mat'
                };
            self.groupFilesList = containers.Map(keySet,valueSet);
            
        end
        % -----------------------------------------------------------------
        function N = getNElGroups(self)
            S = self.loadGroupFile();
            N = length(S.groups);
        end
        % -----------------------------------------------------------------
        function groups = getElGroups(self)
            S = self.loadGroupFile();
            groups = S.groups;
        end
        
        
        % -----------------------------------------------------------------
        function createLocalElectrodeGroups(self, maxElPerGroup)
            MES = self.DS.MultiElectrode.toStruct();
            
            try
                load(self.files.groupFile)
                disp('Groups already exist')
            catch
                disp('Create groups...');
                electrodePositions = MES.electrodePositions;
                electrodeNumbers   = MES.electrodeNumbers;
                [groupsidx nGroupsPerElectrode] = mysort.mea.constructLocalElectrodeGroups(electrodePositions(:,1), electrodePositions(:,2), 'maxElPerGroup', maxElPerGroup);
                disp(['Number of groups created: ' num2str(length(groupsidx))])
                
                groups = {};
                for ii= 1:length(groupsidx)
                    groups{ii} = electrodeNumbers(groupsidx{ii});
                end
                save(self.files.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
            end
        end
        
        %% ----------------------------------------------------------------
        function [R, P] = startSorting(self, varargin)
            % No sorting loadable, start sorting
            %P.runMode = 'grid';
            P.sortingMode = 'BSSE'; % Alternatives are: euler, grid_parfor, grid_for, localHDSorting
            
            P.forceExecution = false;
            %P.parfor = true;
            P.dataPath = self.dPath;
            P.queue = 'regular';
            P.maxElPerGroup = 9;
            P = util.parseInputs(P, varargin, 'error');
            
            try
                % Try loading the results. If the sorting was not run yet,
                % this will fail and we start the sorting
                if P.forceExecution
                    % If the user wants to run the sorter even if the
                    % results are already there yet, jump to the catch.
                    error('Jump to catch statement');
                end
                R = self.loadSortingResult();
                
            catch
                
                self.createLocalElectrodeGroups(P.maxElPerGroup);
                
                if ~strcmp(P.sortingMode, 'localHDSorting')
                    %P.sortingMode = 'BSSE';
                    %if strcmp(P.runMode, 'grid') || strcmp(P.runMode, 'grid_locally')
                    [R, P] = self.sortOnGrid(P);
                    
                    %                     assert(isa(self.DS, 'mysort.mea.CMOSMEA'), 'At the moment, the grid framework can only be run with CMOSMEA objects! If you want to test the sorter with something else, use the flag ''runMode'', ''localHDSorting''');
                    %                     sourceFiles = self.DS.getSourceFileNames();
                    %
                    %                     self.sortjob = grid.SortJob(self.name, P.dataPath, sourceFiles, ...
                    %                         'queue', P.queue, 'onEuler', self.onEuler, 'runtime_hours', 5);
                    %
                    %                     self.sortjob.createBOTMGroups();
                    %                     %self.sortjob.copyDataFiles();
                    %                     self.sortjob.setTaskParameters();
                    %                     self.sortjob.prepareTasks();
                    %
                    %                     if strcmp(P.runMode, 'grid')
                    %                         %self.createAutoSubmitToken();%self.sortjob.scratchlocation.folders.main);
                    %                         [nCompleted, tasksNotCompleted, nErrors, tasksWithErrors] = self.sortjob.summarySnapshot(true);
                    %                         if self.sortjob.nTasks > nCompleted
                    %                             self.sortjob.createAutoSubmitToken();
                    %                         else
                    %                             disp('All tasks seem to be completed...')
                    %                         end
                    %
                    %                         all_tasks_completed = self.sortjob.waitForTasksToFinish(60);
                    %
                    %                     elseif strcmp(P.runMode, 'grid_locally')
                    %
                    %                         if P.parfor
                    %                             parfor t = 1:self.sortjob.nTasks()
                    %                                 self.sortjob.runTaskLocally(t);
                    %                             end
                    %                         else
                    %                             for t = 1:self.sortjob.nTasks()
                    %                                 self.sortjob.runTaskLocally(t);
                    %                             end
                    %                         end
                    %                         all_tasks_completed = true;
                    %
                    %                     else
                    %                         error('Something went wrong here!')
                    %                     end
                    %                     self.sortjob.summarizeReports();
                    %                    % self.sortjob.copyBackResults();
                    %
                elseif strcmp(P.sortingMode, 'localHDSorting')
                    % This function includes the postprocessing step and
                    % saves the results directluy to a file:
                    [gdf_merged, T_merged, localSorting, localSortingID] = ...
                        hdsort.startHDSorting(self.DS, self.dPath, self.name);
                    R.gdf_merged = gdf_merged;
                    R.T_merged = T_merged;
                    R.localSorting = localSorting;
                    R.localSortingID = localSortingID;
                end
            end
            
        end
        
        % -----------------------------------------------------------------
        function R = loadSortingResult(self)
            self.files.sortingResult = fullfile(self.dPath, [self.name '_results.mat']);
            if isempty(self.filesLoaded.sortingResult)
                try
                    self.filesLoaded.sortingResult = load(self.files.sortingResult);
                catch
                    error('Could not load sorting result, run sorter first!');
                end
            end
            R = self.filesLoaded.sortingResult;
        end
        
        % -----------------------------------------------------------------
        % Functions for using the computer grid:
        % -----------------------------------------------------------------
        
        function P = prepareSortJob(self, varargin)
            %P.runMode = 'grid';
            %P.forceExecution = false;
            %P.postProcFunc = '';
            %P.parfor = true;
            %P.dataPath = self.dPath;
            
            P = util.parseInputs(struct(), varargin, 'merge');
            
            assert(isa(self.DS, 'mysort.mea.CMOSMEA'), 'At the moment, the grid framework can only be run with CMOSMEA objects! If you want to test the sorter with something else, use the flag ''sortingMode'', ''localHDSorting''');
            
            if ~strcmp(P.sortingMode, 'BSSE') && ~strcmp(P.sortingMode, 'euler') 
                        gridType = 'BSSE';
                    else
                        gridType = P.sortingMode;
                    end
            
            sourceFiles = self.DS.getSourceFileNames();
            self.sortjob = grid.SortJob(self.name, self.dPath, sourceFiles, ...
                'groupFile', self.files.groupFile, ...
                'gridType', gridType, 'queue', P.queue)
            
            %self.sortjob.createBOTMGroups();
            self.sortjob.setTaskParameters();
            %self.sortjob.prepareTasks();
            %self.sortjob.summarizeReports();
        end
        
        % -----------------------------------------------------------------
        function [R, P] = sortOnGrid(self, varargin)
            P = util.parseInputs(struct(), varargin, 'merge');
            
            P = self.prepareSortJob(P);
            %if strcmp(P.runMode, 'grid') || strcmp(P.runMode, 'grid_locally')
            
            
            %if strcmp(P.runMode, 'grid') || strcmp(P.runMode, 'grid_locally')
            %assert(isa(self.DS, 'mysort.mea.CMOSMEA'), 'At the moment, the grid framework can only be run with CMOSMEA objects! If you want to test the sorter with something else, use the flag ''runMode'', ''localHDSorting''');
            %sourceFiles = self.DS.getSourceFileNames();
            
            %self.sortjob = grid.SortJob(self.name, P.dataPath, sourceFiles, ...
            %    'queue', P.queue, 'onEuler', self.onEuler, 'runtime_hours', 5);
            
            %self.sortjob.createBOTMGroups();
            %self.sortjob.copyDataFiles();
            %self.sortjob.setTaskParameters();
            %self.sortjob.prepareTasks();
            
            if strcmp(P.sortingMode, 'BSSE') || strcmp(P.sortingMode, 'euler')
                assert(~strcmp(P.sortingMode, 'euler'), 'Not implemented yet!')
                %self.createAutoSubmitToken();%self.sortjob.scratchlocation.folders.main);
                
                [nCompleted, tasksNotCompleted, nErrors, tasksWithErrors] = self.sortjob.summarySnapshot(true);
                if self.sortjob.nTasks > nCompleted
                    self.sortjob.createAutoSubmitToken();
                else
                    disp('All tasks seem to be completed...')
                end
                
                all_tasks_completed = self.sortjob.waitForTasksToFinish(60);
                
            elseif strcmp(P.sortingMode, 'grid_parfor')
                %if P.parfor
                parfor t = 1:self.sortjob.nTasks()
                    self.sortjob.runTaskLocally(t);
                end
                all_tasks_completed = true;
                
                %end
            elseif strcmp(P.sortingMode, 'grid_for')
                for t = 1:self.sortjob.nTasks()
                    self.sortjob.runTaskLocally(t);
                end
                all_tasks_completed = true;
                
                
            else
                error('Unknown sortingMode!')
            end
            self.sortjob.summarizeReports();
            
        end
        
        
        % -----------------------------------------------------------------
        function [R, P] = postprocessGridSorting(self, varargin)
            % newPostProcFunc is a string that specifies which
            % prostprocessing function should be used (merging of
            % templates).
            % Default is '', i.e. no string is passed on.
            %P0.newPostProcFunc = '';
            %P.newComputed = true;
            P.forceNewPostprocessing = false;
            %P.postProcessOnGrid = true;
            P.postProcessMode = 'BSSE';
            
            P = util.parseInputs(P, varargin, 'error');
            
            try
                R = self.loadSortingResult();
                %P.newComputed = false;
                assert(~P.forceNewPostprocessing, 'Force new postprocessing')
            catch
                assert( ~isempty(self.sortjob), 'Make sure that a SortJob object is created before starting the postprocessing.')
                
                disp('Start postprocessGridSorting...')
                %try
                %load(self.files.sortingResult);
                %assert(~P.forceNewPostprocessing, 'Force new postprocessing')
                %disp('Postprocessing already finished.')
                %catch
                
                if ~strcmp(P.postProcessMode, 'localPostProcessing')
%                if P.postProcessOnGrid
                    
                    if ~strcmp(P.postProcessMode, 'BSSE') && ~strcmp(P.postProcessMode, 'euler') 
                        gridType = 'BSSE';
                    else
                        gridType = P.postProcessMode;
                    end
                    
                    self.postprocessjob = grid.PostprocessJob(self.name, ...
                        self.sortjob.folders.root, self.files.groupFile, ...
                        self.sortjob.folders.groups, self.files.sortingResult, ...
                        'gridType', gridType, 'runtime_hours', 1)
                    
                    self.postprocessjob.setTaskParameters();
                    %self.postprocessjob.prepareTasks();
                    self.postprocessjob.createAutoSubmitToken();
                    all_tasks_completed = self.postprocessjob.waitForTasksToFinish(10);
                    R = load(self.files.sortingResult);
                    
                else
                    disp('Load group information...');
                    GF = load(self.files.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
                    
                    disp('Start postprocessing...');
                    [R, P_] = mysort.HDSorting.processLocalSortings(...
                        self.sortjob.folders.groups,...
                        self.name, GF.groups, GF.groupsidx, ...
                        'groupPaths', self.sortjob.folders.groups);
                    
                    disp('Check and save data...');
                    units = unique(R.gdf_merged(:,1));
                    nU = length(units);
                    assert(length(R.localSorting) == nU, 'must be identical');
                    assert(length(R.localSortingID) == nU, 'must be identical');
                    assert(size(R.T_merged,3) == nU, 'must be identical');
                    
                    disp('Saving postprocessing results...')
                    save(self.files.sortingResult, '-struct', 'R', '-v7.3');
                end
                
                %end
            end
        end
        
        % -----------------------------------------------------------------
        function [R, P] = reuptakeSorting(self, varargin)
            P.dataPath = self.dPath;
            %P.postProcFunc = '';
            P = util.parseInputs(P, varargin, 'error');
            
            try
                R = self.loadSortingResult();
            catch
                try
                    assert( ~isempty(self.sortjob), 'Create new SortJob object...')
                    %all_tasks_completed = self.sortjob.waitForTasksToFinish(60);
                catch
                    P = self.prepareSortJob(P);
                    
                    %sourceFiles = self.DS.getSourceFileNames();
                    
                    %self.sortjob = grid.SortJob(self.name, P.dataPath, sourceFiles, ...
                    %    'gridType', self.gridType, 'runtime_hours', 5);
                    %self.sortjob.createBOTMGroups();
                    %self.sortjob.copyDataFiles();
                    %self.sortjob.setTaskParameters();
                    %self.sortjob.prepareTasks();
                end
                all_tasks_completed = self.sortjob.waitForTasksToFinish(60);
                %end
                
                self.sortjob.summarizeReports();
                %self.sortjob.copyBackResults();
            end
        end
        
        
        
        %% GROUP ORGA -----------------------------------------------------
        function f = getGroupFolder(self, groupNr)
            %f = fullfile(self.dPath, ['sort_' self.name], 'groups', sprintf('group%04d', groupNr));
            f = fullfile(self.groupFolder, sprintf('group%04d', groupNr));
        end
        function keys = listFilesInGroup(self)
            keys = self.groupFilesList.keys;
            disp(keys)
        end
        
        function [fname, does_exist] = getFileNameInGroup(self, groupNr, key)
            foldername = self.getGroupFolder(groupNr);
            fname = fullfile(foldername, ['sort_' self.name self.groupFilesList(key)]);
            does_exist = exist(fname);
        end
        
        function out = loadFileInGroup(self, groupNr, key, out)
            if nargin < 4
                out = struct();
            end
            fname = self.getFileNameInGroup(groupNr, key);
            x = load(fname);
            out = util.mergeStructs(out, x);
        end
        
        % -----------------------------------------------------------------
        function out = deleteFilesInGroups(self, key)
            out = true;
            for groupNr = 1:self.getNElGroups()
                fname = self.getFileNameInGroup(groupNr, key);
                
                try
                    delete(fname);
                catch
                    out = false;
                end
            end
        end
        function out = deleteResultsFile(self)
            out = true;
            try
                delete(self.files.sortingResult);
            catch
                out = false;
            end
            
        end
        
        % LOADING GROUP FILES --------ALL GROUPS---------------------------
        function G = loadGroupFile(self)
            if isempty(self.filesLoaded.groupFile)
                self.filesLoaded.groupFile = load(self.files.groupFile);
            end
            G = self.filesLoaded.groupFile;
        end
        
        % -----------------------------------------------------------------
        function G = loadGroupStruct(self)
            if isempty(self.filesLoaded.GroupStruct)
                self.filesLoaded.GroupStruct = load(self.files.GroupStruct);
            end
            G = self.filesLoaded.GroupStruct;
        end
        % -----------------------------------------------------------------
        function gdfs = loadGroupGDFs(self)
            G = self.loadGroupStruct();
            gdfs = {G.G.gdf};
        end
        
        % LOADING GROUP FILES --------INDIVIDUAL GROUPS--------------------
        function S = loadDetectedSpikes4Group(self, groupNr)
            groupFolder = self.getGroupFolder(groupNr);
            S = load(fullfile(groupFolder, [self.name '.020spikes_det.mat']));
        end
        
        % -----------------------------------------------------------------
        function [wfs, electrodePositions, channelNumbers] = getAllSpikesFromUnitID(self, unitID)
            
            %% Get all the spikes of one neuron:
            % Get the group and the id within the group:
            groupNumber = floor(unitID/1000);
            number_within_group = mod(unitID, 1000);
            
            %% Get the units in each group:
            G_struct = self.loadGroupStruct();
            gdf_in_group = G_struct.G(groupNumber).gdf;
            units_in_group = G_struct.G(groupNumber).units;
            
            units_idx_in_group = find(G_struct.G(groupNumber).units == number_within_group);
            unit_id_in_group = G_struct.G(groupNumber).units(units_idx_in_group);
            gdf_in_group = G_struct.G(groupNumber).gdf;
            
            % Get the spike indices for the unit:
            spike_idx = gdf_in_group(:,1) == unit_id_in_group;
            
            % Get the spike waveforms:
            spikes_struct = self.loadFileInGroup(groupNumber, 'spikes_cut');
            try
                wfs = mysort.h5.matrix(spikes_struct.spikeCut.wfs.fname, '/wfs', 1)
            catch
                warning('Could not read waveforms!')
                wfs = [];
            end
            
            % Get the electrode positions:
            G = self.loadGroupFile();
            electrodePositions = G.electrodePositions(G.groupsidx{groupNumber},:);
            channelNumbers = G_struct.G(groupNumber).sortedElectrodes;
            
        end
        
        
        %% PLOT FUNCTIONS #################################################
        function P = plotElectrodeGroups(self, varargin)
            P.fh = [];
            P.ah = [];
            P = util.parseInputs(P, varargin, 'error');
            
            S = self.loadGroupFile();
            if isempty(P.ah)
                if isempty(P.fh)
                    P.fh = mysort.plot.figure;
                end
                P.ah = axes;
            end
            set(P.ah, 'nextplot', 'add');
            N = self.getNElGroups();
            for i=1:N
                c = mysort.plot.vectorColor(i);
                m = mysort.plot.markerTypes(i);
                ep = S.electrodePositions(S.groupsidx{i},:);
                plot(P.ah, ep(:,1), ep(:,2), '.', 'marker', m, 'color', c, 'linewidth', 2, 'markersize', 16);
            end
        end
        
        
        %% DEBUGGING FUNCTIONS ############################################
        function runPostProcessing(self)
            %% Process all local sortings into a final sorting
            G = self.loadGroupFile();
            %             load(groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
            
            disp('Postprocessing...');
            [gdf_merged, T_merged, localSorting, localSortingID, NeuronCombinations] =...
                mysort.HDSorting.processLocalSortings(self.dPath, self.name, G.groups, G.groupsidx, ...
                'newPostProcFunc', 'newMergeFelix');
            units = unique(gdf_merged(:,1));
            nU = length(units);
            
            mysort.plot.figure([1200 900])
            ah = axes;
            set(ah, 'nextplot', 'add');
            nG = length(NeuronCombinations);
            for g=1:nG
                nGU = size(NeuronCombinations(g).templates.wfs,3);
                maxCol = 'r';
                notMaxCol = 'k';
                for gu=1:nGU
                    myp = id2pos(g,gu);
                    if NeuronCombinations(g).templates.maxInThisGroup(gu)
                        c = maxCol;
                    else
                        c = notMaxCol;
                    end
                    plot(ah, myp(1), myp(2), 'x', 'color', c, 'linewidth', 2, 'markersize', 16);
                    % Plot my masters
                    for m=1:size(NeuronCombinations(g).templates.masterTemplate{gu},1)
                        masterg  = NeuronCombinations(g).templates.masterTemplate{gu}(m,1);
                        mastergu =  NeuronCombinations(g).templates.masterTemplate{gu}(m,2);
                        masterp = id2pos(masterg, mastergu);
                        if g==1
                            plot(ah, [myp(1) masterp(1)], [myp(2) masterp(2)], '-', 'color', 'r', 'linewidth', 2);
                        else
                            plot(ah, [myp(1) masterp(1)], [myp(2) masterp(2)], '-', 'color', 'b', 'linewidth', 2);
                        end
                    end
                end
            end
            
            assert(length(R.localSorting) == nU, 'must be identical');
            assert(length(R.localSortingID) == nU, 'must be identical');
            assert(size(T_merged,3) == nU, 'must be identical');
            %             save(fullfile(outPath, [runName '_results.mat']), 'gdf_merged', 'T_merged', 'localSorting', 'localSortingID');
            
            function p = id2pos(g,gu)
                p = [g gu];
            end
        end
    end
    
    %     %% STATIC FUNCTIONS ###################################################
    %     methods(Static)
    %
    %         function submit_daemon()
    %             error('function deprecated!')
    %             cd('~/trunk/matlab')
    %             pd = pdefs();
    %             tokenFolder = fullfile( pd.tokenFiles );
    %             log_file = '~/submit_demon.log';
    %             cd(tokenFolder)
    %
    %             while 1
    %                 fnames = dir(fullfile(tokenFolder, 'start_*.mat'));
    %                 if ~isempty(fnames)
    %                     disp('Found token, processing...')
    %                     for i=1:length(fnames)
    %                         token_file = fullfile(tokenFolder, fnames(i).name);
    %                         t = load(token_file);
    %                         submit_token_file = fullfile(tokenFolder, 'submitted', fnames(i).name);
    %
    %                         disp(token_file)
    %
    %                         p = regexp(t.sortingFolder, filesep, 'split')
    %                         while ~strcmp(p{1}, 'Mea1k')
    %                             p = {p{2:end}};
    %                         end
    %                         fullSortingFolder = fullfile(pd.mea1kIntermediate, p{2:end});
    %
    %                         disp(['Sorting Location: ' fullSortingFolder])
    %                         disp(['Sorting Name: ' t.sortingName])
    %
    %                         submit_str = sprintf('qsub %s/%s_job.sh', fullSortingFolder, t.sortingName);
    %                         move_str   = sprintf('mv %s %s', token_file, submit_token_file);
    %
    %                         cd('~')
    %                         mysort.util.logToFile(log_file, submit_str)
    %                         [status, result] = system(submit_str);
    %                         mysort.util.logToFile(log_file, result)
    %                         cd(tokenFolder)
    %
    %                         if status == 0
    %                             disp('Submit successful')
    %                             disp(result)
    %                         else
    %                             disp('Submit failed')
    %                         end
    %
    %                         mysort.util.logToFile(log_file, move_str)
    %                         [status, result] = system(move_str);
    %                         mysort.util.logToFile(log_file, result)
    %                         pause(2.5)
    %                         disp('Done processing. Waiting...')
    %                     end
    %                 end
    %                 pause(10)
    %             end
    %             disp('Done.')
    %         end
    %     end
    
end