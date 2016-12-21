classdef PostprocessJob < grid.GridJob
    properties (SetAccess=private)
    end
    
    properties
        postprocessJobP
        sortingName
        %groupFile
        %groupsFolder
        sortingResult
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = PostprocessJob(sortingName, rootFolder, groupFile, groupsFolder, sortingResult, varargin)
            self = self@grid.GridJob(['postprocess_' sortingName], rootFolder, varargin{:});
            self.taskType = 'PostprocessJob';
            
            p = struct;
            p.newPostProcFunc = '';
            p = hdsort.util.parseInputs(p, self.P_untreated, 'error');
            self.postprocessJobP = p;
            
            if ~isempty(self.postprocessJobP.newPostProcFunc)
                error('newPostProcFunc not implemented!')
            end
            
            self.sortingName  = ['sort_' sortingName];
            self.files.group    = groupFile;
            self.folders.groups = groupsFolder;
            self.sortingResult = sortingResult;
            
            
            self.nWorkersPerNode = 16;
            self.startIndex = 1;
            self.endIndex = 1;
            self.taskIDs = self.startIndex:self.endIndex;
        end
        
        %------------------------------------------------------------------
        function setTaskParameters(self)
            
            taskParameters = struct;
            
            
            taskParameters.runName = self.jobName;
            taskParameters.sortingName  = self.sortingName;
            taskParameters.groupFile    = self.files.group;
            taskParameters.groupsFolder = self.folders.groups;
            taskParameters.sortingResult = self.sortingResult;
            
            %%
            if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                warning('Sorting not started from a linux machine might cause problems!')
                taskParameters.groupFile    = grid.GridJob.convertToLinux(self.files.group);
                taskParameters.groupsFolder = grid.GridJob.convertToLinux(self.folders.groups);
                taskParameters.sortingResult = grid.GridJob.convertToLinux(self.sortingResult);
            end
            
            disp('Create task files from groups...');
            %% Create cell variable allTaskParamters:
            for ii = 1:length(self.taskIDs)
                taskParameters.taskID = self.taskIDs(ii);
                self.allTaskParameters{ii} = taskParameters;
            end
           
            self.prepareTasks();
        end
        
    end
    
    methods(Static)
        %------------------------------------------------------------------
        %------------------------- RUN  FUNCTION ---------------------------
        %------------------------------------------------------------------
        function run(taskFile, debugFlag)
            if nargin < 2
                debugFlag = false;
            end
            
            %% Load taskFile:
            T = load(taskFile);
            taskP = T.taskParameters
            clear T;
           
            %% Check necessary parameters:
            assert( isfield(taskP, 'groupFile'), 'Task aborted: field taskParameters.rawFile not specified!');
            
            if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                taskP.groupFile = grid.GridJob.convertToLocal(taskP.groupFile);
                taskP.groupsFolder = grid.GridJob.convertToLocal(taskP.groupsFolder);
                taskP.reportFile = grid.GridJob.convertToLocal(taskP.reportFile);
            end
            
            assert( exist(taskP.groupFile, 'file') == 2, ['Task aborted: groupFile ' taskP.groupFile ' not found!']);
            assert( exist(taskP.groupsFolder, 'dir') == 7, ['Task aborted: no valid groupsFolder specified! Path given:' taskP.groupsFolder]);
            assert( isfield(taskP, 'taskID'), 'Task aborted: field taskParameters.taskID not specified!');
            
            %% (Re-)Set reporting file:
            rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
            rep(:,:) = [0 0];
         
            if ~debugFlag
                try
                    mainBlock();
                catch ME
                    errorHandling(ME);
                end
            else
                mainBlock();
            end
            
            function mainBlock()
                disp('Load group information...');
                GF = load(taskP.groupFile, 'groups', 'electrodeNumbers', 'electrodePositions', 'nGroupsPerElectrode', 'groupsidx');
                
                disp('Start postprocessing...');
                [R, P] = mysort.HDSorting.processLocalSortings(taskP.groupsFolder,...
                    taskP.sortingName, GF.groups, GF.groupsidx, ...
                    'groupPaths', taskP.groupsFolder);
                
                disp('Check and save data...');
                units = unique(R.gdf_merged(:,1));
                nU = length(units);
                assert(length(R.localSorting) == nU, 'must be identical');
                assert(length(R.localSortingID) == nU, 'must be identical');
                assert(size(R.T_merged,3) == nU, 'must be identical');
                
                disp('Saving postprocessing results...')
                save(taskP.sortingResult, '-struct', 'R', '-v7.3');
                
                %% Write to reporter file:
                disp('Writing results...')
                rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [1 0];
            end
            
            function errorHandling(ME)
                
                disp('Catch error...')
                errStr = mysort.hdsort.util.buildLastErrString(ME);
                disp(errStr)
                
                rep = mysort.ds.binaryFileMatrix(taskP.reportFile, [1 2], 'writable', true);
                rep(:,:) = [0 1];
                rethrow(ME)
            end
        end
        
    end
end
