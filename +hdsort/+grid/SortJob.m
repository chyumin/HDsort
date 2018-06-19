%Example:
% sj = grid.SortJob(jobName, dataFiles, destFolder);
% sj.setTaskParameters();
% sj.prepareTasks();
%
% --> Submit to grid
%
% all_tasks_completed = sj.waitForTasksToFinish();
% if all_tasks_completed
%     sj.copyBackResults();
%     sj.postprocessBOTMSorting();
% end

classdef SortJob < hdsort.grid.GridJob
    properties (SetAccess=private)
    end
    
    properties
        name
        LEGs
        sortingParameters
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = SortJob(name, rootFolder, groupsFolder, LEGs, preprocessedFiles, sortingParameters, varargin)
            
            self = self@hdsort.grid.GridJob(['sort_' name], rootFolder, varargin{:});
            self.taskType = 'SortJob';
            
            self.name = name;
            self.LEGs = LEGs;
            self.sortingParameters = sortingParameters;
            self.folders.groups = groupsFolder;
            
            if ~iscell(preprocessedFiles) preprocessedFiles = {preprocessedFiles}; end
            self.files.preprocessed = preprocessedFiles;
            
            self.startIndex = 1;
            self.endIndex = self.LEGs.N;
            self.taskIDs = self.startIndex:self.endIndex;
        end
        
        function setTaskParameters(self)
            
            taskParameters.sortingParameters = self.sortingParameters;
            taskParameters.name = self.name;
            taskParameters.preprocessedFileList = self.files.preprocessed;
            assert(~isempty(taskParameters.preprocessedFileList), 'You must provide a list of preprocessed files!')
            
            
            disp('Create task files from groups...');
            % Create cell variable allTaskParamters:
            for ii = 1:length(self.taskIDs)
                
                taskParameters.groupFolder = hdsort.util.convertPathToOS(fullfile(self.folders.groups, self.LEGs.name{ii}));
                assert(exist(taskParameters.groupFolder) == 7, 'groupFolder does not exist!')
                
                taskParameters.wfsFile = fullfile(taskParameters.groupFolder, [taskParameters.name '.040spikes_cut.h5']);
                assert(exist(taskParameters.wfsFile) == 2, 'wfsFile does not exist!')
                
                taskParameters.covFile = fullfile(taskParameters.groupFolder, [taskParameters.name '.060cov.mat']);
                assert(exist(taskParameters.covFile) == 2, 'covFile does not exist!')
                
                taskParameters.groupidx = self.LEGs.groupsidx{ii};
                taskParameters.taskID = self.taskIDs(ii);
                
                self.allTaskParameters{ii} = taskParameters;
            end
            
            self.prepareTasks();
        end
    end
    
    %------------------------------------------------------------------
    methods(Static)
        %------------------------------------------------------------------
        %------------------------- RUN  FUNCTION ---------------------------
        %------------------------------------------------------------------
        function run(taskFile, debugFlag)
            if nargin < 2
                debugFlag = false;
            end
            
            %% Load taskParameters:
            T = load(taskFile);
            taskParameters = T.taskParameters;
            
            taskParameters.groupFolder = hdsort.util.convertPathToOS(taskParameters.groupFolder);
            taskParameters.wfsFile = hdsort.util.convertPathToOS(taskParameters.wfsFile);
            taskParameters.covFile = hdsort.util.convertPathToOS(taskParameters.covFile);
            taskParameters.reportFile = hdsort.util.convertPathToOS(taskParameters.reportFile);
            taskParameters.preprocessedFileList = hdsort.util.convertPathToOS(taskParameters.preprocessedFileList);
            
            %% Check that all files exist:
            assert(exist(taskParameters.wfsFile) == 2, 'wfs file does not exist!')
            assert(exist(taskParameters.covFile) == 2, 'cov file does not exist!')
            assert(exist(taskParameters.groupFolder) == 7, 'groupFolder does not exist!')
            
            reportFolder = fileparts(taskParameters.reportFile);
            assert(exist(reportFolder) == 7, 'reportFolder does not exist!')
            
            assert( isfield(taskParameters, 'sortingParameters'), 'Task aborted: field taskParameters.sortingParameters does not exist!')
            
            %% (Re-)Set reporting file:
            rep = hdsort.filewrapper.util.BinaryFileMatrix(taskParameters.reportFile, [1 2], 'writable', true);
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
                
                hdsort.leg.sortOneLEG(taskParameters.name, ...
                                            taskParameters.wfsFile, ...
                                            taskParameters.covFile, ...
                                            taskParameters.groupFolder, ...
                                            taskParameters.sortingParameters, ...
                                            taskParameters.groupidx, ...
                                            taskParameters.preprocessedFileList);
                
                %% Write to reporter file:
                rep = hdsort.filewrapper.util.BinaryFileMatrix(taskParameters.reportFile, [1 2], 'writable', true);
                rep(:,:) = [1 0];
            end
            
            function errorHandling(ME)
                
                disp('Catch error...')
                errStr = hdsort.util.buildLastErrString(ME);
                disp(errStr)
                
                rep = hdsort.filewrapper.util.BinaryFileMatrix(taskParameters.reportFile, [1 2], 'writable', true);
                rep(:,:) = [0 1];
                rethrow(ME)
            end
        end
        
    end
end
