classdef PreprocessJob < grid.GridJob
    properties (SetAccess=private)
    end
    
    properties
        preprocessJobP
        
        oldMea1kFiles
        rawFiles
        
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = PreprocessJob(sortingName, rootFolder, rawFiles, varargin)
            self = self@grid.GridJob(['preprocess_' sortingName], rootFolder, varargin{:});
            self.taskType = 'PreprocessJob';
            
            p = struct;
            p.configFile = ''; % Legacy for old mea1k-files
            p.destinationFolder = '';
            
            p = hdsort.util.parseInputs(p, self.P_untreated, 'error');
            self.preprocessJobP = p;
            
            %% Set the data files:
            if ~iscell(rawFiles) rawFiles = {rawFiles}; end
            self.rawFiles = rawFiles;
            
            if ~isempty(self.preprocessJobP.destinationFolder)
                self.folders.destination = self.preprocessJobP.destinationFolder;
            else
%                self.folders.destination = self.folders.main;
                self.folders.destination = self.folders.root;
            end
            
            
            if ~isempty(self.preprocessJobP.configFile)
                self.oldMea1kFiles = true;
            else
                self.oldMea1kFiles = false;
            end
            
            self.startIndex = 1;
            self.endIndex = numel(self.rawFiles);
            self.taskIDs = self.startIndex:self.endIndex;
        end
        
        function bool = alreadyPreprocessed(self)
             [preprocessedFiles, cmdFiles, rawFiles] = self.getFileNames();
             for preprocessedFile = preprocessedFiles(:)'
                 try
                     DS = mysort.mea.CMOSMEA(preprocessedFile{1});
                 catch
                     bool = false;
                    return;
                 end
             end
             bool = true;
        end
        
        %------------------------------------------------------------------
        % This is also a legacy function for old mea1k-file format:
        function [preprocessedFiles, rawFiles, cmdFiles] = getFileNames(self)
            rawFiles = self.rawFiles;
            preprocessedFiles = cell(length(rawFiles), 1);
            cmdFiles = cell(length(rawFiles), 1);
            
            for k = 1:length(rawFiles)
                rawFile = rawFiles{k};
                rootName = getRootName(rawFile);
                
                preprocessedFiles{k} = fullfile(self.folders.destination, [rootName '.h5']);
                cmdFiles{k} = fullfile(fileparts(rawFile), [rootName '.raw.cmd']);
            end
            function r = getRootName(fileName)
                [pathstr,r,ext] = fileparts(fullfile( fileName ));
                if ~strcmp(ext, '')
                    r = getRootName(r);
                end
            end
        end
        
        %------------------------------------------------------------------
        function setTaskParameters(self)
            
            taskParameters = struct;
            taskParameters.runName = self.jobName;
            taskParameters.destinationFolder = self.folders.destination;
            
            if self.oldMea1kFiles
                taskParameters.configFile = self.preprocessJobP.configFile;
                taskParameters.oldMea1kFiles = true;
            end
            
            if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                warning('Sorting not started from a linux machine might cause problems!')
                taskParameters.destinationFolder = grid.GridJob.convertToLinux(taskParameters.destinationFolder);
                taskParameters.configFile = grid.GridJob.convertToLinux(taskParameters.configFile);             
            end
            
            disp('Create task files from groups...');
            %% Create cell variable allTaskParamters:
            for ii = 1:length(self.taskIDs)
                taskParameters.taskID = self.taskIDs(ii);
                
                if ~isempty(strfind(computer, 'WIN')) | ~isempty(strfind(computer, 'MACI64'))
                	taskParameters.rawFile = grid.GridJob.convertToLinux(self.rawFiles{ii});
                else
                    taskParameters.rawFile = self.rawFiles{ii};
                end
                
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
            assert( isfield(taskP, 'rawFile'), 'Task aborted: field taskParameters.rawFile not specified!');
            assert( exist(taskP.rawFile, 'file') == 2, ['Task aborted: dataFile ' taskP.rawFile ' not found!']);
            assert( exist(taskP.destinationFolder, 'dir') == 7, ['Task aborted: no valid destinationFolder specified! Path given:' taskP.destinationFolder]);
            %assert( isfield(taskP, 'groupidx'), 'Task aborted: field taskParameters.groupidx not specified!');
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
                if isfield(taskP, 'oldMea1kFiles')
                    preprocessedFiles = mysort.mea.preprocessAllMea1kH5Files(taskP.rawFile, ...
                        taskP.configFile, taskP.destinationFolder)
                else
                    fileDS = mysort.mea.MeaRecording(taskP.rawFile);
                    preprocessedFiles = fileDS.preprocessFile(taskP.destinationFolder);
                end
                
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
