function [fileList fullPathFileList] = findFilesInFolder(folder, extension)
if nargin < 2
    extension = '.h5';
end

fileList_ = dir(folder);
fileList = {};
fullPathFileList = {};

for j = 1:length(fileList_)
    [path_ name_ ext] = fileparts(fileList_(j).name);
    
    if strcmp(ext, '.h5')
        fileList = {fileList{:}, [name_ ext]};
        fullPathFileList = {fullPathFileList{:}, fullfile(folder, [name_ ext])};
    end
end
%                     t = data_list(j).name(16:23);
%                     t(t=='_') = [];
%                     t = str2double(t);
%
%                     sameFiles = dir(fullfile(self.folders.recordings, 'data', [name_ '*']));
%
%                     %% Find Stimulus files:
%                     mat_type = 'unknown'; mat_file = '';
%                     for jj = 1:length(sameFiles)
%                         [path_ name_ ext] = fileparts(sameFiles(jj).name);
%                         if strcmp(ext, '.h5')
%                             h5_file = sameFiles(jj).name;
%                         elseif strcmp(ext, '.mat')
%                             mat_file = sameFiles(jj).name;
%                             try
%                                 mat_type = strsplit(sameFiles(jj).name, '_');
%                             catch
%                                 mat_type = strread(sameFiles(jj).name,'%s','delimiter','_');
%                             end
%                             mat_type = mat_type{6};
%                         elseif strcmp(ext, '.cmd')
%                             cmd_file = sameFiles(jj).name;
%                         end
%                     end
%
%                     %% Create new session or add the files to previous session
%                     config_idx = find(t > ct, 1, 'last');
%                     if config_idx > previous_config_idx;
%                         s = s + 1;
%                         sessions(s).configName = cs(config_idx).name(1:8);
%                         sessions(s).fileList = {data_list(j).name};
%
%                         if ~isempty(mat_file)
%                             sessions(s).stimulusFiles = {{mat_type; mat_file; cmd_file}'};
%                         else
%                             sessions(s).stimulusFiles = {};
%                         end
%                     else
%                         sessions(s).fileList = {sessions(s).fileList{:}, data_list(j).name}';
%                         if ~isempty(mat_file)
%                             sessions(s).stimulusFiles = { sessions(s).stimulusFiles{:}, {mat_type ; mat_file; cmd_file}'}';
%                         end
%                     end
%                     previous_config_idx = config_idx;
%                 end

