function c = config()
% This software supports only grid engines that are run on linux operating 
% systems. In order to be able to start a sorting from Windows or Mac 
% computers, you must specify the paths relative to these systems.
% You must further specify the file location path on the linux system and a
% folder name that is common to all operating systems. The GridJob object 
% will then create the path names in the correct linux syntax independent 
% form the computer on which the sorting was started.

user_name = getenv('USER');

if ~isempty(strfind(computer, 'WIN')) 
    %% Windows
    error('Not specified yet!')
    
    c.home = '';
    c.localSortingPath = '';
    c.tokenFilesFolder = 'S:\group\hierlemann_Mea1k\intermediate_data\Mea1k\tokens\';
    
elseif ~isempty(strfind(computer, 'MACI64')) 
    %% Mac
    c.home = fullfile('/', 'Volumes', 'hierlemann');
    c.localSortingPath = fullfile('/', 'Volumes', 'hierlemann');
    c.tokenFilesFolder = fullfile(c.localSortingPath, 'intermediate_data', 'Mea1k', user_name, 'tokens');
    
else
    %% Linux
    c.home = fullfile('/', 'net', 'bs-filesvr02', 'export', 'group', 'hierlemann');
    c.localSortingPath = fullfile('/', 'net', 'bs-filesvr02', 'export', 'group', 'hierlemann');
    c.tokenFilesFolder = fullfile(c.localSortingPath, 'intermediate_data', 'Mea1k', user_name, 'tokens');
end

c.linuxSortingPath = '/net/bs-filesvr02/export/group/hierlemann';
c.commonFolderName = 'Mea1k';  

%c.submit_host = 'bs-submit01';

end


