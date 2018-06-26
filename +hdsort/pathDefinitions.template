function def = pathDefinitions(varargin)
    if nargin == 0
        OStype = computer;
    else
        OStype = varargin{1};
    end

    def.pdefsPath = pwd;
    def.callerPath = [];
    
    user_config = hdsort.global_user_config();
    user_name = user_config.user_name;
    
    if ~isempty(strfind(OStype, 'WIN'))
        
        def.server01 = '\\bs-filesvr01\';
        def.server02 = '\\bs-filesvr02\';
        
    elseif ~isempty(strfind(OStype, 'MACI64'))
        
        def.server01 = fullfile('/Volumes', 'filesvr01'); 
        def.server01 = fullfile('/Volumes', 'filesvr02');

        %if ~exist(def.server01) warning('/Volumes/filesvr01 does not exist!'); end
        %if ~exist(def.server01) warning('/Volumes/filesvr02 does not exist!'); end
        
    elseif ~isempty(strfind(OStype, 'GLNXA64'))
        
        def.server01 = fullfile('/net', 'bs-filesvr01', 'export', 'group', 'groupname');
        def.server02 = fullfile('/net', 'bs-filesvr02', 'export', 'group', 'groupname');
    
    else
        error('OS unknown!')
    end
    
    def.recordings = fullfile(def.server01, 'recordings');
    def.local = '.';
end
