function newPath = convertPathToOS(path, OStype)
if isempty(path)
    newPath = [];
    return;
end

if nargin == 1
    OStype = computer;
elseif ~isempty(strfind(OStype, 'lin')) || ~isempty(strfind(OStype, 'LIN'))
    OStype = 'GLNXA64';
elseif ~isempty(strfind(OStype, 'win')) || ~isempty(strfind(OStype, 'WIN'))
    error('Not supported yet!')
    OStype = 'PCWIN64';
elseif ~isempty(strfind(OStype, 'mac')) || ~isempty(strfind(OStype, 'MAC'))
    OStype = 'MACI64';
end

if iscell(path)
    newPath = {};
    for ii = 1:numel(path)
        newPath{ii} = hdsort.util.convertPathToOS(path{ii}, OStype);
    end
    return
end

if isempty(path) 
    newPath = '';
    return
end
if strcmp(path(1), '.')
    newPath = path;
    return
end

%% Check the path until it finds the correct match
pd_mac = hdsort.pathDefinitions('MACI64');
pd_win = hdsort.pathDefinitions('PCWIN64');
pd_lin = hdsort.pathDefinitions('GLNXA64');

bestMatchRest = []; bestField = []; nBestMatch = 1e5;
for pd_ = {pd_mac, pd_lin, pd_win}
    fields_ = fields(pd_{1});
    
    for f_ = fields_'
        
        
        f = f_{1};
        
        if ~isempty( strfind(path, getfield(pd_{1}, f)) )
            rest = strsplit(path, getfield(pd_{1}, f));
            
            N = length(rest{2});
            
            if N < nBestMatch
                nBestMatch = N;
                bestMatchRest = rest{2};
                bestField = f;
            end
        end
    end
end
assert(~isempty(bestField), 'Path not found, can not be converted!')

pd = hdsort.pathDefinitions(OStype);
newPath = fullfile( getfield(pd, bestField), bestMatchRest);

end