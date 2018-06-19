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

fields_ = [fields(pd_mac); fields(pd_lin); fields(pd_win)];
fields_ = unique(fields_);
for f_ = fields_'
    f = f_{1};
    
    try
        if ~isempty(strfind(path, getfield(pd_lin, f)))
            rest = strsplit(path, getfield(pd_lin, f));
            delim = '/';
            break;
        end
    catch
    end
    
    try
        if ~isempty(strfind(path, getfield(pd_mac, f)))
            rest = strsplit(path, getfield(pd_mac, f));
            delim = '/';
            break;
        end
    catch
    end
    
    try
        if ~isempty(strfind(path, getfield(pd_win, f)))
            rest = strsplit(path, getfield(pd_win, f));
            delim = '\';
            break;
        end
    catch
    end
    f = '';
end
assert(~isempty(f), 'Path not found, can not be converted!')

rest = strsplit(rest{end}, delim);

%%
pd = hdsort.pathDefinitions(OStype);
newPath = fullfile( getfield(pd, f), rest{:});

end