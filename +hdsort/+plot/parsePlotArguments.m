function P = parsePlotArguments(P_in, varargin)

if mod(numel(varargin{:}), 2)
    varargin{:}
    varargin = {P_in, varargin{:}};
    varargin{:}
    P_in = struct();
end

P.ah = [];
P.fh = [];
P.title = '';
P.Fs = 20000;
P.axis = [];
P.ylabel = '';
P.xlabel = '';
P.xTicklabel = [];
P.yTicklabel = [];
P.color = [];
P.nColors = [];

P.subplot = [];
P.spacerX = 0.1;
P.spacerY = 0.1;
P.offsetX = 0.1;
P.offsetY = 0.1;

%P = mysort.util.parseInputs(P, varargin, 'error');

%% Check if there is already a default value:
for f = fields(P)'
    if ~isfield(P_in, f)
        P_in.(f{:}) = P.(f{:});
    end
end

[P untreatedArgs] = mysort.util.parseInputs(P_in, varargin{:}, 'error');

untreatedArgs = [fieldnames(untreatedArgs)'; struct2cell(untreatedArgs)']';

if isempty(P.ah) & isempty(P.fh)
    P.fh = mysort.plot.figure('name',P.title);
end
if ~isempty(P.subplot)
    P.ah = mysort.plot.subplots(P.subplot, 'spacerX', P.spacerX,'spacerY', P.spacerY, 'offsetX', P.offsetX, 'offsetY', P.offsetY);
elseif isempty(P.ah)
    P.ah = axes();
end
if numel(P.ah) == 1
    axes(P.ah);
else
    axes(P.ah(1));
end

%% Take care of the colors:
%if ~isnumeric(P.nColors)
%    assert(isempty(P.color), 'Do not specify nColors if color is not empty!')

if isempty(P.color)
    if ~isempty(P.nColors)
        %% P.color is left empty, but nColors is specified:
        P.color = mysort.plot.vectorColor(1:P.nColors);
    end
elseif isstr(P.color)
    %% P.color is a string specifying a colormap:
    colormap(P.color);
elseif (numel(P.color) == 1) & isnumeric(P.color)
    %% P.color is single numeric value:
    if P.color < 1
        P.color = zeros(P.nColors,3);
    else
        P.color = mysort.plot.vectorColor(1:P.color);
    end
end

