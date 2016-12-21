classdef PlotInterface < handle
    properties (SetAccess=protected)
        titleBox
        
        hdsort.plot.bj
        
        Children
        Parent
    end
    
    properties
        hdsort.plot.ame
        
        ah
        fh
        
        hold_on
        title
        Fs
        axis
        ylabel
        xlabel
        legend
        
        XTick
        YTick
        XTickLabel
        YTickLabel
        XTickMode
        YTickMode
        XTickLabelMode
        YTickLabelMode
        XScale
        YScale
        
        xlim
        ylim
        
        color
        nColors
        
        Transparency
        MarkerFaceAlpha
        MarkerEdgeAlpha
        FaceAlpha
        EdgeAlpha
        
        FontSize
        FontWeight
        LineWidth
        LineSpec
        
        showAh
        %doNotSetAh
        showOnStartup
        
        width
        height
        
        addTo
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = PlotInterface(varargin)
            
            % Treate various combinations of input arguments:
            P_in = struct();
            %parseargs = false;
            if nargin == 1
                assert(isstruct(varargin{1}), 'If there is only one input, it must be a struct!')
                P_in = varargin{1};
                varargin = {};
            elseif nargin > 1
                %   parseargs = true;
                if isstruct(varargin{1})
                    assert(mod(nargin, 2) == 1, 'Input arguments must come as pairs!')
                    P_in = varargin{1};
                    varargin = {varargin{2:end}};
                else
                    assert( ischar(varargin{1}), 'First input arguments must be a string!')
                    assert( mod(nargin, 2) == 0, 'Input arguments must come as pairs!')
                    P_in = struct();
                end
            end
            
            P.addTo = [];
            P.ah = [];
            P.fh = [];
            P.hold_on = true;
            P.title = '';
            P.Fs = 20000;
            P.axis = [];
            P.ylabel = '';
            P.xlabel = '';
            P.legend = [];
            
            P.XTick = [];
            P.YTick = [];
            P.XTickLabel = [];
            P.YTickLabel = [];
            P.XTickMode = 'auto';
            P.YTickMode = 'auto';
            P.XTickLabelMode = 'auto';
            P.YTickLabelMode = 'auto';
            P.XScale = 'linear';
            P.YScale = 'linear';
            
            P.xlim = [];
            P.ylim = [];
            
            P.color = [];
            P.nColors = [];
   
            P.Transparency = 1.0;
            P.MarkerFaceAlpha = 1.0;
            P.MarkerEdgeAlpha = 1.0;
            P.FaceAlpha = 1.0;
            P.EdgeAlpha = 1.0;
            
            P.FontSize = 10;
            P.FontWeight = 'normal';
            P.LineWidth = 0.5;
            P.LineSpec = '-';
            
            P.showAh = true;
            P.showOnStartup = true;
            
            %% Add a default value to P_in if it is not already set in subclass:
            for f = fields(P)'
                if ~isfield(P_in, f)
                    P_in.(f{:}) = P.(f{:});
                end
            end
            
            P = hdsort.util.parseInputs(P_in, varargin, 'error');
            
            self.setProperties(P, varargin);
            
            if isempty(self.title) self.title = ''; end
            
            
            self.Children = {};
            self.Parent = [];
            
            self.setColor(self.color, self.nColors);
            
        end
        
        
        % -----------------------------------------------------------------
        function  setProperties(self, P, overrideSettings)
            
            if iscell(P)
                P_ = P;
                P = struct;
                assert(mod(numel(P_), 2) == 0, 'Cell must be pairs of field values.')
                for ii = 1:2:numel(P_)
                    assert(ischar(P_{ii}), 'Cell must be pairs of field values.')
                    P.(P_{ii}) = P_{ii+1};
                end
            end
            
            for f_ = fields(P)'
                f = f_{1};
                self.(f) = P.(f);
            end
            
            %% Copy Properties of other PlotInterface:
            if ~isempty(P.addTo)
                for f_ = properties(hdsort.plot.PlotInterface)'
                    f = f_{1};
                    self.(f) = P.addTo.(f);
                end
                
                for ii = 1:2:numel(overrideSettings)
                    self.(overrideSettings{ii}) = overrideSettings{ii+1};
                end
            end
            
        end
        
        % -----------------------------------------------------------------
        function fh_ = createFigure(self, varargin)
            P.color = 'w';
            P.position = [];
            P.width = 600;
            P.height = 300;
            P.w = [];
            P.h = [];
            nargin_ = nargin;
            if nargin_ > 0 && isnumeric(varargin{1}) && length(varargin{1}) == 2
                % special case where it is called as figure([w h], ...)
                P.w = varargin{1}(1);
                P.h = varargin{1}(2);
                varargin = varargin(2:end);
                nargin_ = nargin_ -1;
            end
            
            if nargin_==1
                fh = varargin{1};
                figure(fh);
            else
                [P restargs] = hdsort.util.parseInputs(P, varargin, 'split');
                restargs = hdsort.util.deflateP(restargs);
                fh = figure(restargs{:});
            end
            if ~isempty(P.h); P.height = P.h; end
            if ~isempty(P.w); P.width  = P.w; end
            set(fh,'color', P.color);
            if isempty(P.position)
                screen_size = get(0, 'ScreenSize');
                position = get(fh, 'position');
                position(1) = round((screen_size(3) - P.width)/2);
                position(3) = P.width;
                position(2) = round((screen_size(4) - P.height)/2);
                position(4) = P.height;
                set(fh, 'position', position);
            end
            self.dataCursorMode(fh);
            if nargout == 1
                fh_ = fh;
            end
        end
        
        % -----------------------------------------------------------------
        function setFigureTitle(self, title_str)
            % hdsort.plot. a title to a figure not to a single axes
            % inputs:
            %   title - title string if figure handle is provided
            self.setAh(self.ah);
            
            if nargin > 1
                self.title = title_str;
            end
            
            units = get(self.fh, 'Units');
            set(self.fh, 'Units', 'pixel');
            fig_pos_pixel = get(self.fh, 'Position');
            set(self.fh, 'Units', units);
            
            delete(self.titleBox);
            
            self.titleBox = annotation(self.fh, ...
                'textbox', [0.1 .1 .1 .1],...
                'Units', 'pixel',...
                'String',{self.title},...
                'HorizontalAlignment','center',...
                'Interpreter', 'none',...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'fontsize', self.FontSize, ...
                'FontWeight', self.FontWeight);
            set(self.titleBox, 'position', [0 fig_pos_pixel(4)-40 fig_pos_pixel(3) 40]);
            set(self.titleBox, 'Units', 'normalized');
        end
        
        % -----------------------------------------------------------------
        function setAxisTitle(self, title_str)
            % hdsort.plot. a title to a figure not to a single axes
            % inputs:
            %   title - title string if figure handle is provided
            self.setAh(self.ah);
            if nargin > 1
                self.title = title_str;
            end
            title(self.ah, self.title, ...
                'FontSize', self.FontSize, ...
                'FontWeight', self.FontWeight);
        end
        
        % -----------------------------------------------------------------
        function setAxis(self, varargin)
            self.setAh(self.ah);
            axis0 = [self.ah.XLim self.ah.YLim];
            axis_in = axis0;
            
            if nargin == 2
                axis_in = varargin{1};
                assert(numel(axis_in) == 4, 'Axis must be a 4-element vector!')
            elseif nargin > 2
                field = varargin{1};
                val = [varargin{2:end}];
                assert(numel(val) == 2, 'Input must be a 2 element vector!')
                
                self.setAh(self.ah);
                if strcmp(field, 'x')
                    %val(isnan(val)) = self.ah.XLim(isnan(val));
                    axis_in = [val self.ah.YLim];
                elseif strcmp(field, 'y')
                    %val(isnan(val)) = self.ah.YLim(isnan(val)); 
                    axis_in = [self.ah.XLim val];
                else
                    error('Input argument must be either x or y')
                end
            end
            
            axis_in(isnan(axis_in)) = axis0(isnan(axis_in));
            self.axis = axis_in;
            %self.show();
        end
        
        % -----------------------------------------------------------------
        function ah = structAh(self, ah_in)
            if ~isempty(ah_in) & isnumeric(ah_in)
                axes(ah_in);
                ah = gca;
            else
                ah = ah_in;
            end
        end
        
        % -----------------------------------------------------------------
        function setAh(self, ah)
            
            ah = self.structAh(ah);
            
            if isempty(ah) || ~isvalid(ah)
                self.setFh(self.fh);
            else
                figure(ah.Parent);
                axes(ah);
            end
            self.fh = gcf;
            self.ah = gca;
            
            if ~self.showAh
                hidem(self.ah);
            end
        end
        
        % -----------------------------------------------------------------
        function setFh(self, fh)
            if isempty(fh) || ~isvalid(fh)
                self.fh = self.createFigure('name', self.title);
                self.ah = axes();
            else
                figure(fh);
                try
                    self.setAh(fh.Children(1));
                catch
                    self.ah = axes();
                end
            end
            self.fh = gcf;
            self.ah = gca;
        end
        
        % -----------------------------------------------------------------
        function setFontSize(self, fontsize)
            if nargin > 1
                self.FontSize = fontsize;
            end
            
            %axesHandle = findall(figureHandle,'type','axes');
            self.setAh(self.ah);
            handles = self.ah;
            handles = [handles; findall(self.ah,'Type','text')];
            set(handles, 'FontSize', self.FontSize);
        end
        
        % -----------------------------------------------------------------
        function setLineWidth(self, linewidth)
            if nargin > 1
                self.LineWidth = linewidth;
            end
            
            self.setAh(self.ah);
            handles = self.ah;
            handles = [handles; findall(self.ah, 'Type', 'Line')];
            set(handles, 'LineWidth', self.LineWidth);
        end
        
        % -----------------------------------------------------------------
        function showLegend(self, varargin)
            if nargin < 2
                legend_labels = self.legend;
            elseif nargin == 2
                legend_labels = varargin{1};
            else
                legend_labels = varargin;
            end
            
            if ~isempty(legend_labels)
                legend(legend_labels);
            end
        end
        
        % -----------------------------------------------------------------
        function setColor(self, color, nColors)
            if isempty(color)
                if ~isempty(nColors)
                    %% P.color is left empty, but nColors is specified:
                    self.color = self.vectorColor(1:nColors);
                end
            elseif ischar(color)
                %% color is a string specifying a colormap:
                self.color = colormap(color);
            elseif (numel(color) == 1) && isnumeric(color)
                %% color is single numeric value:
                if color < 1
                    self.color = zeros(nColors,3);
                else
                    self.color = self.vectorColor(1:color);
                end
            elseif numel(color) == 3 && nargin == 3 && ~isempty(nColors)
                self.color = repmat(color(:)', nColors, 1);
            elseif size(color, 2) == 3
                self.color = color;
            end
        end
        
        function useDefaultValues(self)
            self.XTick = get(self.ah, 'XTick');
            self.XTickLabel = get(self.ah, 'XTickLabel');
            self.YTick = get(self.ah, 'YTick');
            self.YTickLabel = get(self.ah, 'YTickLabel');
        end
        
        % -----------------------------------------------------------------
        function update(self, varargin)
            
            if nargin > 1
                assert(mod(nargin, 2)== 1, 'Arguments must come in pairs!')
                
                for v = 1:2:numel(varargin)%{1:2:end}
                    try
                        self.(varargin{v}) = varargin{v+1};
                    catch
                        warning(['Property ' varargin{v} ' could not be updated!']);
                    end
                end
                
            end
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            error('This function must be implemented in the derived class!')
        end
        
        % -----------------------------------------------------------------
        function show(self)
            if ~self.showOnStartup
                self.showOnStartup = true;
                return;
            end
            
            % Got to the correct axis handle or create a new one:
            self.setAh(self.ah);
            %assert(self.ah.Parent == self.fh, 'Ah must be a child of fh!')
            
            if self.hold_on
                hold on;
            end
            
            self.show_();
            
            if ~isa(self.ah, 'matlab.graphics.axis.PolarAxes')
                xlabel(self.xlabel);
                ylabel(self.ylabel);
                
                if isvalid(self.ah)
                    
                    set(self.ah, 'XTick', self.XTick);
                    set(self.ah, 'XTickLabel', self.XTickLabel);
                    
                    set(self.ah, 'YTick', self.YTick);
                    set(self.ah, 'YTickLabel', self.YTickLabel);
                    
                    set(self.ah, 'XTickMode', self.XTickMode);
                    set(self.ah, 'YTickMode', self.YTickMode);
                    
                    set(self.ah, 'XTickLabelMode', self.XTickLabelMode);
                    set(self.ah, 'YTickLabelMode', self.YTickLabelMode);
                    
                    set(self.ah, 'XScale', self.XScale);
                    set(self.ah, 'YScale', self.YScale);
                    
                end
                
                if ~isempty(self.axis) axis(self.axis); end
                if ~isempty(self.xlim) xlim(self.xlim); end
                if ~isempty(self.ylim) ylim(self.ylim); end
            end
            
            self.setAxisTitle();
            self.setTransparency();
        end
        
        % -----------------------------------------------------------------
        function resizeFigure(self, width, height)
            
            if ~isvalid(self.fh)
                self.show();
            end
            set(self.fh, 'Units', 'pixels');
            
            if nargin == 2
                assert(numel(width) == 2, 'Specify width and height as a 2 element vector!')
                self.width = width(1);
                self.height = width(2);
            elseif nargin == 3
                self.width = width;
                self.height = height;
            else
                error('Not implemented yet!')
                
                
%                 
%                 maxS = self.fh.Position(3:4);
%                 %maxSY = self.fh.Position(4);
%                 for ii = 1:numel(self.fh.Children)
%                     try
%                         ch = self.fh.Children(ii);
%                         set(ch, 'Units', 'pixels');
%                         %p = self.fh.Children(ii).Position(1:2) + self.fh.Children(ii).Position(3:4);
%                         %self.fh.Children(ii)
%                         %self.fh.Children(ii).Position
%                         [min_ max_] = self.getMaxLabelPosition(ch)
%                         p = max_ - min_;
%                         p = p(1:2);
%                     catch
%                         p = [0,0];
%                     end
%                     maxS
%                     p
%                     maxS = max(maxS, p);
%                 end
%                 
%                 self.width = maxS(1);
%                 self.height = maxS(2);

                self.width = self.ah.OuterPosition(3);
                self.height = self.ah.OuterPosition(4);
                
            end
            
            set(self.fh, 'position', [self.fh.Position(1) self.fh.Position(2) self.width self.height]);%*11.0/12.0])
        end
        
        % -----------------------------------------------------------------
        function closed = closeFigure(self)
            closed = false;
            try
                closed = close(self.fh);
            catch
                closed = true;
            end
        end
        
        % -----------------------------------------------------------------
        function closed = closeAxis(self)
            % closed = false;
            try
                delete(self.ah);
                closed = true;
            catch
                closed = false;
            end
        end
        
        % -----------------------------------------------------------------
        function saveFig(self, fname, varargin)
            %From mysort.hdsort.plot.savefig(fig, fname):
            P.dpi = 300;
            P.fig = 1;
            P.png = 1;
            P.pdf = 0;
            P.eps = 1;
            P.emf = 0;
            P.ai = 0;
            P.special = 0;
            P.render_manual = true;
            P = hdsort.util.parseInputs(P, varargin, 'error');
            
            if nargin == 1
                fname = self.hdsort.plot.ame;
                warning(['No filename specified. Figure saved as ' fname '.fig'])
            end
            
            if isempty(self.fh) || ~isvalid(self.fh)
                self.show();
            end
            
            if P.render_manual
                self.fh.RendererMode = 'manual';
                self.fh.Renderer = 'painters';
            end
            
            
            set(self.fh, 'PaperPositionMode', 'auto');   % Use screen size
            visibility = get(self.fh, 'visible');
            set(self.fh, 'visible', 'on');
            
            % From mysort.hdsort.plot.figureChildrenSet(self.fh, 'Units', 'normalized'):
            count = 0;
            try
                set(self.fh, 'Units', 'normalized');
                count = count + 1;
            catch
                % Ignore this
            end
            c = get(self.fh, 'Children');
            for i=1:length(c)
                count = count + self.figureChildrenSet(c(i), 'Units', 'normalized');
            end
            
            % Save:
            if P.png
                print(self.fh, ['-r' num2str(P.dpi)], '-dpng', [fname '.png']);
            end
            if P.emf
                saveas(self.fh, [fname '.emf'], 'emf');
            end
            if P.fig
                saveas(self.fh, [fname '.fig'], 'fig');
            end
            if P.eps
                saveas(self.fh, [fname '.eps'], 'eps');
            end
            if P.ai
                saveas(self.fh, [fname '.ai'], 'ai');
            end
            if P.special || P.pdf
                print(self.fh, [fname '.pdf'], ['-d' 'pdfwrite'] )
            end
            set(self.fh, 'visible', visibility);
        end
        
        
        
        % -----------------------------------------------------------------
        function setTransparency(self, varargin)
            if nargin == 2
                self.Transparency = varargin{1};
                assert(isnumeric(self.Transparency), 'Input must either be double of key-value pairs!')
                P.FaceAlpha = self.Transparency;
                P.EdgeAlpha = self.Transparency;
                P.MarkerFaceAlpha = self.Transparency;
                P.MarkerEdgeAlpha = self.Transparency;
            else
                P.FaceAlpha = self.FaceAlpha;
                P.EdgeAlpha = self.EdgeAlpha;
                P.MarkerFaceAlpha = self.MarkerFaceAlpha;
                P.MarkerEdgeAlpha = self.MarkerEdgeAlpha;
                P = hdsort.util.parseInputs(P, varargin, 'error');
            end
            
            if ~isempty(self.hdsort.plot.bj) && any(ismember(fieldnames(self.hdsort.plot.bj), {'FaceAlpha', 'EdgeAlpha', 'MarkerFaceAlpha', 'MarkerEdgeAlpha'}))
            %if any([isfield(self.hdsort.plot.bj, 'FaceAlpha'), isfield(self.hdsort.plot.bj, 'EdgeAlpha'),  isfield(self.hdsort.plot.bj, 'MarkerFaceAlpha'), isfield(self.hdsort.plot.bj, 'MarkerEdgeAlpha')])
                for ii = 1:numel(self.hdsort.plot.bj)
                    for f = fields(P)'
                        try
                            set(self.hdsort.plot.bj(ii), f{1}, P.(f{1}));
                        catch
                        end
                    end
                end
            end
        end
        
        % -----------------------------------------------------------------
        function count = figureChildrenSet(self, handle, parameter, value)
            count = 0;
            try
                set(handle, parameter, value);
                count = count +1;
            catch
                % Ignore this
            end
            c = get(handle, 'Children');
            for i=1:length(c)
                count = count + self.figureChildrenSet(c(i), parameter, value);
            end
        end
        
        % -----------------------------------------------------------------
        function dataCursorMode(self, figHandle)
            persistent old_pos;
            persistent old_cval;
            if nargin == 0
                figHandle = self.fh;
            end
            set(datacursormode(figHandle),'UpdateFcn',{@mydataCursor, '%10.1f'});
            
            % ---------------------------------------------------------------------
            function txt = mydataCursor(empt, event_obj, str)
                %global dtip;
                pos = get(event_obj,'Position');
                % check if this is an imagesc plot
                target = get(event_obj, 'Target');
                try
                    cdata = get(target, 'CData');
                    xdata = get(target, 'XData');
                    ydata = get(target, 'YData');
                    xidx = find(abs(xdata-pos(1))<eps);
                    yidx = find(abs(ydata-pos(2))<eps);
                    cval = cdata(yidx, xidx);
                catch
                    cval = [];
                end
                
                px = sprintf(str,pos(1));
                py = sprintf(str,pos(2));
                txt = {['X: ', px],['Y: ', py]};
                if ~isempty(cval)
                    txt = [txt {['C: ', sprintf(str, cval)]}];
                end
                dtip = struct();
                dtip.x = pos(1);
                dtip.y = pos(2);
                assignin('base','dtip',struct());
                fprintf('x: %10.2f y: %10.2f', dtip.x, dtip.y);
                if exist('old_pos','var')
                    try
                        dx = sprintf(str,pos(1)-old_pos(1));
                        dy = sprintf(str,pos(2)-old_pos(2));
                        txt = {txt{:}, ['dX: ', dx], ['dY: ', dy]};
                        dtip.dx = pos(1)-old_pos(1);
                        dtip.dy = pos(2)-old_pos(2);
                        fprintf(' dx: %10.2f dy: %10.2f\n', dtip.dx, dtip.dy);
                    catch
                        % old_pos is not properly initialized. Ignore this.
                        fprintf('\n');
                    end
                else
                    fprintf('\n');
                end
                old_pos = pos;
                old_cval = cval;
                dtip.old_pos = pos;
                assignin('base','dtip',dtip);
            end
        end
        
    end
    
    methods (Access = protected)
        
        % -----------------------------------------------------------------
        function setChild(self, child, idx)
            if any(cellfun(@(x) isequal(x,child) , self.Children)) %any( ismember(self.Children, child) )
                return;
            end
            if nargin == 2
                self.Children = {self.Children child};
            else
                try
                    self.removeChild(self.Children(idx));
                catch
                    % Let pass
                end
                self.Children{idx} = child;
            end
            child.setParent(self);
        end
        
        function removeChild(self, child)
            idx = cellfun(@(x) isequal(x,child) , self.Children); %ismember(self.Children, child);
            if ~any(idx)
                return;
            end
            self.Children{idx} = [];
        end
        
        function setParent(self, parent)
            if ismember(self.Parent, parent)
                return;
            end
            if ~isempty(self.Parent)
                self.Parent.removeChild(self);
            end
            self.Parent = parent;
            self.Parent.setChild(self);
        end
        
    end
    
    methods(Static)
        % -----------------------------------------------------------------
        function [colorVec markerSet C] = vectorColor(nIdx)
            % from: markerSet = mysort.hdsort.plot.lineTypes():
            markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
            %markerSet = {markers{1:(mod(N,length(markers))+1)}};
            
            % Just returns colors
            C =[0.0 0.0 1.0 % 01 blue
                0.0 1.0 0.0 % 02 lime
                1.0 0.0 0.0 % 03 red    
                0.9 0.9 0.0 % 04 yellow
                1.0 0.0 1.0 % 05 fuchsia
                0.0 1.0 1.0 % 06 aqua
                0.5 0.5 1.0 % 07
                1.0 0.5 0.5 %.6 .3 .4      % 08
                0.3 0.7 0.3 % 09
                0.0 0.0 0.5 % 10 navy
                0.0 0.5 0.0 % 11 green
                0.5 0.0 0.0 % 12 maroon
                0.0 0.5 0.5 % 13 teal
                0.5 0.0 0.5 % 14 purple
                0.5 0.5 0.0 % 15 olive
                0.2 0.5 0.8 % 16
                0.8 0.3 0.2
                0.3 0.3 0.3
                0.0 0.0 0.0 % 19 black
                0.7 0.3 0.8
                0.2 0.9 0.9
                0.8 0.8 0.3
                0.5 0.5 0.5 % 07 gray
                0.1 0.3 0.1];
                
            
            assert( size(nIdx, 1) == 1, 'nIdx must be a row vector!')
            N = numel(nIdx);
            nRep = ceil(N/size(C,1));
            colorVec_ = repmat(C, nRep, 1);
            colorVec = colorVec_(nIdx, :);
            
            nRep_markers = ceil(N/numel(markers));
            markerSet_ = repmat(markers, 1, nRep_markers);
            markerSet = {markerSet_{nIdx}};
        end
        
        % -----------------------------------------------------------------
        function [minS, maxS] = getMaxLabelPosition(parent)
            minS = [0, 0, 0];
            maxS = [0, 0, 0];
            for f = fieldnames(parent)'
                try
                    set(parent.(f{:}), 'Units', 'pixels');
                    p = parent.(f{:}).Position;
                catch
                    continue;
                end
                if numel(p) == 3
                    maxS = max(maxS, p);
                    minS = min(minS, p);
                end
            end
        end
        % -----------------------------------------------------------------
        
    end
    
end