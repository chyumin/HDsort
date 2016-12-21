classdef Subplots < plot.PlotInterface
    properties (SetAccess=protected)
        nSubplots
        nRows
        nColumns
        subplotHandles
        previousIdx
    end
    
    properties
        offsetX
        offsetY
        spacerTop
        spacerRight
        spacerX
        spacerY
        labels
        preferVertical
        
        matrix
        holdOn
        %parentFh
        upperTriangle
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Subplots(nSubplots, varargin)
            
            P.offsetX = .06;
            P.offsetY = .06;
            P.spacerTop   = .08;
            P.spacerRight = .0;
            P.spacerX = .05;
            P.spacerY = .08;
            P.labels  = {};
            P.preferVertical = 1;
            
            P.matrix = 0;
            P.holdOn = 1;
            P.upperTriangle = false;
            self = self@plot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Subplots';
            
            self.nSubplots = nSubplots;
            self.previousIdx = 0;
            
            %% Set the subplot sizes:
            if length(self.nSubplots) == 2
                self.nRows = self.nSubplots(1);
                self.nColumns = self.nSubplots(2);
                self.nSubplots = self.nRows*self.nColumns;
            elseif self.preferVertical
                self.nColumns = floor(sqrt(self.nSubplots));
                self.nRows = ceil(self.nSubplots/self.nColumns);
            else
                self.nRows = floor(sqrt(self.nSubplots));
                self.nColumns = ceil(self.nSubplots/self.nRows);
            end
            
            self.showAh = false;
            self.show();
        end
        
        % -----------------------------------------------------------------
        function ah = getSubplotHandle(self, idx)
            if nargin == 2
                ah = self.structAh(self.subplotHandles(idx));
            else
                for ii = 1:length(self.subplotHandles);
                    ah(ii) = self.structAh(self.subplotHandles(ii));
                end
            end
        end
        
        % -----------------------------------------------------------------
        function mode = linkaxes(self, mode)
            if nargin == 1 mode = 'xy'; end
            linkaxes(self.getSubplotHandle(), mode); 
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            %hidem(self.ah);
            
            %if isempty(self.parentFh)
            %    self.parentFh = gcf;
            %end
            
            assert(isempty(self.labels) || length(self.labels) >= self.nSubplots, 'If labels is provided it must contain one label per subplot!');
            assert(self.nRows < 37 && self.nColumns <37, 'Too many rows or columns to plot!');
            
            w = max(0.01, (.99-self.offsetX-(self.nColumns-1)*self.spacerX-self.spacerRight)/(self.nColumns));
            h = max(0.01, (.99-self.offsetY-(self.nRows-1)*self.spacerY-self.spacerTop)/(self.nRows));
            
            % check if spacers are too big
            maxHeight = self.offsetY + (self.nRows-1)*(self.spacerY + h);
            if maxHeight > 1
                self.spacerY = 0;
                h = max(0.01, (1-2*self.offsetY)/(self.nRows));
            end
            maxWidth = self.offsetX + (self.nColumns-1)*(self.spacerX + w);
            if maxWidth > 1 || self.spacerX > w
                self.spacerX = 0;
                w = max(0.01, (1-1.5*self.offsetX-(self.nColumns-1)*self.spacerX)/(self.nColumns));
            end
            
          %  if ~isempty(self.title)
          %      % Correct hight, to leave room for title
          %      h = max(0.01, h - .05/self.nRows);
          %      self.setFigureTitle();
          %  end
            
            % calculate XCoords and YCoords:
            count = 1; xC = zeros(self.nColumns,1); yC = zeros(self.nRows,1);
            self.subplotHandles = []; %zeros(self.nSubplots,1);
            %self.subplotIdx = [];
            for j=self.nRows:-1:1
                for i=1:self.nColumns
                    if count > self.nSubplots
                        continue
                    end
                    if ~self.upperTriangle || i >= (self.nRows-j+1)
                        xC(i) = self.offsetX + (i-1)*(self.spacerX + w);
                        yC(j) = self.offsetY + (j-1)*(self.spacerY + h);
                        
                        
                        if isempty(self.subplotHandles)
                            self.subplotHandles = axes('Units','normalized','position',[xC(i) yC(j) w h]);
                        else
                            self.subplotHandles(count) = axes('Units','normalized','position',[xC(i) yC(j) w h]);
                        end
                        
                        %self.subplotIdx(count) = count;
                            %self.subplotHandles(count) = axes('Units','normalized','position',[xC(i) yC(j) w h], 'parent', self.parentFh);
                        
                        if self.holdOn
                            set(self.subplotHandles(count), 'nextplot', 'add');
                        end
                        
                        if ~isempty(self.labels)
                            set(self.subplotHandles(count),'Units', 'pixel');
                            pos = get(self.subplotHandles(count), 'Position' );
                            set(self.subplotHandles(count),'Units', 'normalized');
                            
                            ann = annotation('textbox',[0 .1 0 .1],...
                                'Units','pixel',...
                                'FitBoxToText','off',...
                                'String',{['\textbf{' self.labels{count} '}']},...
                                'FontSize', 14,...
                                'Interpreter', 'latex',...
                                'LineStyle','none');%, 'parent', self.parentFh);
                            
                            set(ann, 'Position', [max(0,pos(1)-48) pos(2)+pos(4)+5 30 30]);
                            set(ann, 'Unit', 'normalized');
                        end
                        count = count +1;
                    end
                end
            end
            if self.upperTriangle
                %self.subplotHandles = self.subplotHandles(self.subplotHandles>0);
            end
            if self.matrix
                if self.upperTriangle
                    tmp = self.subplotHandles;
                    self.subplotHandles = nan([self.nColumns self.nRows]);
                    count = 0;
                    for i=1:self.nColumns
                        for j=i:self.nRows
                            count = count+1;
                            self.subplotHandles(i,j) = tmp(count);
                        end
                    end
                else
                    self.subplotHandles = reshape(self.subplotHandles, [self.nColumns self.nRows])';
                end
            end
            
            n = 1;
            for c = self.Children
                self.add(c{1}, n);
                n = n + 1;
            end
            
            self.setFigureTitle();
        end
        
        
        function add(self, PI, idx)
            if isempty(PI)
                return
            end
            %%
            if nargin == 2
                idx = self.previousIdx + 1;
            end
            
            if any(size(idx) ~= [1, 1])
                n = sub2ind([self.nRows, self.nColumns], idx(1), idx(2));
            else
                n = idx;
            end
            
            if ~isequal(PI.fh, self.fh)
                PI.closeFigure();
            else
                PI.closeAxis();
            end
            
            PI.setAh(self.subplotHandles(n));
            PI.show();
            self.setChild(PI, n);
            
            self.previousIdx = n;
        end
        
    end
end