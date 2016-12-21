classdef Subhdsort.plot. < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        nSubhdsort.plot.
        nRows
        nColumns
        subhdsort.plot.andles
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
        function self = Subhdsort.plot.(nSubhdsort.plot., varargin)
            
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
            self = self@hdsort.plot.PlotInterface(P, varargin{:});
            
            self.hdsort.plot.ame = 'Subhdsort.plot.';
            
            self.nSubhdsort.plot. = nSubhdsort.plot.;
            self.previousIdx = 0;
            
            %% Set the subhdsort.plot.sizes:
            if length(self.nSubhdsort.plot.) == 2
                self.nRows = self.nSubhdsort.plot.(1);
                self.nColumns = self.nSubhdsort.plot.(2);
                self.nSubhdsort.plot. = self.nRows*self.nColumns;
            elseif self.preferVertical
                self.nColumns = floor(sqrt(self.nSubhdsort.plot.));
                self.nRows = ceil(self.nSubhdsort.plot./self.nColumns);
            else
                self.nRows = floor(sqrt(self.nSubhdsort.plot.));
                self.nColumns = ceil(self.nSubhdsort.plot./self.nRows);
            end
            
            self.showAh = false;
            self.show();
        end
        
        % -----------------------------------------------------------------
        function ah = getSubhdsort.plot.andle(self, idx)
            if nargin == 2
                ah = self.structAh(self.subhdsort.plot.andles(idx));
            else
                for ii = 1:length(self.subhdsort.plot.andles);
                    ah(ii) = self.structAh(self.subhdsort.plot.andles(ii));
                end
            end
        end
        
        % -----------------------------------------------------------------
        function mode = linkaxes(self, mode)
            if nargin == 1 mode = 'xy'; end
            linkaxes(self.getSubhdsort.plot.andle(), mode); 
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            %hidem(self.ah);
            
            %if isempty(self.parentFh)
            %    self.parentFh = gcf;
            %end
            
            assert(isempty(self.labels) || length(self.labels) >= self.nSubhdsort.plot., 'If labels is provided it must contain one label per subhdsort.plot.');
            assert(self.nRows < 37 && self.nColumns <37, 'Too many rows or columns to hdsort.plot.');
            
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
            self.subhdsort.plot.andles = []; %zeros(self.nSubhdsort.plot.,1);
            %self.subhdsort.plot.dx = [];
            for j=self.nRows:-1:1
                for i=1:self.nColumns
                    if count > self.nSubhdsort.plot.
                        continue
                    end
                    if ~self.upperTriangle || i >= (self.nRows-j+1)
                        xC(i) = self.offsetX + (i-1)*(self.spacerX + w);
                        yC(j) = self.offsetY + (j-1)*(self.spacerY + h);
                        
                        
                        if isempty(self.subhdsort.plot.andles)
                            self.subhdsort.plot.andles = axes('Units','normalized','position',[xC(i) yC(j) w h]);
                        else
                            self.subhdsort.plot.andles(count) = axes('Units','normalized','position',[xC(i) yC(j) w h]);
                        end
                        
                        %self.subhdsort.plot.dx(count) = count;
                            %self.subhdsort.plot.andles(count) = axes('Units','normalized','position',[xC(i) yC(j) w h], 'parent', self.parentFh);
                        
                        if self.holdOn
                            set(self.subhdsort.plot.andles(count), 'nexthdsort.plot., 'add');
                        end
                        
                        if ~isempty(self.labels)
                            set(self.subhdsort.plot.andles(count),'Units', 'pixel');
                            pos = get(self.subhdsort.plot.andles(count), 'Position' );
                            set(self.subhdsort.plot.andles(count),'Units', 'normalized');
                            
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
                %self.subhdsort.plot.andles = self.subhdsort.plot.andles(self.subhdsort.plot.andles>0);
            end
            if self.matrix
                if self.upperTriangle
                    tmp = self.subhdsort.plot.andles;
                    self.subhdsort.plot.andles = nan([self.nColumns self.nRows]);
                    count = 0;
                    for i=1:self.nColumns
                        for j=i:self.nRows
                            count = count+1;
                            self.subhdsort.plot.andles(i,j) = tmp(count);
                        end
                    end
                else
                    self.subhdsort.plot.andles = reshape(self.subhdsort.plot.andles, [self.nColumns self.nRows])';
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
            
            PI.setAh(self.subhdsort.plot.andles(n));
            PI.show();
            self.setChild(PI, n);
            
            self.previousIdx = n;
        end
        
    end
end