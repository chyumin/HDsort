classdef Boxplot < hdsort.plot.PlotInterface
    properties (SetAccess=protected)
        boxplot
        
        nGroups
        nDataPoints
        nSubGroups
        
        grouplabels
        subgrouplabels
        
        groupIdx
        groups
        data
        reshaped_data
    end
    
    properties % Display settings
        displayN
        subgroupnames
        groupnames
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Boxplot(data, groupIdx, varargin)
            % Data a NxS matrix containing N datapoints and S dimensions (subgroups)
            % each.
            if nargin == 1
                groupIdx = [];
                varargin = {};
            elseif ischar(groupIdx)
                varargin = {groupIdx, varargin{:}};
                groupIdx = [];
            end
            %assert(nargin > 1, 'Specify groupIdx, [] is possible.')
            
            
            P.displayN = true;
            P.subgroupnames = {};
            P.groupnames = {};
            self = self@hdsort.plot.PlotInterface(P, varargin{:})
            
            self.plotName = 'Boxplot';
            
            self.data = data;
            if size(self.data, 1) == 1
                self.data = self.data';
            end
            
            if isempty(groupIdx)
                assert(isempty(self.groupnames), 'You can not specify groupnames without grouplabels!')
                self.groupIdx = ones(1, size(self.data,1));
            else
                self.groupIdx = groupIdx;
            end
            
            self.groups = unique(self.groupIdx);
            self.nGroups = length(self.groups);
            [nDataPoints_ nSubGroups_] = size(self.data);
            self.nDataPoints = nDataPoints_;
            self.nSubGroups = nSubGroups_;
            
            %assert(self.nGroups == numel(self.groupnames), 'Number of groups must be the same as the number of groupnames!')
            
            %% Create labels for groups:
            if isempty(self.groupnames)
                for g = 1:self.nGroups
                    self.groupnames{g} = num2str(self.groups(g));
                end
            end
            self.grouplabels = repmat(self.groupnames,1,self.nSubGroups);
            
            try
                assert(numel(self.groupnames) == self.nGroups, 'Number of grouplabels must be the same as the number of groups!')
            catch ME
                try
                    %% When the groups miss one element but their indices are otherwise in order:
                    self.groupnames = self.groupnames(self.groups);
                    self.grouplabels = {self.grouplabels{self.groups}};
                catch
                    rethrow(ME)
                end
            end
            
            %% Create labels for subgroups:
            if isempty(self.subgroupnames)
                for sg = 1:self.nSubGroups
                    self.subgroupnames{sg} = num2str(sg);
                end
            end
            assert(numel(self.subgroupnames) == self.nSubGroups, 'Number of subgroupnames must be the same as the number of subgroups!')
            
            self.subgrouplabels = [];
            for sg = 1:self.nSubGroups
                self.subgrouplabels = [self.subgrouplabels, repmat({self.subgroupnames{sg}},1,self.nGroups)];
            end
            
            self.color = hdsort.plot.vectorColor(1:self.nSubGroups);
            
            %% Reshape data:
            self.reshaped_data = zeros(self.nDataPoints, self.nGroups*self.nSubGroups)*NaN;
            for g = 1:self.nGroups
                dataGroup = self.data(self.groupIdx == self.groups(g), :);
                [L nSubGroups_] = size(dataGroup);
                assert(self.nSubGroups==nSubGroups_, 'error')
                
                for sg = 1:self.nSubGroups
                    ii =(sg-1)*self.nGroups + g;
                    self.reshaped_data(1:L,ii) = dataGroup(1:L,sg);
                end
            end
            
            self.show();
        end
        
        % -----------------------------------------------------------------
        function show_(self)
            
            if self.nSubGroups == 1
                self.boxplot = boxplot(self.reshaped_data, self.grouplabels);
            else
                self.boxplot = boxplot(self.reshaped_data, {self.grouplabels, self.subgrouplabels}, ...
                    'colors', self.color, 'factorgap',[5 2], 'labelverbosity', 'minor');
            end
            
            if self.displayN
                M = nanmean(self.reshaped_data);
                N = sum(~isnan(self.reshaped_data));
                y = max(self.reshaped_data(:));
                
                z_ = findobj( self.ah, 'Tag', 'Box'); z = z_(end:-1:1);
                
                for k1 = 1:size(self.reshaped_data,2)
                    text( z(k1).XData(1), y, sprintf('N = %d', N(k1)), 'FontSize', self.FontSize);
                    %text( z(size(self.reshaped_data,2)+1 - k1).XData(1), y, sprintf('N = %d', N(k1)), 'FontSize', self.FontSize);
                end
            end
            
            self.XTick = get(self.ah, 'XTick');
            self.XTickLabel = get(self.ah, 'XTickLabel');
            self.XTickMode = 'manual';
            self.XTickLabelMode = 'manual';
            
            if isempty(self.YTick)
                self.YTick = get(self.ah, 'YTick');
            else
                self.YTickMode = 'manual';
            end
            
            if isempty(self.YTickLabel)
                self.YTickLabel = get(self.ah, 'YTickLabel');
            else
                self.YTickLabelMode = 'manual';
            end
            
        end
    end
end
