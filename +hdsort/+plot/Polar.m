classdef Polar < plot.PlotInterface
    properties (SetAccess=protected)
        compass
        polar
        
        
    end
    
    properties
        vectors
        theta
        rho
        
        %colormap
        normalize
        closeCurve
        %yShift
        
        %plotAll
        %plotMean
        %plotMedian
        %plotStd
        %meanColor
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = Polar(THETA, RHO, varargin)
            % Convert row vector to column vector:
            %if size(THETA, 1) == 1
            %   THETA = THETA(:);
            %end
            VECTORS = []; 
            if nargin == 1 && size(THETA,2) == 2
                VECTORS = THETA;
                THETA = [];
                RHO = [];
                varargin = {};
            elseif nargin == 1 && size(THETA,2) == 1
                VECTORS = [];
                RHO = THETA;
                THETA = (1:size(RHO,1))';
                varargin = {};
            elseif ischar(RHO) && size(THETA,2) == 2
                VECTORS = THETA;
                varargin = {RHO, varargin{:}};
                THETA = [];
                RHO = [];
            elseif ischar(RHO) && size(THETA,2) == 1
                RHO = THETA;
                varargin = {RHO, varargin{:}};
                VECTORS = [];
            end
            
            if size(THETA, 1) == 1
                THETA = THETA(:);
            end
            if size(RHO, 1) == 1
                RHO = RHO(:);
            elseif size(RHO, 2) == size(THETA, 1)
                RHO = RHO';
            end
            
            % 
            %P.yShift = 0.0;
            %P.plotAll = true;
            %P.plotMean = false;
            %P.plotMedian = false;
            %P.plotStd = false;
            %P.meanColor = [0,0,0];
            P.vectors = VECTORS;
            P.normalize = true;
            P.closeCurve = true;
            self = self@plot.PlotInterface(P, varargin{:});
            
            self.plotName = 'Polar';
            
            if any(THETA > 2*pi)
                THETA = deg2rad(THETA);
            end
            
            self.theta = THETA';
            self.rho = RHO';
            if ~isempty(self.vectors) & size(self.vectors, 2) ~= 2
                assert(size(self.vectors, 1) == 2, 'Vectors must be a 2 column vector!')
                self.vectors = self.vectors';
            end
            
            
            %self.doNotSetAh = true;
            
            self.show();
        end
        
        function show_(self)
            self.setColor(self.color, size(self.rho, 1));
            %cla(self.ah)
            %close all
            %figure;
            hold off
            maxR = 0;
            if ~isempty(self.rho) && ~isempty(self.theta)
                %self.polar = polar(self.ah, self.theta, self.rho);
                theta = self.theta;
                rho = self.rho;
                
                if self.closeCurve
                    theta = [theta, theta(1)];
                    rho   = [rho, rho(:, 1)];
                end
                
                for ii = 1:size(rho, 1)
                    fallback = true;
                    if isempty(self.vectors)
                        try
                            self.polar = ezpolar(theta, rho(ii,:), 'Color', self.color(ii, :) );
                            %self.polar = polarplot(theta, rho(ii,:), 'Color', self.color(ii, :) );
                            hold on;
                            fallback = false;
%                            self.polar = polarplot(theta, rho(ii,:), 'r');
                        catch
                            fallback = true
                        end
                    end
                    if fallback
                        self.polar = polar(theta, rho(ii,:));
                        hold on;
%                        self.polar = polar(theta, rho(ii,:), 'r');
                    end            
                end
                maxR = max(rho(:));
                self.ah = gca;
            end
            if ~isempty(self.vectors)
                if ~isempty(self.rho) && ~isempty(self.theta)
                    hold on;
                end
                
                self.compass = compass(self.ah, self.vectors(:,1), self.vectors(:,2));
                maxR = max([sqrt(self.vectors(:,1).^2 + self.vectors(:,2).^2); maxR]);

                self.setAxis([-1 1 -1 1]*maxR)
            end
            
            if fallback
                self.useDefaultValues();
            end
        end
    end
end