function ah = ellipse(mu, Cov, varargin)

P.ah = gca;
P.scaling = [1,1];
P.color = [1, 0, 0];
P = hdsort.util.parseInputs(P, varargin, 'error');
            
Rpos = hdsort.util.ellipse(mu, Cov);
%hdsort.plot.ah, Rpos(:,1)*self.SI.squareSize(1)+P.zeroPoint(1), Rpos(:,2)*self.SI.squareSize(2)+P.zeroPoint(2), 'color', P.color);
hdsort.plot.P.ah, Rpos(:,1)*P.scaling(1), Rpos(:,2)*P.scaling(2), 'color', P.color);
end