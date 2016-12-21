classdef UiTable < myplot.PlotInterface
    properties (SetAccess=protected)
        
    end
    
    properties % Display settings
        varP
        uiTable
        uiTableExtent
        uiTablePosition
        
        xPos
        yPos
        xExtent
        yExtent
        
        Data
        ColumnName
        ColumnFormat
        ColumnEditable
        ColumnWidth
        RowName
        
    end
    
    methods
        
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = UiTable(varargin)
            P.varP = {};
            P.xPos = 20;
            P.yPos = 20;
            P.xExtent = 200;
            P.yExtent = 200;
            
            P.Data = []; 
            P.ColumnName = '';
            P.ColumnWidth = 'auto';
            P.ColumnFormat = {};
            P.ColumnEditable = [];
            P.RowName = 'numbered';
            
            self = self@myplot.PlotInterface(P, varargin{:})
            
            self.plotName = 'UiTable';
            self.show();
        end
        
        function show_(self)
            
            self.uiTable = uitable(self.fh, ...
                'Data', self.Data, ...
            'ColumnName', self.ColumnName, ...
            'ColumnFormat', self.ColumnFormat, ...
            'ColumnEditable', self.ColumnEditable, ...
            'ColumnWidth', self.ColumnWidth, ...
            'RowName', self.RowName);
            
            self.uiTable.Position(3) = self.uiTable.Extent(3);
            self.uiTable.Position(4) = self.uiTable.Extent(4);
            
            self.uiTableExtent = get(self.uiTable,'Extent');
            self.uiTablePosition = get(self.uiTable,'Position');
            set(self.uiTable, 'Position', [self.xPos self.yPos round(1.0*self.uiTableExtent(3)) round(1.0*self.uiTableExtent(4))]);
            set(self.fh, 'position', [self.fh.Position(1) self.fh.Position(2) 25+round(1.05*self.uiTableExtent(3)) 25+round(1.2*self.uiTableExtent(4))]);
            
            hidem(self.ah)
        end
    end
end