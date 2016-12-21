function resizeTableFigure(uiTable,fig)
if nargin < 1
    % When no inputs are provided, it will try to reformat the current window and
    % its table.
    fig = gcf;
    uiTable = get(gcf,'children');
elseif nargin < 2
    fig = get(uiTable,'parent');
end
try
      uiTableExtent = get(uiTable,'Extent');
      uiTablePosition = get(uiTable,'Position');
      set(uiTable,'Position',[20 20 round(1.0*uiTableExtent(3)) round(1.0*uiTableExtent(4))]);
      set(fig,'position',[200 200 round(1.2*uiTableExtent(3)) round(1.2*uiTableExtent(4))]);
catch
      warning('UiTable figure could not be resized!')
end