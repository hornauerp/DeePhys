function setAllFontSizes(figure_handle,fontsz)
% SETALLFONTSZIZES  Set font size on all axes, legends, and colorbars in a figure.
%
% INPUTS:
%   figure_handle - Handle to the target figure
%   fontsz        - Desired font size (numeric)
allPolarInFigure = findall(figure_handle,'type','polaraxes');
allAxesInFigure = findall(figure_handle,'type','axes');
allLegendInFigure = findall(figure_handle,'type','legend');
allColorbarInFigure = findall(figure_handle,'type','colorbar');
arrayfun(@(x) set(x,'FontSize',fontsz),allPolarInFigure)
arrayfun(@(x) set(x,'FontSize',fontsz),allAxesInFigure)
arrayfun(@(x) set(x,'FontSize',fontsz),allLegendInFigure)
arrayfun(@(x) set(x,'FontSize',fontsz),allColorbarInFigure)
arrayfun(@(x) set(x,'TitleFontSizeMultiplier',1),allAxesInFigure)
arrayfun(@(x) set(x,'LabelFontSizeMultiplier',1),allAxesInFigure)

arrayfun(@(x) set(x,'TitleFontSizeMultiplier',1),allPolarInFigure)
end