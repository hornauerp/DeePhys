function setAllFontSizes(figure_handle,fontsz)
allAxesInFigure = findall(figure_handle,'type','axes');
allLegendInFigure = findall(figure_handle,'type','legend');
allColorbarInFigure = findall(figure_handle,'type','colorbar');
arrayfun(@(x) set(x,'FontSize',fontsz),allAxesInFigure)
arrayfun(@(x) set(x,'FontSize',fontsz),allLegendInFigure)
arrayfun(@(x) set(x,'FontSize',fontsz),allColorbarInFigure)
arrayfun(@(x) set(x,'TitleFontSizeMultiplier',1),allAxesInFigure)
arrayfun(@(x) set(x,'LabelFontSizeMultiplier',1),allAxesInFigure)
end