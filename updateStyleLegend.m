function fig = updateStyleLegend(fig,plotName,legendtxt)
% updateStyleLegend.m
% updates the style and legends of a graph
% INPUT 
%   f   figure object
%   plotName  string of the saved graphs
%   legendtxt cell array of strings containing legend texts
% OPTIONAL INPU
%   'landscape' boolean to make portrait (default: true)
% OUTPUT
%   fig   figure object
% Created:
%   Antti Ripatti 30.8.2012
% 
% (c) Antti Ripatti 2012
%
%% style structure for saved graphs
styl = struct();
styl.figure.PaperOrientation = 'Landscape';
styl.figure.PaperType = 'A4';
styl.figure.PaperPositionMode = 'manual';
styl.figure.PaperUnits = 'normalized';
styl.figure.PaperPosition = [0,0,1,1];
styl.figure.Units = 'normalized';
styl.figure.Position = [0,0,1,1];
styl.line.LineWidth = 2;
styl.title.FontSize = 14;
% adding legend
for i=fig;
  figure(i)
  legend(legendtxt,'location','southeast');
end;
for i=1:numel(fig)
  qstyle(fig(i),styl);
  saveas(fig(i),[plotName int2str(i) '.pdf']);
end;
