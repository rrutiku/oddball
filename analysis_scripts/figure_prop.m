% general properties
iFontSize       = 20;
strFontUnit     = 'points';     % [{points} | normalized | inches | centimeters | pixels]
strFontName     = 'Times';      % [Times | Courier | ] TODO complete the list
strFontWeight   = 'normal';     % [light | {normal} | demi | bold]
strFontAngle    = 'normal';     % [{normal} | italic | oblique] ps: only for axes
strInterpreter  = 'latex';      % [{tex} | latex]
fLineWidth      = 1.0;          % width of the line of the axes

% note: it is quite difficult to use the "latex" interpreter for the ticks;
% if absolutely needed google for "format_ticks.m" by Alexander Hayes

set(gca, ...
... 'Position', [1 1 20 10], ... TODO
... 'OuterPosition', [1 1 20 10], ... TODO
...
... 'XGrid', 'on', ... [on | {off}]
... 'YGrid', 'on', ... [on | {off}]
... 'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
... 'XMinorGrid', 'off' , ... [on | {off}]
... 'YMinorGrid', 'off', ... [on | {off}]
... 'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
...
'XTick', -2.0:0.1:2.0, ... ticks of x axis
'YTick', -10:1:10, ... ticks of y axis
'XTickLabel', '', ...
'YTickLabel', '', ...
'XMinorTick', 'off' , ... [on | {off}]
'YMinorTick', 'off', ... [on | {off}]
'TickDir', 'in', ... [{in} | out] inside or outside (for 2D)
'TickLength', [.01 .01], ... length of the ticks
...
'XColor', [.1 .1 .1], ... color of x axis
'YColor', [.1 .1 .1], ... color of y axis
'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
... 'XLim', [0 100], ... limits for the x-axis
... 'YLim', [-10 10], ... limits for the y-axis
...
'FontName', strFontName, ... kind of fonts of labels
'FontSize', iFontSize, ... size of fonts of labels
'FontUnits', strFontUnit, ... units of the size of fonts
'FontWeight', strFontWeight, ... weight of fonts of labels
'FontAngle', strFontAngle, ... inclination of fonts of labels
...
'LineWidth', fLineWidth); % width of the line of the axes

% 
% strXLabel = 'label of x axis';
% strYLabel = 'label of y axis';
% %
% fXLabelRotation = 0.0;
% fYLabelRotation = 90.0;
% 
% xlabel( strXLabel, ...
% 'FontName', strFontName, ...
% 'FontUnit', strFontUnit, ...
% 'FontSize', iFontSize, ...
% 'FontWeight', strFontWeight, ...
% 'Interpreter', strInterpreter);
% %
% ylabel( strYLabel, ...
% 'FontName', strFontName, ...
% 'FontUnit', strFontUnit, ...
% 'FontSize', iFontSize, ...
% 'FontWeight', strFontWeight, ...
% 'Interpreter', strInterpreter);
% %
% set(get(gca, 'XLabel'), 'Rotation', fXLabelRotation);
% set(get(gca, 'YLabel'), 'Rotation', fYLabelRotation);

% in order to make matlab to do not "cut" latex-interpreted axes labels
% set(gca, 'Units', 'normalized', ...
% 'Position', [0.15 0.2 0.75 0.7]);


% note: the xPosition and yPosition are referred in the "axes" units of measure
text( ... xPosition, ...
... yPosition, ...
... 'this is the text plotted', ...
'FontName', strFontName, ...
'FontSize', iFontSize, ...
'FontWeight', strFontWeight, ...
'Interpreter', strInterpreter);