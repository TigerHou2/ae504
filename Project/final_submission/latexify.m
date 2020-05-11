function latexify(w,h)
%LATEXIFY Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 2
    w = 18;
    h = 18;
end

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','ticklabelinterpreter'),...
                            'ticklabelinterpreter','latex')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')

set( gca, 'Color', [1 1 1] )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [w h])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0.2 1.2 w h])
set(gcf, 'PaperPosition', [0.2 1.2 w h])

end

