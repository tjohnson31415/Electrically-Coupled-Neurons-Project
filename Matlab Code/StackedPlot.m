function StackedPlot(X, Y, ScaleLegend)

if nargin == 2
    ScaleLegend = {};
end

% Separation between adjacent plots
sep = 0;
sep = std(Y(:));

maxes = max(Y);
mins  = min(Y);

maxes = maxes(1:end-1);
mins  = mins(2:end);

offsets = maxes - mins + sep;
offsets = [0 offsets];
offsets = cumsum(offsets);

Y = bsxfun(@plus, Y, offsets);

plot(X, Y, 'LineWidth', 2);
set(gca, 'yticklabel', [])
set(gca, 'ytick', [])

hold on
if ~isempty(ScaleLegend)
    position = ScaleLegend{1};
    length = ScaleLegend{2};
    label = ScaleLegend{3};
    
    Xs = [position(1) position(1)];
    Ys = [position(2) position(2) + length];
    
    plot(Xs, Ys, 'k');
    text(Xs(1)+.2, Ys(1)+length/2, [num2str(length) ' ' label]);
end
hold off