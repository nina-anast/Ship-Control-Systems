function Corrcolormap(Data,myLabel)
%% This function plots correlation coefficients in a color map
%% Input Data and names of the variables myLabel

% Check if Data is table and if so convert to array
if istable(Data) == 1
    Data = table2array(Data);
end
% Calculate correlation coefficients of th
Cor = corrcoef(Data);
Cor = tril(Cor,-1);
Cor(logical(eye(size(Cor)))) = 1;
% Set [min,max] value of C to scale colors
clrLim = [-1,1];
% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];
% Compute center of each circle
% This assumes the x and y values were not entered in imagesc()
x = 1 : 1 : size(Cor,2); % x edges
y = 1 : 1 : size(Cor,1); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(Cor==0)=nan; % eliminate cordinates for zero correlations
% Set color of each rectangle
% Set color scale
cmap = jet(256);
Cscaled = (Cor - clrLim(1))/range(clrLim); % always [0:1]
colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size between [0 1]
Cscaled = (abs(Cor) - 0)/1;
diamSize = Cscaled * range(diamLim) + diamLim(1);
% Create figure
fh = figure();
ax = axes(fh);
hold(ax,'on')
colormap(ax,'jet');
tickvalues = 1:length(Cor);
x = zeros(size(tickvalues));
text(x, tickvalues, myLabel, 'HorizontalAlignment', 'right');
x(:) = length(Cor)+1;
text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',90);
% Create circles
theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:),'LineStyle','none'),1:numel(xAll));
axis(ax,'equal')
axis(ax,'tight')
set(ax,'YDir','Reverse')
colorbar()
clim(clrLim);
axis off
end