testX = [0.5,0,-0.5,0,0,0.25,0,-0.25,0];
testY = [0,-0.5,0,0.5,0,0,-0.25,0,0.25];
testZ = ones(9,1)';
% Create polar data
[t,r] = meshgrid((0:1:360)*pi/180,[0 50 100]);
% Convert to Cartesian
[X,Y] = pol2cart(t,r);
F = scatteredInterpolant(testX',testY',testZ');
Z = F(X,Y);
% Plot a graph.
figure
contourf(X,Y,Z,30,'LineColor', 'none');
daspect([1 1 1])
% Colormap
c = hsv;
c = c(1:45,:);
colormap(c);