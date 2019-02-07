function [error]=calcErrorPoints(X, points)
error=norm(calcDistance(points*points')-calcDistance(X), 'fro')/norm(calcDistance(points*points'), 'fro');
end