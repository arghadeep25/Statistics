function [X Y]=ellipse(a,b,x0,y0,phi)

theta = linspace(0, 360, 300)' .* (pi / 180);
% Parametric equation of the ellipse
%----------------------------------------
 x = a*cos(theta);
 y = b*sin(theta);

% Coordinate transform 
%----------------------------------------
 X = cos(phi)*x - sin(phi)*y;
 Y = sin(phi)*x + cos(phi)*y;
 X = X + x0;
 Y = Y + y0;
 if nargout==1, X = [X Y]; 
 end

