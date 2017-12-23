% Function to calculate the direction out of coordinates
% direction ik
function direction=direction(yk, xk, yi, xi, w)

  direction=atan2((yk-yi),(xk-xi));
  
  if direction<0
      direction=direction+2*pi;
  end
  direction=direction-w;
  
end