
function dr=dr_dy_target(xi,yi,xk,yk)

s=calcdist(yi,xi,yk,xk);

dr=(xk-xi)/s^2;