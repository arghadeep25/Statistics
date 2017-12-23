
function dr=dr_dx_target(xi,yi,xk,yk)

s=calcdist(yi,xi,yk,xk);

dr=-(yk-yi)/s^2;
end