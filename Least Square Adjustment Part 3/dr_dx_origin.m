function dr=dr_dx_origin(xi,yi,xk,yk)

s=calcdist(xi,yi,xk,yk);

dr=(yk-yi)/s^2;

end