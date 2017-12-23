function dsdx=ds_dy_target(xi,yi,xk,yk)

s=calcdist(xi,yi,xk,yk);

dsdx=(yk-yi)/s;