function dsdx=ds_dy_origin(xi,yi,xk,yk)

s=calcdist(xi,yi,xk,yk);

dsdx=-(yk-yi)/s; 