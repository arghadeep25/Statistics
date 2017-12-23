function dsdx=ds_dx_target(xi,yi,xk,yk)

s=calcdist(xi,yi,xk,yk);

dsdx=(xk-xi)/s;