function dsdx=ds_dx_origin(xi,yi,xk,yk)

s=calcdist(xi,yi,xk,yk);

dsdx=-(xk-xi)/s; 