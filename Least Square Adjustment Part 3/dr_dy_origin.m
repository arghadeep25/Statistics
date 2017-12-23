function dr=dr_dy_origin(xi,yi,xk,yk)

s=calcdist(xi,yi,xk,yk);

dr=-(xk-xi)/s^2;

end