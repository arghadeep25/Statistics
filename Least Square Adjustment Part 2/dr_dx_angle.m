function dr=dr_dx_angle(xi,yi,xk,yk,xl,yl)

sk=dis(xi,yi,xk,yk);
sl=dis(xi,yi,xl,yl);

dr=(yk-yi)/sk^2-(yl-yi)/sl^2;

end