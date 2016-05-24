function val = accel2(x,y,t)
    global epsilon w A;
    %x = r(1);
    %y = r(2);
    %St = 0.2;
    %global A epsilon w;
    %St_inverse = 1/St;
    %u_accel = (velocity(x,y,t + Delta) - velocity(x,y,t - Delta))/(2*Delta);
    f = epsilon*sin(w*t)*x*x + (1-2*epsilon*sin(w*t))*x; 
    ax = A*pi^2*epsilon*w*cos(pi*y)*cos(pi*f)*cos(w*t)*(x^2-2*x);
    ay = A*pi^2*sin(pi*f)*sin(pi*y)*(2*epsilon*sin(w*t)*x + 1 - 2*epsilon*sin(w*t))*epsilon*w*cos(w*t)*(x^2-2*x);
    ay = ay - 2*A*pi*w*epsilon*cos(w*t)*cos(pi*f)*sin(pi*y)*(x-1);
    %ax = -1*w*sin(w*t);
    %ay = w*cos(w*t);

    val = [ax;ay];
end
