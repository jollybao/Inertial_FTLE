function val = accel(x,y,t)

    global A epsilon w R Delta;
    u_accel = (velocity(x,y,t + Delta) - velocity(x,y,t - Delta))/(2*Delta);
    v = velocity(x,y,t);
    
    u_ax = (velocity(x+Delta,y,t) - velocity(x-Delta,y,t))/(2*Delta);
    u_ay = (velocity(x,y+Delta,t) - velocity(x,y-Delta,t))/(2*Delta);
    u_accel = u_accel + u_ax*v(1) + u_ay*v(2);

    %f = epsilon*sin(w*t)*x*x + (1-2*epsilon*sin(w*t))*x; 
    %ax = A*pi^2*epsilon*w*cos(pi*y)*cos(pi*f)*cos(w*t)*(x^2-2*x);
    %ay = A*pi^2*sin(pi*f)*sin(pi*y)*(2*epsilon*sin(w*t)*x + 1 - 2*epsilon*sin(w*t))*epsilon*w*cos(w*t)*(x^2-2*x);
    %ay = ay - 2*A*pi*w*epsilon*cos(w*t)*cos(pi*f)*sin(pi*y)*(x-1);
    
    val = 1.5*R*u_accel;

end
