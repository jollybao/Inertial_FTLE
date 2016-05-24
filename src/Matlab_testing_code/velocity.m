function val = velocity(x,y,t)
    global w A epsilon;
    f = epsilon*sin(w*t)*x*x + (1-2*epsilon*sin(w*t))*x;
    vx = A*pi*sin(pi*f)*cos(pi*y);
    vy = -1*A*pi*cos(pi*f)*sin(pi*y)*(2*epsilon*sin(w*t)*x + (1-2*epsilon*sin(w*t)));
    %vx = (phi(x,y+Delta,t)-phi(x,y-Delta,t))/(2*Delta);
    %vy = (phi(x+Delta,y,t)-phi(x-Delta,y,t))/(2*Delta);
    
    val = [vx; vy];
    %val = [cos(w*t); sin(w*t)];
end   