function val = update(t,u)
    %global w;
    
    %vx = cos(w*t);
    %vy = sin(w*t);
    
    %ax = -1*w*sin(w*t);
    %ay = w*cos(w*t);
    x = u(1);
    y = u(2);
    a = accel(x,y,t);
    v = velocity(x,y,t);

    %val = [vx; vy; 1; 1];
    %val = [u(3);u(4);ax;ay];
    val = [u(3);u(4); a];
    %val = [v;0;0];
    
end    