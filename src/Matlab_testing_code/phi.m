
function val = phi(x,y,t)
    global w A epsilon;
    f = epsilon*sin(w*t)*x*x + (1-2*epsilon*sin(w*t))*x;
    val = A*sin(pi*f)*sin(pi*y);
end