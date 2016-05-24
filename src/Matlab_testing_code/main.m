global A epsilon w Delta dt T L H St R St_inverse;
A = 0.1;
epsilon = 0.25;
w = 0.6*pi;
Delta = 0.000001;
dt = 0.1;
T = 15;
L = 40;
H = 20;
St = 0.2;
R = 1;
St_inverse = 1/St;


final_pos = zeros(H*L,2);
i = 1;
tic;
for y = linspace(0.01,0.99,H)
    for x = linspace(0.01,1.99,L)
        state = [x;y;velocity(x,y,0)];
        %[t,state] = ode45(@update, [0 15], state);
        
        for t = 0:dt:15
            k1 = dt*update(t,state);
            k2 = dt*update(t+0.5*dt,state+0.5*k1);
            k3 = dt*update(t+0.5*dt,state+0.5*k2);
            k4 = dt*update(t+dt,state+k3);
            state = state +(k1+2*k2+2*k3+k4)/6;
        end
        final_pos(i,:) = state(1:2);
        %final_pos(i,:) = state(end,1:2);
        i = i+1;
    end
end   
scatter(final_pos(:,1),final_pos(:,2),5,'filled')
toc;
%{
option = odeset('RelTol',1e-10,'AbsTol',1e-8);
x = 0.4;
y = 0.4;
v = velocity(x,y,0);
u = [x;y;v];
%display(u);
time = linspace(0,15,150);
[t,U]=ode45(@update,time,u,option);
figure
plot(U(:,1),U(:,2),'b-o');
%}

%{
u = [x;y;v];
pts = zeros(4,1501);
i = 1;
for t = 0:dt:15
    k1 = dt*update(t,u);
    k2 = dt*update(t+0.5*dt,u+0.5*k1);
    k3 = dt*update(t+0.5*dt,u+0.5*k2);
    k4 = dt*update(t+dt,u+k3);
    u = u +(k1+2*k2+2*k3+k4)/6;
    pts(1,i) = u(1);
    pts(2,i) = u(2);
    i = i + 1;
end
figure
plot(pts(1,:),pts(2,:));




ax = zeros(20,1);
ay = zeros(20,1);
ax2 = zeros(20,1);
ay2 = zeros(20,1);
i = 1;
for x = linspace(0.1,1.9,20)
    a = accel(x,0.6,0);
    ax(i,1) = a(1);
    ay(i,1) = a(2);
    
    a2 = accel2(x,0.6,0);
    ax2(i,1) = a2(1);
    ay2(i,1) = a2(2);
    
    i = i + 1;
end

%{
figure
subplot(2,1,1)
plot(ax);
subplot(2,1,2)
plot(ay);

figure
subplot(2,1,1)
plot(ax2);
subplot(2,1,2)
plot(ay2);
%}

figure
subplot(2,1,1)
plot(abs(ax-ax2))
subplot(2,1,2)
plot(abs(ay-ay2))
%}