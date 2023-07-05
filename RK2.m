function [u_final,v_final,w_final] = RK2(a,b,N,u_0,v_0,w_0)
%implements the Runge-Kutta 2 scheme to solve a system of first order IVPs

h = (b-a)/N;
t=a;

u = zeros(1,N+1);
u(1) = u_0;
v = zeros(1,N+1);
v(1) = v_0;
w = zeros(1,N+1);
w(1) = w_0;

f1 = @(t,u,v,w) v;
f2 = @(t,u,v,w) w;
f3 = @(t,u,v,w) (8-(2/(t^3)))+(2/(t^2))*v-(1/t)*w;

for i=1:N
    k1 = h*f1(t,u(i),v(i),w(i));
    l1 = h*f2(t,u(i),v(i),w(i));
    m1 = h*f3(t,u(i),v(i),w(i));
    k2 = h*f1(t+h,u(i)+k1,v(i)+l1,w(i)+m1);
    l2 = h*f2(t+h,u(i)+k1,v(i)+l1,w(i)+m1);
    m2 = h*f3(t+h,u(i)+k1,v(i)+l1,w(i)+m1);
    u(i+1) = u(i)+0.5*k1+0.5*k2;
    v(i+1) = v(i)+0.5*l1+0.5*l2;
    w(i+1) = w(i)+0.5*m1+0.5*m2;
    t = a+i*h;


end

u_final = u(N+1);
v_final = v(N+1);
w_final = w(N+1);

end