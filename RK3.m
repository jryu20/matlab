function [u_n] = RK3(A,B,delta,u_0)
% performs a third-order Runge-Kutta method for an IVP

n = (B-A)/delta;
u_n = zeros(1,n+1);
u_n(1) = u_0;

t = A:delta:B;

F = @(u_n, t) 2*(u_n);

for i=1:n
    k_1 = delta*F(u_n(i),t(i));
    k_2 = delta*F((u_n(i)+0.5*k_1),(t(i)+0.5*delta));
    k_3 = delta*F((u_n(i)+0.5*k_1+0.5*k_2),(t(i)+delta));
    u_n(i+1) = u_n(i)+(1/6)*(k_1+4*k_2+k_3);

end