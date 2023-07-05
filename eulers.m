function [u_n] = eulers(A,B,delta,u_0)

% uses forward Euler's method for an IVP

n = (B-A)/delta;

u_n = zeros(1,n+1);   
u_n(1) = u_0;
t = A:delta:B;

F = @(u_n,t) 2*(u_n);


for i=1:n
    u_n(i+1) = u_n(i) + delta*F(u_n(i), t(i));

end