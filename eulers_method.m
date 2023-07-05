function [u_0_final, u_1_final] = eulers_method(A,B,delta)

% solves a system of 2 IVPs using Euler's method

n = (B-A)/delta;

u_0 = zeros(1,n+1);   
u_0(1) = 4;
u_1 = zeros(1,n+1);   
u_1(1) = 1;

F_0 = @(u_0,u_1) u_1;
F_1 = @(u_0,u_1) 4*(u_0) - cos(u_1) + 7*n*delta;


for i=1:n
    u_0(i+1) = u_0(i) + delta*F_0(u_0(i), u_1(i));
    u_1(i+1) = u_1(i) + delta*F_1(u_0(i), u_1(i));

u_0_final = u_0(n);
u_1_final = u_1(n);
end