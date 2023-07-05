function [u_11,u_21] = RK4(A,B,delta)
% performs a fourth-order Runge-Kutta method for an IVP

n = (B-A)/delta;
u_11 = zeros(1,n+1);
u_11(1) = 1;
u_12 = zeros(1,n+1);
u_12(1) = 0;

u_21 = zeros(1,n+1);
u_21(1) = 0;
u_22 = zeros(1,n+1);
u_22(1) = 1;
x = A:delta:B;

F_11 = @(u_11, u_12) u_12;
F_12 = @(u_11, u_12) (x^2-3)*u_11;
F_21 = @(u_21, u_22) u_22;
F_22 = @(u_21, u_22) (x^2-3)*u_21;

for i=1:n
    k_1 = delta*F_11(u_11(i),x(i));
    k_2 = delta*F_11((u_11(i)+0.5*k_1),(x(i)+0.5*delta));
    k_3 = delta*F_11((u_11(i)+0.5*k_2),(x(i)+0.5*delta));
    k_4 = delta*F_11((u_11(i)+k_3),(x(i)+delta));
    u_11(i+1) = u_11(i)+(1/6)*(k_1+2*k_2+2*k_3+k_4);

    k_11 = delta*F_21(u_21(i),x(i));
    k_22 = delta*F_21((u_21(i)+0.5*k_1),(x(i)+0.5*delta));
    k_33 = delta*F_21((u_21(i)+0.5*k_2),(x(i)+0.5*delta));
    k_44 = delta*F_21((u_21(i)+k_3),(x(i)+delta));
    u_21(i+1) = u_21(i)+(1/6)*(k_11+2*k_22+2*k_33+k_44);

end