function [solplot] = lin_shooting(A,B,alpha,beta,h)
%uses the linear shooting method
%the algorithm for Runge-Kutta was mainly motivated by the pseudocode
%in the Burden and Faires textbook, Section 11.1

p = @(x) 0;
q = @(x) x^2-3;
r = @(x) 0;
n = (B-A)/h;

u = zeros(2,n+1);
v = zeros(2,n+1);
W1 = zeros(1,n+1);
W2 = zeros(1,n+1);

u(1,1) = 1;
u_1 = alpha;
u_2 = 0;
v_1 = 0;
v_2 = 1;

for i=1:n    
    x = A+(i-1)*h;
    k11 = h*u_2;
    k12 = h*(p(x)*u_2+q(x)*u_1+r(x));
    k21 = h*(u_2+0.5*k12);
    k22 = h*(p(x+0.5*h)*(u_2+0.5*k12)+q(x+0.5*h)*(u_1+0.5*k11)+r(x+0.5*h));
    k31 = h*(u_2+0.5*k22);
    k32 = h*(p(x+0.5*h)*(u_2+0.5*k22)+q(x+0.5*h)*(u_1+0.5*k21)+r(x+0.5*h));
    k41 = h*(u_2+k32);
    k42 = h*(p(x+h)*(u_2+k32)+q(x+h)*(u_1+k31)+r(x+h));
    u_1 = u_1+(k11+2*k21+2*k31+k41)/6;
    u_2 = u_2+(k12+2*k22+2*k32+k42)/6;

    k11 = h*v_2;
    k12 = h*(p(x)*v_2+q(x)*v_1);
    k21 = h*(v_2+0.5*k12);
    k22 = h*(p(x+0.5*h)*(v_2+0.5*k12)+q(x+0.5*h)*(v_1+0.5*k11));
    k31 = h*(v_2+0.5*k22);
    k32 = h*(p(x+0.5*h)*(v_2+0.5*k22)+q(x+0.5*h)*(v_1+0.5*k21));
    k41 = h*(v_2+k32);
    k42 = h*(p(x+h)*(v_2+k32)+q(x+h)*(v_1+k31));
    v_1 = v_1+(k11+2*k21+2*k31+k41)/6;
    v_2 = v_2+(k12+2*k22+2*k32+k42)/6;

    u(1,i+1) = u_1;
    u(2,i+1) = u_2;
    v(1,i+1) = v_1;
    v(2,i+1) = v_2;
end

w_1 = alpha;
w_2 = (beta-u(1,n+1))/(v(1,n+1));
x_list = zeros(1,n+1);

for i=1:n+1
    W1(i) = u(1,i)+w_2*v(1,i);
    W2(i) = u(2,i)+w_2*v(2,i);
    x_list(i) = A+(i-1)*h;
end

solplot = plot(x_list,W1);