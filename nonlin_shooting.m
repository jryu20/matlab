function [solplot] = nonlin_shooting(A,B,alpha,beta,h)
%uses the nonlinear shooting method
%the algorithm for Runge-Kutta was mainly motivated by the pseudocode
%in the Burden and Faires textbook, Section 11.2

f = @(u,u_prime,x) -4*u+cos(x);
f_u = @(u,u_prime,x) -4;
f_u_prime = @(u,u_prime,x) 0;

n = (B-A)/h;
TK = (beta-alpha)/(B-A);
k=1;

w1 = zeros(1,n+1);
w2 = zeros(1,n+1);
x_list = zeros(1,n+1);

while k<=10
    w1(1) = alpha;
    w2(1) = TK;
    u_1 = 0;
    u_2 = 1;

    for i=1:n
        x = A+(i-1)*h;
        x_list(i) = x;
        k11 = h*w2(i);
        k12 = h*f(x,w1(i),w2(i));
        k21 = h*(w2(i)+0.5*k12);
        k22 = h*f(x+0.5*h,w1(i)+0.5*k11,w2(i)+0.5*k12);
        k31 = h*(w2(i)+0.5*k22);
        k32 = h*f(x+0.5*h,w1(i)+0.5*k21,w2(i)+0.5*k22);
        k41 = h*(w2(i)+k32);
        k42 = h*f(x+h,w1(i)+k31,w2(i)+k32);
        w1(i+1) = w1(i)+(k11+2*k21+2*k31+k41)/6;
        w2(i+1) = w2(i)+(k12+2*k22+2*k32+k42)/6;

        k11 = h*u_2;
        k12 = h*(f_u(x,w1(i),w2(i))*u_1+f_u_prime(x,w1(i),w2(i))*u_2);
        k21 = h*(u_2+0.5*k12);
        k22 = h*(f_u(x+0.5*h,w1(i),w2(i))*(u_1+0.5*k11)+f_u_prime(x+0.5*h,w1(i),w2(i))*(u_2+0.5*k21));
        k31 = h*(u_2+0.5*k22);
        k32 = h*(f_u(x+0.5*h,w1(i),w2(i))*(u_1+0.5*k21)+f_u_prime(x+0.5*h,w1(i),w2(i))*(u_2+0.5*k22));
        k41 = h*(u_2+k32);
        k42 = h*(f_u(x+h,w1(i),w2(i))*(u_1+k31)+f_u_prime(x+h,w1(i),w2(i))*(u_2+k32));
        u_1 = u_1+(k11+2*k21+2*k31+k41)/6;
        u_2 = u_2+(k12+2*k22+2*k32+k42)/6;
    end


    TK = TK-(w1(n+1)-beta)/u_1;
    k = k+1;
    x_list(n+1) = B;
    solplot = plot(x_list,w1);
    
   
end

