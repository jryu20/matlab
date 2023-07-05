function [solplot] = finite_diff(a,b,alpha,beta,h)
%uses finite difference method to solve a linear BVP
%the code for Jacobi's method has been mainly pulled from: 
%www.mathworks.com/matlabcentral/fileexchange/72563-gauss-jacobi-method


n = (b-a)/h;
x = zeros(n-1,1);
for i=1:n-1
    x(i) = a+i*h;
end

p = 2;
q = -1;
r = @(x) x*(exp(1)^x)-x;

A = zeros(n-1,n-1); 
for i=1:n-1
    A(i,i) = 2+(h^2)*polyval(q,x(i));
end
for i=2:n-1
    A(i,i-1) = -1-(h*polyval(p,x(i))/2);
    A(i-1,i) = -1+(h*polyval(p,x(i-1))/2);
end


B = zeros(n-1,1);
B(1) = -(h^2)*r(x(1))+(1+(h*polyval(p,x(1))/2))*alpha;  
B(n-1) = -(h^2)*r(x(n-1))+(1-(h*polyval(p,x(n-1))/2))*beta;
for i=2:n-2
    B(i) = -(h^2)*r(x(i));

end

P=[A B]; 
[row, col] = size(P); 
U=zeros(row,1); 
C=zeros(row,1);
Err=ones(row,1); 

merr=max(Err);
while merr>0.01
    for m=1:1:row       
       C(m,1)=(1/P(m,m))*(P(m,col)-sum(A(m,:)*U(:,1))+A(m,m)*U(m,1));
       Err(m,1)= abs(C(m,1)-U(m,1));
    end 
    U(:,1)=C(:,1);
    merr=max(Err);
end

solplot = plot(x,U);