function [Y, t] = RK42d(f, t0, T, y0, N)

h = (T - t0)/(N - 1); %Calculate and store the step size
Y = zeros(N,2); %Initialize the X and Y vector
t = linspace(t0,T,N); % A vector to store the time values
Y(1,:) = y0; % Start Y vector at the intial values.
for i = 1:(N-1)
    y = Y(i,:)';
    k1 = f(t(i),y);
    k2 = f(t(i) +0.5*h, y +0.5*h*k1);
    k3 = f(t(i) +0.5*h, y +0.5*h*k2);
    k4 = f(t(i) +h, y +h*k3);
    Y(i+1,:) = y + (h/6)*(k1+ 2.*k2 + 2*k3 + k4);
end
end