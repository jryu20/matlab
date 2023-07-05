function [v, lambda] = normalized_power(A,K)
% implements the normalized power method
%   this finds the leading eigenvalue and eigenvector of a square matrix

if size(A,1) ~= size(A,2)
    error("Your matrix A must be square")
end


x_0 = rand(size(A,1),1);
x_k = x_0;

for k = 1:K
    y_k = A*x_k;
    x_k = (1/norm(y_k,2))*y_k;
    lambda = ((x_k')*A*x_k)/((x_k')*x_k);

    v = x_k;

end