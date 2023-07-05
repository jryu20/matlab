function B = HTriDiag(A)
% ==========================================
% This function implements an efficient way to tridiagonalize A via
% repeated conjugation with Householder transformations.
% ==========================================

[m,n] = size(A);
B = zeros(n,n);

for i=1:n-1
    B(i,i) = A(1,1);
    atilde = A(2:end,1);
    vtemp = zeros(size(atilde));
    vtemp(1) = norm(atilde,2);
    B(i,i+1) = norm(atilde,2);
    B(i+1,i) = norm(atilde,2);
    v = atilde - vtemp;
    A = HMatMult(v,HMatMult(v,A(2:end,2:end)')');
end

B(n,n) = A(1,1);

end