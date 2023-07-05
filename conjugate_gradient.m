function [xk] = conjugate_gradient(A, b, K, x)
%implements the conjugate gradient algorithm with K iterations

    if size(A,1) < size(A,2) || size(A,1) ~= size(b,1) 
        error("Matrix dimensions do not match.")
    end
    
    Q_A = (A.')*A;
    rk = Q_A*x - (A.')*b;
    dk = rk;
    xk = x;

    for i = 0:K
        alpha_k = ((rk.')*dk)/((dk.')*Q_A*dk);
        xk = xk-(alpha_k*dk);
        rk = Q_A*xk - (A.')*b;
        beta_k = ((rk.')*Q_A*dk)/((dk.')*Q_A*dk);
        dk = rk-(beta_k*dk);
    end
 end  
   
    