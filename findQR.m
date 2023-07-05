function [Q,R] = findQR(A)
%finds the QR decomposition of a matrix using Householder transformations

    if size(A,1) < size(A,2)
        error("Matrix dimensions do not match.")
    end
    
    n = size(A,1);
    m = size(A,2);
    P_matrix_cell = cell(m,1);
    a_new = A;
    for i = 1:m
        switch i
            case 1
                a = a_new(1:end,1:end);
            otherwise
                a = a_new(2:end,2:end);
        end

        I = eye(size(a,1));
        e1 = I(:,1); 
        a_i = a(:,1);
        v = a_i - norm(a_i)*e1;
        u = (1/norm(v))*v;
        P_matrix = I - 2*(u)*(u.');

           
        a_new = P_matrix*a;

        if (i > 1)
            matrix_n = zeros(n);
            matrix_n(i:end,i:end) = P_matrix;
            P_matrix = matrix_n;
            for j = 1:i-1
                P_matrix(j,j) = 1;
            end  

        end

        P_matrix_cell{i} = P_matrix;
    end
    
    Q_1 = P_matrix_cell{1};
    for k = 2:m
        Q_1 = Q_1 * P_matrix_cell{k};
        Q = Q_1;
    end
    
    R = Q^(-1)*A;

end