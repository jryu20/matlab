function [solutions] = solveLSQR(A,varargin)

    if size(A,1) <= size(A,2)
        error("Matrix dimensions do not match.")
    end


    [Q,R] = findQR(A);
    R1 = R(1:size(A,2),:);
    varargout = cell(1,nargin-1);
    for i = 1:(nargin-1)
        
        if size(varargin{i}) ~= size(A,1)
        error("Matrix dimensions do not match.")
        end

        c = (Q.')*varargin{i};
        c1 = c(1:size(A,2));

        x = zeros(size(c1,1),1);
        for j = size(c1,1):-1:1   
            x(j) = c1(j)/R1(j,j);
            c1(1:j-1) = c1(1:j-1)-(R1(1:j-1,j)*x(j));
            x_k = x;
        end

        varargout{i} = x_k;
        solutions = varargout;

    end
 end
 
