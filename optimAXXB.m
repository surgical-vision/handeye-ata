function E = optimAXXB(x, A, B)
    N = size(A, 3);
    X = [rodrigues(x(1:3)), x(4:6);0 0 0 1];
    E = 0;    
    for i = 1:N
        E = E + norm((A(:, :, i)*X)\X*B(:, :, i) - eye(4));
    end
end