function [X, C] = HandEye_ST_test_convergence(A, B, FLAG_ALGO)

    max_ite = 20;

    C = zeros(6, 31);
    
    N = size(A, 3);

    v_a = zeros(3, N);
    v_b = zeros(3, N);
    om_a = zeros(3, N);
    om_b = zeros(3, N);
    K = zeros(8*N, 4);
    switch FLAG_ALGO
        case 1
            X = HandEye_IDQ(A, B);
        case 2
            X = HandEye_Tsai(A, B);
        case 3
            X = HandEye_DQ(A, B);
        case 4
            X = eye(4);
        otherwise
            error('Invalid input.\n\n');
    end
    
    
    LHS = zeros(3*N, 3);
    RHS = zeros(3*N, 1);
    
    plane = zeros(N, 3);
    %Construct L and Lp matrix
    L = zeros(4*N, 4);
    Lp = zeros(4*N, 4);
    
    for j = 1:N
        a = logm(A(:, :, j));
        b = logm(B(:, :, j));
        v_a(:, j) = a(1:3 ,4);
        v_b(:, j) = b(1:3, 4);
        om_a(:, j) = rodrigues(A(1:3, 1:3, j));
        om_b(:, j) = rodrigues(B(1:3, 1:3, j));
        LHS(3*j - 2:3*j, :) = skew3(om_a(:, j));
        [a, ap] = getDualQ(A(1:3, 1:3, j), A(1:3, 4, j));
        [b, bp] = getDualQ(B(1:3, 1:3, j), B(1:3, 4, j));
        x = a - b;
        y = a + b;
        z = ap - bp;
        w = ap + bp;
        L(4*j-3:4*j, :) = [x(4), -x(1:3)';x(1:3), skew3(y(1:3))+x(4)*eye(3)];
        Lp(4*j-3:4*j, :) = [z(4), -z(1:3)';z(1:3), skew3(w(1:3))+z(4)*eye(3)];

        K(8*j - 7:8*j - 4, 1:4) = [x(4), -(x(1:3))'; x(1:3), skew3(y) + x(4)*eye(3)];

    end

    
    
    counter = 1;
    
    C(1:3, counter) = rodrigues(X(1:3, 1:3));
    C(4:6, counter) = X(1:3, 4);
    
    while (counter <= max_ite)
        
        R_init = X(1:3, 1:3);
        t_init = X(1:3, 4);
        
        for j = 1:N
            vec_buf = v_a(:, j) - cross(t_init, om_a(:, j));

            K(8*j - 3:8*j, 1) = [vec_buf - v_b(:, j); 0];
            K(8*j - 3:8*j, 2:4) = [skew3(vec_buf + v_b(:, j)); (-vec_buf + v_b(:, j))'];
        end
        
        [~, ~, v_basis] = svd(K);
        v_basis = v_basis(:, 4);

        qR = [v_basis(2:4); v_basis(1)];
        
        X(1:3, 1:3) = q2dcm(qR)';

        for j = 1:N
            RHS(3*j - 2:3*j, 1) = X(1:3, 1:3)*v_b(:, j) - v_a(:, j);
        end


        X = [X(1:3, 1:3) LHS\RHS;0 0 0 1];
        
        diff = X/[R_init, t_init;0 0 0 1];
        

        C(1:3, counter + 1) = rodrigues(X(1:3, 1:3));
        C(4:6, counter + 1) = X(1:3, 4);

        counter = counter + 1;
        
        
    end

end