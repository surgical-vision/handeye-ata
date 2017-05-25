function X = HandEye_ST(A, B)
    converge_thres = 10;
    counter_rotation = 0;
    counter_translation = 0;
    max_ite = 50;
    
    N = size(A, 3);

    v_a = zeros(3, N);
    v_b = zeros(3, N);
    om_a = zeros(3, N);
    om_b = zeros(3, N);
    K = zeros(8*N, 4);

    X = eye(4);
    
    LHS = zeros(3*N, 3);
    RHS = zeros(3*N, 1);
    
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

        K(8*j - 7:8*j - 4, 1:4) = [x(4), -(x(1:3))'; x(1:3), skew3(y) + x(4)*eye(3)];
    end

    counter = 0;
    
    while (((counter_rotation < converge_thres) || (counter_translation < converge_thres)) && (counter < max_ite))
        
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
            LHS(3*j - 2:3*j, :) = skew3(om_a(:, j));
            RHS(3*j - 2:3*j, 1) = X(1:3, 1:3)*v_b(:, j) - v_a(:, j);
        end


        X = [X(1:3, 1:3) LHS\RHS;0 0 0 1];
        
        diff = X/[R_init, t_init;0 0 0 1];
        
        if (norm(rodrigues(diff(1:3, 1:3))) < 1e-3)
            counter_rotation = counter_rotation + 1;
        else
            counter_rotation = 0;
        end
        
        if (norm(diff(1:3, 4)) < 1e-3)
            counter_translation = counter_translation + 1;
        else
            counter_translation = 0;
        end
        
        counter = counter + 1;
        ccc = rodrigues(X(1:3, 1:3));
        C(1:3, counter) = ccc;
        C(4:6, counter) = X(1:3, 4);
    end
        
end