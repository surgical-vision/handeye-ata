function X = HandEye_IDQ(A, B)
    %Estimate Camera to robot gripper transformation from the equation 
    %Ta: Camera to grid transformation in the form of (4x4xN)
    %Tb: Gripper to robot base transformation in the form on (4x4xN)
    %, where N is the number of captured movement.
    %AX = XB
    N = size(A, 3);
    %Construct L and Lp matrix
    L = zeros(4*N, 4);
    Lp = zeros(4*N, 4);
    for j = 1:N
        [a, ap] = getDualQ(A(1:3, 1:3, j), A(1:3, 4, j));
        [b, bp] = getDualQ(B(1:3, 1:3, j), B(1:3, 4, j));
        x = a - b;
        y = a + b;
        z = ap - bp;
        w = ap + bp;
        L(4*j-3:4*j, :) = [x(4), -x(1:3)';x(1:3), skew3(y(1:3))+x(4)*eye(3)];
        Lp(4*j-3:4*j, :) = [z(4), -z(1:3)';z(1:3), skew3(w(1:3))+z(4)*eye(3)];
    end

    %Determine real part and dual part of matrix X
    [~, ~, V] = svd(L);
    q = V(:, 4);
    qprime = L\(-Lp*q);

    q = [q(2:4); q(1)];
    qprime = [qprime(2:4); qprime(1)];
    
    %Determine eye to hand tranformation
    R = q2dcm(q)';
    t = 2*qmult(qprime,qconj(q));
    t = t(1:3);
    X = [R, t;0 0 0 1];
end