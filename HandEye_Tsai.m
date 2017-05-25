function X = HandEye_Tsai(A, B)
    %Estimate Camera to robot gripper transformation from the equation 
    %TaL: Left camera to grid transformation in the form of (4x4xN)
    %Tb: Gripper to robot base transformation in the form of (4x4xN)
    %, where N is the number of captured movement.
    X = eye(4);
    N = size(A, 3);
    RLHS = zeros(3*N, 3);
    RRHS = zeros(3*N, 1);
    for j = 1:N
        rA = rodrigues(A(1:3, 1:3, j));
        rB = rodrigues(B(1:3, 1:3, j));
        thetaA = norm(rA);
        thetaB = norm(rB);
        rAn = rA/thetaA;
        rBn = rB/thetaB;
        pA = 2*sin(thetaA/2)*rAn;
        pB = 2*sin(thetaB/2)*rBn;
        RLHS(3*j - 2:3*j, :) = skew3(pA + pB);
        RRHS(3*j - 2:3*j, :) = pB - pA;
    end
    
    axisRR = RLHS\RRHS;
    axisRR = 2*axisRR/sqrt(1 + norm(axisRR)^2);
    X(1:3, 1:3) = (1-norm(axisRR)*norm(axisRR)/2)*eye(3)+0.5*(axisRR*axisRR'+sqrt(4-norm(axisRR)*norm(axisRR))*skew3(axisRR));
    LHS = zeros(3*N, 3);
    RHS = zeros(3*N, 1);
    for j = 1:N
        LHS(3*j - 2:3*j, :) = A(1:3, 1:3, j) - eye(3);
        RHS(3*j - 2:3*j, :) = X(1:3, 1:3)*B(1:3, 4, j) - A(1:3, 4, j);
    end
    X(1:3, 4) = LHS\RHS;
end