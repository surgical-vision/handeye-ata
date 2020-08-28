function [XL, XR, C] = HandEye_test_convergence(camera_pose, robot_pose, FLAG_STEREO, FLAG_ALGO, stereoT)

    N = size(robot_pose, 3);
    if (FLAG_STEREO)  %%Use stereo information (calibrate both camera depenently)%%

        A = zeros(4, 4, 4*N);
        B = zeros(4, 4, 4*N);
        for i = 1:N
            if (i == N)
                A(:, :, i) = camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 1);
                A(:, :, i + N) = stereoT*camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 2)/stereoT;
                A(:, :, i + 2*N) = stereoT*camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 1);
                A(:, :, i + 3*N) = camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 2)/stereoT;
                B(:, :, i) = robot_pose(:, :, 1)\robot_pose(:, :, i);
            else
                A(:, :, i) = camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 1);
                A(:, :, i + N) = stereoT*camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 2)/stereoT;
                A(:, :, i + 2*N) = stereoT*camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 1);
                A(:, :, i + 3*N) = camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 2)/stereoT;
                B(:, :, i) = robot_pose(:, :, i + 1)\robot_pose(:, :, i);
            end
        end
        B(:, :, N + 1:2*N) = B(:, :, 1:N);
        B(:, :, 2*N + 1:3*N) = B(:, :, 1:N);
        B(:, :, 3*N + 1:4*N) = B(:, :, 1:N);

    else %%Not Use stereo information%%

        A = zeros(4, 4, N);
        B = zeros(4, 4, N);
        for i = 1:N
            if (i == N)
                A(:, :, i) = camera_pose(:, :, 1)/camera_pose(:, :, i);
                B(:, :, i) = robot_pose(:, :, 1)\robot_pose(:, :, i);
            else
                A(:, :, i) = camera_pose(:, :, i + 1)/camera_pose(:, :, i);
                B(:, :, i) = robot_pose(:, :, i + 1)\robot_pose(:, :, i);
            end
        end

    end

    [XR, C] = HandEye_ST_test_convergence(A, B, FLAG_ALGO);

    if (size(camera_pose, 4) == 1)
        XL = eye(4);
    else
        XL = stereoT\XR;
    end
end