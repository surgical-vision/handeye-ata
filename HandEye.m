function [XL, XR] = HandEye(camera_pose, robot_pose, FLAG_STEREO, FLAG_ALGO, stereoT)
    N = size(camera_pose, 3);
    
    camRTcamL = stereoT;

    if (FLAG_STEREO)  %%Use stereo information (calibrate both camera depenently)%%
        
            AR = zeros(4, 4, 4*N);
            AL = zeros(4, 4, 4*N);
            B = zeros(4, 4, 4*N);
            for i = 1:N
                if (i == N)
                    AR(:, :, i) = camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 1);
                    AL(:, :, i) = camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 2);
                    AR(:, :, i + N) = camRTcamL*camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 2)/camRTcamL;
                    AL(:, :, i + N) = camRTcamL\camera_pose(:, :, 1, 1)/(camera_pose(:, :, i, 1))*camRTcamL;
                    AR(:, :, i + 2*N) = camRTcamL*camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 1);
                    AL(:, :, i + 2*N) = camRTcamL\camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 2);
                    AR(:, :, i + 3*N) = camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 2)/camRTcamL;
                    AL(:, :, i + 3*N) = camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 1)*camRTcamL;
                    B(:, :, i) = robot_pose(:, :, 1)\robot_pose(:, :, i);
                else
                    AR(:, :, i) = camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 1);
                    AL(:, :, i) = camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 2);
                    AR(:, :, i + N) = camRTcamL*camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 2)/camRTcamL;
                    AL(:, :, i + N) = camRTcamL\camera_pose(:, :, i + 1, 1)/(camera_pose(:, :, i, 1))*camRTcamL;
                    AR(:, :, i + 2*N) = camRTcamL*camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 1);
                    AL(:, :, i + 2*N) = camRTcamL\camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 2);
                    AR(:, :, i + 3*N) = camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 2)/camRTcamL;
                    AL(:, :, i + 3*N) = camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 1)*camRTcamL;
                    B(:, :, i) = robot_pose(:, :, i + 1)\robot_pose(:, :, i);
                end
            end
            B(:, :, N + 1:2*N) = B(:, :, 1:N);
            B(:, :, 2*N + 1:3*N) = B(:, :, 1:N);
            B(:, :, 3*N + 1:4*N) = B(:, :, 1:N);

        switch FLAG_ALGO
            case 1 %% Our algorithm %%
                XR = HandEye_ST(AR, B);
                XL = HandEye_ST(AL, B);
            case 2 %% Barreto's algorithm (Improved Dual Quaternion, linear) %%
                XR = HandEye_IDQ(AR, B);
                XL = HandEye_IDQ(AL, B);
            case 3 %% Tsai's algorithm (Conventional method) %%
                XR = HandEye_Tsai(AR, B);
                XL = HandEye_Tsai(AL, B);
            case 4 %% Konstantinos' algorithm (Dual Quaternion with quadratic eq.) %%
                XR = HandEye_DQ(AR, B);
                XL = HandEye_DQ(AL, B);
            otherwise
                error('Invalid ALGORITHM FLAG.');
        end
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
        switch FLAG_ALGO
            case 1 %% Our algorithm %%
                XR = HandEye_ST(A, B);
            case 2 %% Barreto's algorithm (Improved Dual Quaternion, linear) %%
                XR = HandEye_IDQ(A, B);
            case 3 %% Tsai's algorithm (Conventional method) %%
                XR = HandEye_Tsai(A, B);
            case 4 %% Konstantinos' algorithm (Dual Quaternion with quadratic eq.) %%
                XR = HandEye_DQ(A, B);
            otherwise
                error('Invalid ALGORITHM FLAG.');
        end
        XL = stereoT\XR;
    end

%     The non-linear optimisation is for real data where a lot of local
%     minima can be found.
%     h = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 30000);
%         
%     
%     x_init = [rodrigues(XR(1:3, 1:3)); XR(1:3, 4)];
%     [refine_X, ~, ~, ~, ~, ~, J] = lsqnonlin(@(x) optimAXXB(x, AR, B), x_init, [], [], h);
%     
%     XR = [rodrigues(refine_X(1:3)), refine_X(4:6);0 0 0 1];
%     
%     x_init = [rodrigues(XL(1:3, 1:3)); XL(1:3, 4)];
%     [refine_X, ~, ~, ~, ~, ~, J] = lsqnonlin(@(x) optimAXXB(x, AL, B), x_init, [], [], h);
%     
%     XL = [rodrigues(refine_X(1:3)), refine_X(4:6);0 0 0 1];
end