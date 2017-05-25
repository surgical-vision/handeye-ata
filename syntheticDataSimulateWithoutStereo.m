%%%%SYNTHETIC DATA EXPERIMENT
%%%A comparison between hand-eye calibration with stereo information and 
%%%without, checking whether it is improved with stereo information.

clc;
clear all;
close all;

R_w = [0 0 0]'; %Robot origin
G_w = [-0.3, 0.5,  0.7]'; %Calibration grid position

N_motion = 6;

using_stereo = 1;
using_motion_selection = 0;

XL = cell(4, 1);
XR = cell(4, 1);
XL_S = cell(4, 1);
XR_S = cell(4, 1);

%Noise coefficient
NoiseCoeffT = 0:2e-4:2e-3;
NoiseCoeffR = (0:0.1:1)*pi/180;
%Initialize error matrix
DegreeErrorLR = zeros(2, length(NoiseCoeffT), 5);
DegreeSTDLR = zeros(2, length(NoiseCoeffT), 5);
TranErrorLR = zeros(2, length(NoiseCoeffT), 5);
TranSTDLR = zeros(2, length(NoiseCoeffT), 5);
DegreeErrorLR_S = zeros(2, length(NoiseCoeffT), 5);
DegreeSTDLR_S = zeros(2, length(NoiseCoeffT), 5);
TranErrorLR_S = zeros(2, length(NoiseCoeffT), 5);
TranSTDLR_S = zeros(2, length(NoiseCoeffT), 5);



DQ_produce_imag = zeros(2, length(NoiseCoeffT));

%State how many times the program is simulated per one noise coefficient
nTimes = 100;



%Number of Point
for noise = 1:length(NoiseCoeffT)
    DegreeErrorL = zeros(1, nTimes);
    TranErrorL = zeros(1, nTimes);
    DegreeErrorR = zeros(1, nTimes);
    TranErrorR = zeros(1, nTimes);
    
    DegreeErrorL_S = zeros(1, nTimes);
    TranErrorL_S = zeros(1, nTimes);
    DegreeErrorR_S = zeros(1, nTimes);
    TranErrorR_S = zeros(1, nTimes);
    
    
    
    for j = 1:nTimes
        
        %Setting offset respected to end effector
        %(Transformation from robot to camera)
        X_off = 0.2*randn(1);
        Y_off = 0.27*randn(1);
        Z_off = randn(1);


        %Create eye to hand transformation.
        theta = 6*pi/180*randn(1);
        alpha = 6*pi/180*randn(1);
        phi = 6*pi/180*randn(1);
        Rx = Rot(theta, 'y')*Rot(phi, 'x')*Rot(alpha, 'z');
        Tgripper2camL = [eye(3), [X_off; Y_off; Z_off].*randn(3,1); 0 0 0 1]*[Rx, zeros(3, 1);0 0 0 1];

        x_l2r = 0.1*randn(1);
        y_l2r = 0.05*randn(1);
        z_l2r = 0.5*randn(1);

        theta = 6*pi/180*randn(1);
        alpha = 6*pi/180*randn(1);
        phi = 6*pi/180*randn(1);
        Rx = Rot(theta, 'y')*Rot(phi, 'x')*Rot(alpha, 'z');
        Tleft2right = [eye(3), [x_l2r; y_l2r; z_l2r].*randn(3,1); 0 0 0 1]*[Rx, zeros(3, 1);0 0 0 1];
        Tgripper2camR = Tleft2right*Tgripper2camL;
        
        %Calculate robot base to calibration grid transformation
        Tbase2grid = [1 0 0 G_w(1);0 1 0 G_w(2);0 0 1 G_w(3);0 0 0 1]*[Rot(-pi/2, 'z'), zeros(3, 1);0 0 0 1];
        Tgrid2base = inv(Tbase2grid);
        
        fprintf('noise number = %d,, trial number = %d\n\n', noise, j);         

        %Initialize gripper transformation and camera transformation
        Tgrid2camL = zeros(4, 4, N_motion);
        Tgrid2camR = zeros(4, 4, N_motion);
        Tgripper2base = zeros(4, 4, N_motion);
        Tbase2gripper = zeros(4, 4, N_motion);

        ranT = 0.5*randn(3, N_motion);
        ranR = 15*pi*randn(3, N_motion)/180;
        noiseT = NoiseCoeffT(noise)*randn(3, N_motion);
        noiseR = NoiseCoeffR(noise)*randn(3, N_motion);
        noiseTcamL = NoiseCoeffT(noise)*randn(3, N_motion);
        noiseRcamL = NoiseCoeffR(noise)*randn(3, N_motion);
        noiseTcamR = NoiseCoeffT(noise)*randn(3, N_motion);
        noiseRcamR = NoiseCoeffR(noise)*randn(3, N_motion);
        
        Tgripper2baseNoise = zeros(size(Tbase2gripper)); %robot 2 world (from forward kin)
        Tbase2gripperNoise = zeros(size(Tbase2gripper)); %world 2 robot
        Tgrid2camLNoise = zeros(size(Tgrid2camL));
        Tgrid2camRNoise = zeros(size(Tgrid2camR));

        for i = 1:N_motion     
            Tgripper2base(:, :, i) = [rodrigues(ranR(:, i)), ranT(:, i);0 0 0 1];

            Tbase2gripper(:, :, i) = inv(Tgripper2base(:, :, i));

            Tgrid2camL(:, :, i) = Tgripper2camL*Tbase2gripper(:, :, i)*Tgrid2base;
            Tgrid2camR(:, :, i) = Tgripper2camR*Tbase2gripper(:, :, i)*Tgrid2base;

            Tgripper2baseNoise(:, :, i) = Tgripper2base(:, :, i)*[rodrigues(noiseR(:, i)), noiseT(:, i);0 0 0 1];
            Tbase2gripperNoise(:, :, i) = inv(Tgripper2baseNoise(:, :, i));
            Tgrid2camLNoise(:, :, i) = Tgrid2camL(:, :, i)*[rodrigues(noiseRcamL(:, i)), noiseTcamL(:, i);0 0 0 1];
            Tgrid2camRNoise(:, :, i) = Tgrid2camR(:, :, i)*[rodrigues(noiseRcamR(:, i)), noiseTcamL(:, i);0 0 0 1];
        end
        
        TA(:, :, :, 1) = Tgrid2camRNoise;
        TA(:, :, :, 2) = Tgrid2camLNoise;
        TB = Tgripper2baseNoise;
        
        for i = 1:4
            [XL, XR] = HandEye(TA, TB, ~using_stereo, i, Tleft2right);
            [XL_S, XR_S] = HandEye(TA, TB, using_stereo, i, Tleft2right);
            
            if ((i == 4) && (~isreal(XL)) || (~isreal(XR)))
                DQ_produce_imag(noise) = DQ_produce_imag(noise) + 1;
                XL = Tgripper2camL;
                XR = Tgripper2camR;
            end
            
            if ((i == 4) && (~isreal(XL_S)) || (~isreal(XR_S)))
                DQ_produce_imag(noise) = DQ_produce_imag(noise) + 1;
                XL = Tgripper2camL;
                XR = Tgripper2camR;
            end
            
            dL = Tgripper2camL/XL;
            dR = Tgripper2camR/XR;
                
            DegreeErrorL(i, j) = norm(rodrigues(dL(1:3, 1:3)))*180/pi;
            DegreeErrorR(i, j) = norm(rodrigues(dR(1:3, 1:3)))*180/pi;
            TranErrorL(i, j) = norm(dL(1:3, 4))*1000;
            TranErrorR(i, j) = norm(dR(1:3, 4))*1000;
            
            dL = Tgripper2camL/XL_S;
            dR = Tgripper2camR/XR_S;
            
            DegreeErrorL_S(i, j) = norm(rodrigues(dL(1:3, 1:3)))*180/pi;
            DegreeErrorR_S(i, j) = norm(rodrigues(dR(1:3, 1:3)))*180/pi;
            TranErrorL_S(i, j) = norm(dL(1:3, 4))*1000;
            TranErrorR_S(i, j) = norm(dR(1:3, 4))*1000;
        end
        
    end
    
    DegreeError_cell = cell(2, 4);
    TranError_cell = cell(2, 4);
    
    DegreeError_cell_S = cell(2, 4);
    TranError_cell_S = cell(2, 4);
    
    %Record the data
    for i = 1:4
        
        %Delete outlier results
        DegreeError_cell{1, i} = deleteoutliers(DegreeErrorL(i, :));
        DegreeError_cell{2, i} = deleteoutliers(DegreeErrorR(i, :));
        TranError_cell{1, i} = deleteoutliers(TranErrorL(i, :));
        TranError_cell{2, i} = deleteoutliers(TranErrorR(i, :));
        
        %Delete outlier results
        DegreeError_cell_S{1, i} = deleteoutliers(DegreeErrorL_S(i, :));
        DegreeError_cell_S{2, i} = deleteoutliers(DegreeErrorR_S(i, :));
        TranError_cell_S{1, i} = deleteoutliers(TranErrorL_S(i, :));
        TranError_cell_S{2, i} = deleteoutliers(TranErrorR_S(i, :));
        
        DegreeErrorLR(1, noise, i) = mean(DegreeError_cell{1, i});
        DegreeSTDLR(1, noise, i) = std(DegreeError_cell{1, i});

        TranErrorLR(1, noise, i) = mean(TranError_cell{1, i});
        TranSTDLR(1, noise, i) = std(TranError_cell{1, i});

        DegreeErrorLR(2, noise, i) = mean(DegreeError_cell{2, i});
        DegreeSTDLR(2, noise, i) = std(DegreeError_cell{2, i});

        TranErrorLR(2, noise, i) = mean(TranError_cell{2, i});
        TranSTDLR(2, noise, i) = std(TranError_cell{2, i});
        
        DegreeErrorLR_S(1, noise, i) = mean(DegreeError_cell_S{1, i});
        DegreeSTDLR_S(1, noise, i) = std(DegreeError_cell_S{1, i});

        TranErrorLR_S(1, noise, i) = mean(TranError_cell_S{1, i});
        TranSTDLR_S(1, noise, i) = std(TranError_cell_S{1, i});

        DegreeErrorLR_S(2, noise, i) = mean(DegreeError_cell_S{2, i});
        DegreeSTDLR(2, noise, i) = std(DegreeError_cell_S{2, i});

        TranErrorLR_S(2, noise, i) = mean(TranError_cell_S{2, i});
        TranSTDLR_S(2, noise, i) = std(TranError_cell_S{2, i});
    end
end

for i = 1:4
    figure(i), 
    subplot(2, 1, 1)
    errorbar(NoiseCoeffT*1000, TranErrorLR(1, :, i), TranSTDLR(1, :, i), 'r'), hold on, 
    errorbar(NoiseCoeffT*1000, TranErrorLR_S(1, :, i), TranSTDLR_S(1, :, i), 'g'), 
    xlabel('noise coeff (mm)'), ylabel('error (mm)'),
    title('translational error')
    legend('without stereo', 'with stereo'), grid on, 
    hold off
    subplot(2, 1, 2)
    errorbar(NoiseCoeffR*180/pi, DegreeErrorLR(1, :, i), DegreeSTDLR(1, :, i), 'r'), hold on,
    errorbar(NoiseCoeffR*180/pi, DegreeErrorLR_S(1, :, i), DegreeSTDLR_S(1, :, i), 'g'), 
    xlabel('noise coeff (degree)'), ylabel('error (degree)')
    title('rotational error')
    legend('without stereo', 'with stereo'), grid on, 
    hold off
end