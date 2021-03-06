%%%%SYNTHETIC DATA EXPERIMENT
%%%Test the performance of each algorithm when they are tested with
%%%an increasing number of motion included in the calibration.

clc;
clear;

%Generate grid
points.gt = zeros(4, 196);
points.gt(4, :) = 1;
points.gt(1, :) = repmat(0:3:39, [1, 14]);
for i = 1:14
    points.gt(2, 14*i-13:14*i) = 3*(i - 1);
end

%Generate camera data
K.fcL = [800, 900]';
K.fcR = [800, 900]';
%Intrinsic matrix for left and right camera
K.K0L = [850, 0, 0; 0, 920, 0;0, 0, 1];
K.K0R = [850, 0, 0; 0, 920, 0;0, 0, 1];
%Distortion parameters
K.kcbpL = zeros(5, 1);
K.kcbpR = zeros(5, 1);
%Camera offset
K.ccL = [400, 200]';
K.ccR = [400, 200]';
KMat = [850, 0, 400, 0; 0, 920, 200, 0;0 0 1, 0];

%Initialise optim options
h = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 1500);

R_w = [0 0 0]'; %Robot origin
G_w = [-0.3, 0.5,  0.7]'; %Calibration grid position

% 0 for using Algorithm 1 (monocular)
% 1 for using Algorithm 2 (stereo)
using_stereo = 1;

%State how many times the program is simulated per one noise coefficient
nTimes = 100;

%Noise coefficient
NoiseCoeffT = 5e-4;
NoiseCoeffR = pi/180;

% Number of motion included in calibration
N_motion = 3:12;

%Initialize error matrix
DegreeErrorLR = zeros(2, 4, length(N_motion));
DegreeSTDLR = zeros(2, 4, length(N_motion));
TranErrorLR = zeros(2, 4, length(N_motion));
TranSTDLR = zeros(2, 4, length(N_motion));

%Run for every included number of motion
for N = 1:length(N_motion)
    
    DegreeErrorL = zeros(5, nTimes);
    TranErrorL = zeros(5, nTimes);
    DegreeErrorR = zeros(5, nTimes);
    TranErrorR = zeros(5, nTimes);
    
    %Repeat the experiment for (nTimes) times.
    for j = 1:nTimes

        %Calculate robot base to calibration grid transformation
        Tbase2grid = [1 0 0 G_w(1)*randn(1, 1);0 1 0 G_w(2)*randn(1, 1);0 0 1 G_w(3)*randn(1, 1);0 0 0 1]*[Rot(-pi/2, 'z'), zeros(3, 1);0 0 0 1];
        Tgrid2base = inv(Tbase2grid);

        %Setting offset respected to end effector
        %(Transformation from robot to camera)
        X_off = -0.5*randn(1);
        Y_off = 0.27*randn(1);
        Z_off = randn(1);

        %Create eye to hand transformation.
        theta = 5*pi/180*randn(1);
        alpha = 6*pi/180*randn(1);
        phi = 6*pi/180*randn(1);
        Rx = Rot(theta, 'y')*Rot(phi, 'x')*Rot(alpha, 'z');
        %The hand-eye transformation for left camera: XL
        Tgripper2camL = [eye(3), [X_off; Y_off; Z_off].*randn(3,1); 0 0 0 1]*[Rx, zeros(3, 1);0 0 0 1];

        theta = 6*pi/180*randn(1);
        alpha = 4*pi/180*randn(1);
        phi = 3*pi/180*randn(1);
        Rx = Rot(theta, 'y')*Rot(phi, 'x')*Rot(alpha, 'z');
        %The transformation between left and right camera
        Tleft2right = [eye(3), [X_off; Y_off; Z_off].*randn(3,1); 0 0 0 1]*[Rx, zeros(3, 1);0 0 0 1];
        %The hand-eye transformation for left camera: XR
        Tgripper2camR = Tleft2right*Tgripper2camL;

        fprintf('motion number = %d ,, trial number = %d\n\n', N, j);         

        %Initialize gripper transformation and camera transformation
        Tgrid2camL = zeros(4, 4, N_motion(N) );
        Tgrid2camR = zeros(4, 4, N_motion(N) );
        Tgripper2base = zeros(4, 4, N_motion(N) );
        Tbase2gripper = zeros(4, 4, N_motion(N) );

        ranT = 0.5*randn(3, N_motion(N) );
        ranR = 30*pi*randn(3, N_motion(N) )/180;
        noiseT = NoiseCoeffT*randn(3, N_motion(N) );
        noiseR = NoiseCoeffR*randn(3, N_motion(N) );
        noiseTcamL = NoiseCoeffT*randn(3, N_motion(N) );
        noiseRcamL = NoiseCoeffR*randn(3, N_motion(N) );
        noiseTcamR = NoiseCoeffT*randn(3, N_motion(N) );
        noiseRcamR = NoiseCoeffR*randn(3, N_motion(N) );

        Tgripper2baseNoise = zeros(size(Tbase2gripper)); %robot 2 world (from forward kin)
        Tbase2gripperNoise = zeros(size(Tbase2gripper)); %world 2 robot
        Tgrid2camLNoise = zeros(size(Tgrid2camL));
        Tgrid2camRNoise = zeros(size(Tgrid2camR));

        %Generate robot and camera motions (with and without noise)
        for i = 1:N_motion(N)
            %Robot and camera motions without noise
            Tgripper2base(:, :, i) = [rodrigues(ranR(:, i)), ranT(:, i);0 0 0 1];
            Tbase2gripper(:, :, i) = inv(Tgripper2base(:, :, i));
            Tgrid2camL(:, :, i) = Tgripper2camL*Tbase2gripper(:, :, i)*Tgrid2base;
            Tgrid2camR(:, :, i) = Tgripper2camR*Tbase2gripper(:, :, i)*Tgrid2base;

            %Robot and camera motions with noise
            Tgripper2baseNoise(:, :, i) = Tgripper2base(:, :, i)*[rodrigues(noiseR(:, i)), noiseT(:, i);0 0 0 1];
            Tbase2gripperNoise(:, :, i) = inv(Tgripper2baseNoise(:, :, i));
            Tgrid2camLNoise(:, :, i) = Tgrid2camL(:, :, i)*[rodrigues(noiseRcamL(:, i)), noiseTcamL(:, i);0 0 0 1];
            Tgrid2camRNoise(:, :, i) = Tgrid2camR(:, :, i)*[rodrigues(noiseRcamR(:, i)), noiseTcamL(:, i);0 0 0 1];
            
            %Generate points on images basing on camera pose
            points.points_on_images(:, :, i, 1) = KMat*Tgrid2camL(:, :, i)*points.gt;
            points.points_on_images(:, :, i, 2) = KMat*Tgrid2camR(:, :, i)*points.gt;
            points.points_on_images(:, :, i, 1) = points.points_on_images(:, :, i, 1)./repmat(points.points_on_images(3, :, i, 1), [3, 1, 1, 1]);
            points.points_on_images(:, :, i, 2) = points.points_on_images(:, :, i, 2)./repmat(points.points_on_images(3, :, i, 2), [3, 1, 1, 1]);
        end
        
        clear TA TB;
        
        TA(:, :, :, 1) = Tgrid2camRNoise;
        TA(:, :, :, 2) = Tgrid2camLNoise;
        TB = Tgripper2baseNoise;

        %Run hand-eye calibration for every calibration methods
        for i = 1:4

            [XL, XR] = HandEye(TA, TB, using_stereo, i, Tleft2right);
            
            XR = Tleft2right*XL;
            
                                    
            %Calculate error for both left and right hand-eye matrix
            dL = Tgripper2camL/XL;
            dR = Tgripper2camR/XR;

            DegreeErrorL(i, j) = norm(rodrigues(dL(1:3, 1:3)))*180/pi;
            DegreeErrorR(i, j) = norm(rodrigues(dR(1:3, 1:3)))*180/pi;
            TranErrorL(i, j) = norm(Tgripper2camL(1:3, 4) - XL(1:3, 4))*1000;
            TranErrorR(i, j) = norm(Tgripper2camR(1:3, 4) - XR(1:3, 4))*1000;
            
        end
    end
    
    DegreeError_cell = cell(2, 4);
    TranError_cell = cell(2, 4);
    
    %Record the data
    for i = 1:4
        %Delete outlier results
        DegreeError_cell{1, i} = deleteoutliers(DegreeErrorL(i, :));
        DegreeError_cell{2, i} = deleteoutliers(DegreeErrorR(i, :));
        TranError_cell{1, i} = deleteoutliers(TranErrorL(i, :));
        TranError_cell{2, i} = deleteoutliers(TranErrorR(i, :));
        
        DegreeErrorLR(1, i, N) = mean(DegreeError_cell{1, i});
        DegreeSTDLR(1, i, N) = std(DegreeError_cell{1, i});

        TranErrorLR(1, i, N) = mean(TranError_cell{1, i});
        TranSTDLR(1, i, N) = std(TranError_cell{1, i});

        DegreeErrorLR(2, i, N) = mean(DegreeError_cell{2, i});
        DegreeSTDLR(2, i, N) = std(DegreeError_cell{2, i});

        TranErrorLR(2, i, N) = mean(TranError_cell{2, i});
        TranSTDLR(2, i, N) = std(TranError_cell{2, i});
    end
end

figure(1), 
errorbar(3:12, squeeze(TranErrorLR(1, 1, :)), squeeze(TranSTDLR(1, 1, :)), 'r'), hold on, 
errorbar(3:12, squeeze(TranErrorLR(1, 2, :)), squeeze(TranSTDLR(1, 2, :)), 'g'), 
errorbar(3:12, squeeze(TranErrorLR(1, 3, :)), squeeze(TranSTDLR(1, 3, :)), 'b'), 
errorbar(3:12, squeeze(TranErrorLR(1, 4, :)), squeeze(TranSTDLR(1, 4, :)), 'c'), 
xlabel('Number of motion'), ylabel('error (mm)'),
title('translational error (left)')
legend('ATA', 'Barreto', 'TSAI', 'DQ'), grid on, 
hold off
figure(2), 
errorbar(3:12, squeeze(DegreeErrorLR(1, 1, :)), squeeze(DegreeSTDLR(1, 1, :)), 'r'), hold on,
errorbar(3:12, squeeze(DegreeErrorLR(1, 2, :)), squeeze(DegreeSTDLR(1, 2, :)), 'g'), 
errorbar(3:12, squeeze(DegreeErrorLR(1, 3, :)), squeeze(DegreeSTDLR(1, 3, :)), 'b'), 
errorbar(3:12, squeeze(DegreeErrorLR(1, 4, :)), squeeze(DegreeSTDLR(1, 4, :)), 'c'), 
xlabel('Number of motion'), ylabel('error (degree)')
title('rotational error (left)')
legend('ATA', 'Barreto', 'TSAI', 'DQ'), grid on, 
hold off

% figure(3),
% plot(3:12, DQ_produce_imag, 'LineWidth', 2), xlabel('Number of motion'), ylabel('Number of imaginary result produced from DQ'), grid on, axis tight

% save resultBNMotion.mat