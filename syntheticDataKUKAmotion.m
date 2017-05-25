%%%%SYNTHETIC DATA EXPERIMENT
%%%Check if the error in the robot positioning system is larger than the
%%%error in the camera calibraton.

clc;
clear;
close all;

%Generate grid
points.gt = zeros(4, 196);
points.gt(4, :) = 1;
points.gt(1, :) = repmat(0:3:39, [1, 14]);
for i = 1:14
    points.gt(2, 14*i-13:14*i) = 3*(i - 1);
end

KUKAdata = [589.24, 311.12, 439.30, 139.17, 72.24, 132.82;
            608.15, 304.09, 433.27, 139.14, 72.23, 108.08;
            600.10, 321.45, 438.77, 139.21, 72.24, 157.23;
            628.47, 308.94, 427.41, 140.86, 68.86, 167.94;
            475.69, 307.76, 436.63, 79.36, 75.34, 107.65;
            403.64, 304.35, 391.27, 47.75, 70.44, 49.62;
            399.65, 419.44, 441.20, -9.62, 79.49, -5.96;
            358.56, 427.92, 438.26, -18.08, 73.00, -39.72;
            390.97, 364.67, 434.63, 25.67, 77.50, 27.63;
            426.69, 386.39, 457.47, 38.27, 83.50, 64.66;
            673.06, 287.33, 427.41, 143.70, 59.50, 168.38;
            663.18, 309.38, 403.39, 151.66, 64.85, 128.30;
            663.19, 309.41, 403.51, 149.34, 64.20, 144.89;
            700.96, 314.31, 388.96, 154.51, 60.47, 149.46;
            739.61, 320.30, 369.87, 159.96, 55.03, 154.10;
            498.46, 270.77, 408.07, 95.37, 72.11, 91.70;
            513.73, 260.95, 406.47, 95.32, 72.12, 65.09;
            505.08, 282.56, 413.70, 95.31, 72.12, 111.59;
            491.83, 210.73, 387.86, 91.91, 63.58, 97.11;
            521.78, 583.97, 401.50, -101.85, 66.79, -98.01;
            527.80, 632.24, 386.56, -99.67, 58.73, -112.63;
            523.74, 644.50, 371.93, -99.68, 58.72, -82.68;
            524.79, 640.21, 392.80, -99.75, 58.74, -98.26;
            384.26, 569.31, 388.12, -58.12, 64.70, -58.08;
            392.06, 333.74, 421.94, 29.96, 74.22, 25.86;
            378.42, 349.31, 426.58, 29.92, 74.22, 48.96;
            384.80, 350.04, 428.25, 26.45, 76.26, 19.13;
            507.92, 564.25, 419.06, -95.84, 70.21, -96.90;
            509.73, 565.19, 419.13, -90.63, 67.99, -123.07;
            523.33, 564.22, 418.29, -105.73, 71.50, -78.23;
            544.29, 630.01, 395.19, -101.07, 60.39, -105.62;
            642.73, 198.52, 339.22, 124.65, 52.78, 120.63;
            645.31, 184.47, 329.30, 124.66, 52.77, 96.10;
            644.41, 184.47, 329.70, 122.74, 50.63, 129.07;
            650.60, 175.88, 355.50, 123.19, 50.63, 141.93;
            648.66, 173.15, 341.76, 122.29, 54.90, 71.84];
        
% Number of motion included in calibration
N_motion = size(KUKAdata, 1);
        
%Generate camera data
K.fcL = [846.7686, 920.1253]';
K.fcR = [850.0666, 920.7560]';
%Intrinsic matrix for left and right camera
K.K0L = [846.7686, 0, 0; 0, 920.1253, 0;0, 0, 1];
K.K0R = [850.0666, 0, 0; 0, 920.7560, 0;0, 0, 1];
%Distortion parameters
K.kcbpL = [-0.3616, 0.5989, 0.0038,  -0.015, 0];
K.kcbpR = [-0.3743, 0.5831, 0.0008, -0.0012, 0];
%Camera offset
K.ccL = [433.0905, 191.2167]';
K.ccR = [339.7955, 177.9082]';

R_w = [0 0 0]'; %Robot origin
G_w = [80, 50,  0]'; %Calibration grid position

% 0 for using Algorithm 1 (monocular)
% 1 for using Algorithm 2 (stereo)
using_stereo = 1;

% 0 for adding noise to only robot motion
% 1 for adding noise to both robot and camera motions
both_robot_and_camera_motion = 1;

%Noise coefficient
NoiseCoeffT = 0.025;
NoiseCoeffR = 0.1*pi/180;

%Calculate robot base to calibration grid transformation
Tbase2grid = [1 0 0 -G_w(1);0 1 0 G_w(2);0 0 1 G_w(3);0 0 0 1]*[Rot(-pi/2, 'z'), zeros(3, 1);0 0 0 1];
Tgrid2base = inv(Tbase2grid);

%Setting offset respected to end effector
%(Transformation from robot to camera)
X_off = -10;
Y_off = 10;
Z_off = 150;

%Create eye to hand transformation.
theta = 0.1*pi/180*randn(1);
alpha = 0.1*pi/180*randn(1);
phi = 0.1*pi/180*randn(1);
Rx = Rot(theta, 'y')*Rot(phi, 'x')*Rot(alpha, 'z');
%The hand-eye transformation for left camera: XL
Tgripper2camL = [eye(3), [X_off; Y_off; Z_off].*randn(3, 1); 0 0 0 1]*[Rx, zeros(3, 1);0 0 0 1];

theta = 0.1*pi/180*randn(1);
alpha = 0.1*pi/180*randn(1);
phi = 0.1*pi/180*randn(1);
Rx = Rot(theta, 'y')*Rot(phi, 'x')*Rot(alpha, 'z');
%The transformation between left and right camera
Tleft2right = [eye(3), [X_off; Y_off; Z_off].*randn(3, 1); 0 0 0 1]*[Rx, zeros(3, 1);0 0 0 1];
%The hand-eye transformation for left camera: XR
Tgripper2camR = Tleft2right*Tgripper2camL;       

noiseT = NoiseCoeffT*randn(3, N_motion);
noiseR = NoiseCoeffR*randn(3, N_motion);
noiseTcamL = NoiseCoeffT*randn(3, N_motion);
noiseRcamL = NoiseCoeffR*randn(3, N_motion);
noiseTcamR = NoiseCoeffT*randn(3, N_motion);
noiseRcamR = NoiseCoeffR*randn(3, N_motion);

%Generate robot and camera motions (with and without noise)
for i = 1:N_motion     
    %Robot and camera motions without noise
    Tgripper2base = KUKAOP_transform(KUKAdata(i, :));
    Tbase2gripper = inv(Tgripper2base);
    Tgrid2camL = Tgripper2camL*Tbase2gripper*Tgrid2base;

    %Robot and camera motions with noise
    Tgripper2baseNoise = Tgripper2base*[rodrigues(noiseR(:, i)), noiseT(:, i);0 0 0 1];
    Tbase2gripperNoise = inv(Tgripper2baseNoise);
    Tgrid2camLNoise = Tgripper2camL*Tbase2gripperNoise*Tgrid2base;

    TcamL2grid = inv(Tgrid2camL);
    TcamL2gridNoise = inv(Tgrid2camLNoise);
    
    [xd_L(:, :, i), ~, ~, ~, ~, ~, ~] = project_points2(points.gt(1:3, :), rodrigues(TcamL2grid(1:3, 1:3)), TcamL2grid(1:3, 4), K.fcL, K.ccL, K.kcbpL, 0);
    [xd_L_noise(:, :, i), ~, ~, ~, ~, ~, ~] = project_points2(points.gt(1:3, :), rodrigues(TcamL2gridNoise(1:3, 1:3)), TcamL2gridNoise(1:3, 4), K.fcL, K.ccL, K.kcbpL, 0);
end

err = abs(xd_L - xd_L_noise);

err_avg = zeros(N_motion, 1);

for i = 1:N_motion
    err_avg(i) = mean(deleteoutliers(sqrt(err(1, :, i).^2 + err(2, :, i).^2)));
end
err_avg = deleteoutliers(err_avg);

[~, ind_min] = sort(err_avg);

figure, set(gca, 'FontSize', 14), hold on,
plot(1:length(err_avg), err_avg, 'LineWidth', 2), xlabel('Motion number'), ylabel('Reprojection error (pixel)'), xlim([1, 33]), grid on
hold off

for i = 1:3
    figure,
    set(gca, 'FontSize', 14), hold on
    plot(abs(xd_L(1, :, ind_min(i + 9)))/10, abs(xd_L(2, :, ind_min(i + 9)))/10, 'ro', 'LineWidth', 2), grid on, xlabel('X'), ylabel('Y'), 
    plot(abs(xd_L_noise(1, :, ind_min(i + 9)))/10, abs(xd_L_noise(2, :, ind_min(i + 9)))/10, 'bo', 'LineWidth', 2), hold off
end

