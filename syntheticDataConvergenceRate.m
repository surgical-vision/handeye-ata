%%%%SYNTHETIC DATA EXPERIMENT
%%%Test the convergence rate of our method when different initialisation
%%%methods are applied. Note that this code does not involve the use of
%%%Levenberg Marquadt at the end of solution as we focus on the solution
%%%before the refinement.

clc;
% clear;

R_w = [0 0 0]'; %Robot origin
G_w = [-0.3, 0.5,  0.7]'; %Calibration grid position

% 0 for using Algorithm 1 (monocular)
% 1 for using Algorithm 2(stereo)
using_stereo = 1;

% 0 for using small displacement between camera and robot's gripper.
% 1 for using big displacement between camera and robot's gripper.
small_or_big = 0;

% Maximum iteration for our algorithm
max_ite = 15;

% Number of motion included in calibration
N_motion = 6;

%State how many times the program is simulated per one noise coefficient
nTimes = 100;

%Noise coefficient
NoiseCoeffT = 1e-4;
NoiseCoeffR = pi/180;

%Initialize error matrix
DegreeError = zeros(max_ite, 4);
DegreeSTD = zeros(max_ite, 4);
TranError = zeros(max_ite, 4);
TranSTD = zeros(max_ite, 4);

DegreeError_buf = zeros(max_ite, nTimes, 4);
TranError_buf = zeros(max_ite, nTimes, 4);

for j = 1:nTimes
    
    %Calculate robot base to calibration grid transformation
    gridTbase = [1 0 0 G_w(1);0 1 0 G_w(2);0 0 1 G_w(3);0 0 0 1]*[Rot(-pi/2, 'z'), zeros(3, 1);0 0 0 1];
    baseTgrid = inv(gridTbase);
    
    %stereo matrix
    Rx = rodrigues(0.5*pi/180*randn(3, 1));
    %The transformation between left and right camera
    rightTleft = [Rx, 0.005*randn(3,1); 0 0 0 1];
    
    if (small_or_big)
        %Create eye to hand transformation.
        Rx = rodrigues(pi/18*randn(3, 1));
        %The hand-eye transformation for left camera: XL
        camLTgripper = [Rx, 10*randn(3,1); 0 0 0 1];
    else
        %Create eye to hand transformation.
        Rx = rodrigues(pi/180*randn(3, 1));
        %The hand-eye transformation for left camera: XL
        camLTgripper = [Rx, 0.001*randn(3,1); 0 0 0 1];
    end
    
    %The hand-eye transformation for right camera: XR
    camRTgripper = rightTleft*camLTgripper;

    fprintf('trial number = %d\n\n', j);         

    %Initialize gripper transformation and camera transformation
    camLTgrid = zeros(4, 4, N_motion);
    camRTgrid = zeros(4, 4, N_motion);
    baseTgripper = zeros(4, 4, N_motion);
    gripperTbase = zeros(4, 4, N_motion);

    ranT = 0.4*randn(3, N_motion);
    ranR = 30*pi*randn(3, N_motion)/180;
    noiseT = NoiseCoeffT*randn(3, N_motion);
    noiseR = NoiseCoeffR*randn(3, N_motion);
    noiseTcamL = NoiseCoeffT*randn(3, N_motion);
    noiseRcamL = NoiseCoeffR*randn(3, N_motion);
    noiseTcamR = NoiseCoeffT*randn(3, N_motion);
    noiseRcamR = NoiseCoeffR*randn(3, N_motion);
    
    gripperTbaseNoise = zeros(size(gripperTbase)); %robot 2 world (from forward kin)
    baseTgripperNoise = zeros(size(baseTgripper)); %world 2 robot
    camLTgridNoise = zeros(size(camLTgrid));
    camRTgridNoise = zeros(size(camRTgrid));
    
    %Generate robot and camera motions (with and without noise)
    for i = 1:N_motion     
        %Robot and camera motions without noise
        gripperTbase(:, :, i) = [rodrigues(ranR(:, i)), ranT(:, i);0 0 0 1];
        baseTgripper(:, :, i) = inv(gripperTbase(:, :, i));
        camLTgrid(:, :, i) = camLTgripper*gripperTbase(:, :, i)*baseTgrid;
        camRTgrid(:, :, i) = camRTgripper*gripperTbase(:, :, i)*baseTgrid;
        
        %Robot and camera motions with noise
        gripperTbaseNoise(:, :, i) = gripperTbase(:, :, i)*[rodrigues(noiseR(:, i)), noiseT(:, i);0 0 0 1];
        baseTgripperNoise(:, :, i) = inv(gripperTbaseNoise(:, :, i));
        camLTgridNoise(:, :, i) = camLTgrid(:, :, i)*[rodrigues(noiseRcamL(:, i)), noiseTcamL(:, i);0 0 0 1];
        camRTgridNoise(:, :, i) = camRTgrid(:, :, i)*[rodrigues(noiseRcamR(:, i)), noiseTcamR(:, i);0 0 0 1];
    end

    TA(:, :, :, 1) = camRTgridNoise;
    TA(:, :, :, 2) = camLTgridNoise;
    TB = baseTgripperNoise;
    
    %Run hand-eye calibration with different initialisation methods
    for i = 1:4
        %Matrix C is a 6x(max_ite). Each column consists of 6-DOF rigid
        %transformation: 3 for rotation and 3 for translation. The first
        %three elements need to be converted using rodrigues formula.
        [~, ~, C] = HandEye_test_convergence(TA, TB, using_stereo, i, rightTleft);
        
        %Since Dual Quaternion method sometimes produces imaginary number
        %in the solution and we do not have proper solution for this yet,
        %we force them to be real.
        if (i == 3)
            C = real(C);
        end
        
        %Calculate error
        for k = 1:max_ite
            estimate_x = [rodrigues(C(1:3, k)), C(4:6, k);0 0 0 1];
            dT = camRTgripper/estimate_x;
            
            DegreeError_buf(k, j, i) = norm(rodrigues(dT(1:3, 1:3)))*180/pi;
            TranError_buf(k, j, i) = norm(camRTgripper(1:3, 4) - estimate_x(1:3, 4))*1000;
        end
    end



end

%Calculate mean and standard deviation
for i = 1:4
    for j = 1:max_ite
        DegreeError(j, i) = mean(deleteoutliers(DegreeError_buf(j, :, i)));
        TranError(j, i) = mean(deleteoutliers(TranError_buf(j, :, i)));
        DegreeSTD(j, i) = std(deleteoutliers(DegreeError_buf(j, :, i)));
        TranSTD(j, i) = std(deleteoutliers(TranError_buf(j, :, i)));
    end
end

%Plot the results
figure(1),
set(gca, 'FontSize', 16),
semilogy(1:max_ite, DegreeError(:, 1), 'r', 'LineWidth', 2), hold on,
semilogy(1:max_ite, DegreeError(:, 2), 'g', 'LineWidth', 2),
semilogy(1:max_ite, DegreeError(:, 3), 'b', 'LineWidth', 2),
semilogy(1:max_ite, DegreeError(:, 4), 'k', 'LineWidth', 2),
set(gca, 'FontSize', 16), 
xlabel('Iteration Number'), ylabel('Error in rotation vector (degree)'),
legend('IDQ', 'TSAI', 'DQ', 'Identity', 'Orientation', 'horizontal'), grid on, hold off

figure(2),
set(gca, 'FontSize', 16),
semilogy(1:max_ite, TranError(:, 1), 'r', 'LineWidth', 2), hold on,
semilogy(1:max_ite, TranError(:, 2), 'g', 'LineWidth', 2),
semilogy(1:max_ite, TranError(:, 3), 'b', 'LineWidth', 2),
semilogy(1:max_ite, TranError(:, 4), 'k', 'LineWidth', 2),
set(gca, 'FontSize', 16), 
xlabel('Iteration Number'), ylabel('Error in translation vector (mm)'),
legend('IDQ', 'TSAI', 'DQ', 'Identity'), grid on, hold off

figure(3),
plot(1:max_ite, DegreeSTD(:, 1), 'r', 'LineWidth', 2), hold on,
plot(1:max_ite, DegreeSTD(:, 2), 'g', 'LineWidth', 2),
plot(1:max_ite, DegreeSTD(:, 3), 'b', 'LineWidth', 2),
plot(1:max_ite, DegreeSTD(:, 4), 'k', 'LineWidth', 2),
xlabel('Iteration Number'), ylabel('standard deviation in rotation vector (degree)'),
legend('IDQ', 'TSAI', 'DQ', 'Identity'), grid on, hold off

figure(4),
plot(1:max_ite, TranSTD(:, 1), 'r', 'LineWidth', 2), hold on,
plot(1:max_ite, TranSTD(:, 2), 'g', 'LineWidth', 2),
plot(1:max_ite, TranSTD(:, 3), 'b', 'LineWidth', 2),
plot(1:max_ite, TranSTD(:, 4), 'k', 'LineWidth', 2),
xlabel('Iteration Number'), ylabel('standard deviation in translation vector (mm)'),
legend('IDQ', 'TSAI', 'DQ', 'Identity'), grid on, hold off

%Save results for plotting-gui
% if (small_or_big)
%     save resultConvergeHigh.mat
% else
%     save resultConvergeLow.mat
% end