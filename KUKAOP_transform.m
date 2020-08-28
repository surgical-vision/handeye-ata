function T = KUKAOP_transform(data)
    data(:, 4:6) = data(:, 4:6)*pi/180;
    
    N = size(data, 1);
    
    T = zeros(4, 4, N);
    
    for i = 1:N
        
        angleZ = data(i, 4);
        angleY = data(i, 5);
        angleX = data(i, 6);
        
        Rot_Z = [cos(angleZ) -sin(angleZ) 0;
                sin(angleZ) cos(angleZ) 0;
                0 0 1];
        Rot_Y = [cos(angleY) 0 sin(angleY);
                0 1 0;
                -sin(angleY) 0 cos(angleY)];
        Rot_X = [1 0 0;
                0 cos(angleX) -sin(angleX);
                0 sin(angleX) cos(angleX)];
        tran = [eye(3), data(i, 1:3)';0 0 0 1];
        T(:, :, i) = tran*[Rot_Z, zeros(3, 1); 0 0 0 1]*[Rot_Y, zeros(3, 1); 0 0 0 1]*[Rot_X, zeros(3, 1); 0 0 0 1];
    
    end
end