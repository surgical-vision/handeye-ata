function R = Rot(x, option)

    R = eye(3);
    
    if (strcmp(option, 'x'))
        R = [1, 0, 0;
             0, cos(x), -sin(x);
             0, sin(x),  cos(x)];
    elseif (strcmp(option, 'y'))
        R = [cos(x), 0, sin(x);
            0, 1, 0;
            -sin(x), 0, cos(x)];
    elseif (strcmp(option, 'z'))
        R = [cos(x), -sin(x), 0;
             sin(x),  cos(x), 0;
             0, 0, 1];
    else
        error('invalid input.');
    end

end