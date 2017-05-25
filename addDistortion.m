function [xd, xp] = addDistortion(x, cc, kc)    
    
    m = length(x);
    
    r2 = x(1, :).^2 + x(2, :).^2;
    a = x(1, :).*x(2, :);
    b = 2*x(1, :).^2;
    c = 2*x(2, :).^2;
    
    dx = [2*kc(3)*a + kc(4)*(r2 + b);
          kc(3)*(r2 + c) + 2*kc(4)*a];
    xd = repmat((1 + kc(1)*r2 + kc(2)*r2.^2 + kc(5)*r2.^3), 2, 1).*x(1:2, :) + dx + repmat([cc(1); cc(2)], 1, m);
    xp = x + repmat([cc(1); cc(2); 1], 1, m);
    
end