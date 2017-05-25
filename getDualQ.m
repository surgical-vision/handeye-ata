function [q,qprime] = getDualQ(R,t)    
    %Conversion from rigid transformation R,t to the dual quaternion
    r = rodrigues(R);
    theta = norm(r);
    l = r/norm(theta);
    %Get real part
    q = [sin(theta/2)*l; cos(theta/2)];
    %Get dual part
    qprime = [.5*(q(4)*t+cross(t,q(1:3)));-.5*q(1:3)'*t]; 
%     qprime = 0.5*qmult([0; t], q);
%     qprime3 = [qprime2(4); qprime2(1:3)];
%     disp(q);