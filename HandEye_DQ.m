function X = HandEye_DQ(A, B)
%Our quaternions are like this (q1 q2 q3,s )
    %Get n
    n = size(A, 3);
    %Make movements (a,B) which are interposition transformations
    %(marker2wordl and cam2grid)
    %transform A,B into dual quaternions        
    for i=1:n
        [q,qprime] = getDualQ(A(1:3,1:3,i),A(1:3,4,i));               
        Qa(i).q = q;
        Qa(i).qprime = qprime;
        [q,qprime] = getDualQ(B(1:3,1:3,i),B(1:3,4,i));                
        Qb(i).q = q;
        Qb(i).qprime = qprime;
    end
    
    %The dual quaternion is (Q.q + epsilon*Q.prime)
    %a = Qa.q, a' = Qa.prime  idem for b
    S = [];
    for i=1:n
        S(:,:,i) = [Qa(i).q(1:3)-Qb(i).q(1:3)   skew3(Qa(i).q(1:3)+Qb(i).q(1:3)) zeros(3,1) zeros(3,3);...
                    Qa(i).qprime(1:3)-Qb(i).qprime(1:3)   skew3(Qa(i).qprime(1:3)+Qb(i).qprime(1:3)) Qa(i).q(1:3)-Qb(i).q(1:3)   skew3(Qa(i).q(1:3)+Qb(i).q(1:3))];                                
    end  
    
    %Construct T
    T = [];    
    for i=1:n
        T = [T  S(:,:,i)'];      
    end
    
    T = T';
    %SVD 
    [U,S,V] = svd(T);
    
    %Solution, right null vectors of T
    v7 = V(:,7);
    v8 = V(:,8);
    
    u1 = v7(1:4);
    v1 = v7(5:8);
    
    u2 = v8(1:4);
    v2 = v8(5:8);
    %Now lambda1*v7+lambda2*v8 = [q;qprime]
    %
    %or other:
    %
    %lambda1^2*u1'*u1+2*lambda1*lambda2*u1'*u2+lambda2^2*u2'*u2 = 1   
    %and
    %lambda1^2*u1'*v1 + lambda1*lambda2*(u1'*v2+u2'*v1)+lambda2*u2'*v1 = 0
    %Setting lambda1/lambda2 = s
    %lambda1^2/lambda2^2*u1'*v1 + lambda1*lambda2/lambda2^2*(u1'*v2+u2'*v1)+lambda2^2/lambda2^2*u2'*v1 = 0
    %s^2*u1'*v1 + s*(u1'*v2+u2'*v1)+u2'*v1 = 0
    %s^2*u1'*v1 + s*(u1'*v2+u2'*v1)+u2'*v1 = 0
    a = u1'*v1;
    b = (u1'*v2+u2'*v1);
%     c = u2'*v1;
    c = u2'*v2 ;
    s = roots([a b c]);
    
    s = real(s);
    
    %insert into equation
    val1 = s(1)^2*u1'*u1+2*s(1)*u1'*u2+u2'*u2;
    val2 = s(2)^2*u1'*u1+2*s(2)*u1'*u2+u2'*u2;
    %Take bigger value
    if(val1>val2)
        s = s(1);
        val = val1;
    else
        s = s(2);
        val = val2;
    end
    %Get lambdas
    lambda2 = sqrt(1/val);
    lambda1 = s*lambda2;
    
    %This algorithm gives quaternion with the form of (s, q1 q2
    %q3)->contrary to the notation we used above (q1 q2 q3,s )
    %Therefore we must rearrange the elements!        
    qfinal = lambda1*v7+lambda2*v8;    
    q = [qfinal(2:4);qfinal(1)];
    qprime = [qfinal(6:8);qfinal(5)];
    
    %Extract transformation
    R = q2dcm(q)';    
    t = 2*qmult(qprime,qconj(q));
    t = t(1:3);
    
    %Assign output arguments
    X = [R t;[0 0 0 1]];
    err=[];
end