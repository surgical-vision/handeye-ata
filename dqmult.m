function dq = dqmult(p, q)

    real_p = p(1:4);
    dual_p = p(5:8);
    
    real_q = q(1:4);
    dual_q = q(5:8);
    
    dq(1:4) = qmult(real_p, real_q);
    dq(5:8) = qmult(real_p, dual_q) + qmult(dual_p, real_q);

end