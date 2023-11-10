% this is the unit test for the orbital paramter calculation.

M2E()
E2T()


function [] = M2E()
    M = 27;
    e = 0.5;
    

    expect_E = 48.43417991;
    E_anom = mean_to_eccentric_anomaly(M,e,0.000000001);
    
    assert(abs(expect_E-E_anom)<0.00000001)

end

function [] = E2T()
    E = 48.43417991;
    e = 0.5;
    

    expect_M = 75.83971719;
    M_anom = eccentric_to_true_anomaly(E,e);
    disp(expect_M)
    disp(M_anom)
    assert(abs(expect_M-M_anom)<0.00000001)

end