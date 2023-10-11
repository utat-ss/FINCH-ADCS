function f = eccentric_to_true(E,e)
    ee = sqrt((1+e)/(1-e));
    f = 2*atan( ee*tan(E/2) );
end