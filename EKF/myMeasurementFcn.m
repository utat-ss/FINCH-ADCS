function y = myMeasurementFcn(x)
% x = [q0123;p;q;r]
q0123 =x(1:4);
pqr = x(5:7);
y = [q0123;pqr]; 
end