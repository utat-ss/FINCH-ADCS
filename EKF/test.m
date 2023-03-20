% Initialize the AEKF parameters
dt = 0.01;  % time step
q = [1; 0; 0; 0];  % initial attitude quaternion
P = eye(6);  % initial state covariance matrix
Q = diag([0.01, 0.01, 0.01, 0.01, 0.01, 0.01]);  % process noise covariance
R = diag([0.1, 0.1, 0.1]);  % measurement noise covariance

% Simulation loop
for i = 1:N
    % Generate random measurements (accelerometer, gyroscope, magnetometer)
    [acc, gyro, mag] = get_measurements();
    
    % AEKF prediction step
    w = gyro - bias;  % angular rate measurement with bias compensation
    omega = [0, -w(1), -w(2), -w(3); w(1), 0, w(3), -w(2); w(2), -w(3), 0, w(1); w(3), w(2), -w(1), 0];  % skew-symmetric matrix
    F = eye(6) + 0.5*dt*[omega, zeros(4,2); zeros(2,4), zeros(2,2)];  % state transition matrix
    q = q + 0.5*quaternion_multiply(q, [0; w])*dt;  % quaternion update
    Qd = Q*dt;  % discrete process noise covariance
    P = F*P*F' + Qd;  % state covariance prediction
    
    % AEKF update step
    H = get_measurement_jacobian(q);  % measurement Jacobian matrix
    z = [acc; mag];  % measurement vector
    y = z - get_measurement_model(q);  % measurement residual
    S = H*P*H' + R;  % innovation covariance
    K = P*H'/S;  % Kalman gain
    q = quaternion_multiply(q, quaternion_exp(0.5*K*y));  % quaternion update
    P = (eye(6) - K*H)*P;  % state covariance update
end