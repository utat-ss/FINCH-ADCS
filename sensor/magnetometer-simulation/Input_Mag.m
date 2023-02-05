N = 100;
Ts = 0.1;
m.time = Ts*[0:N]';
c1 = int8(ones([101,1]));
c2 = int8(ones([101,1]));
c3 = int8(ones([101,1]));
m.signals(1).values = [c1 c2 c3];
m.signals(1).dimensions = 3;