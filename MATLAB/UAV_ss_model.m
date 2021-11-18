% UAV Longutdinal State Space Model -- Generated From FlightSimulator Project
% Linearization is perfrormed for steady state flight with 
% velocity of 17 m/sec, flight angle of 0 rad, altitude of 0 m and orbit radius of 0 m.

%% trimmed states and control inputs : 
% states = [V, alpha, beta, p, q, r, phi, theta, psi, Xe, Ye, h]
x_trim = [17, 0.158617, 0, 0, 0, 0, 0, 0.158617, 0, 0, 0, 0]';
% controls = [throttle, elevator, aileron, rudder]
u_trim = [0.512061, -0.425364, 0.00509932, 0]';

%% Longitudinal Linear Dynamics
A_lon = [
-0.191467	9.40567	2.2799e-11	-9.81	0;
-0.0660877	-3.04424	0.976062	2.82685e-10	0;
4.56259e-11	-46.2157	-3.60042	0	0;
0	0	1	0	0;
-4.88498e-10	-17	0	17	0;
];
B_lon = [
6.32923	-0.123697;
-0.0595545	-0.070068;
0	-16.6984;
0	0;
0	0;
];
C_lon = eye(5);
D_lon = 0;

uav_lon = ss(A_lon, B_lon, C_lon, D_lon,...
'StateName', { 'V','alpha','q','theta','h' }, ...
'InputName', { 'throttle','elevator' }, ...
'OutputName', { 'V','alpha','q','theta','h' });
disp('eig(A_lon) = '); disp(eig(A_lon))

%% Lateral-Directional Linear Dynamics
A_lat = [
-0.528205	0.157953	-0.987447	0.569815	0;
-44.6996	-15.3876	7.41543	0	0;
9.05235	-0.0782623	-0.834805	0	0;
0	1	0.159961	0	0;
0	0	1.01271	0	0;
];
B_lat = [
0.0404239	0.102407;
60.5206	-0.830644;
2.31743	-11.5051;
0	0;
0	0;
];
C_lat = eye(5);
D_lat = 0;

uav_lat = ss(A_lat, B_lat, C_lat, D_lat,...
'StateName', { 'beta','p','r','phi','psi' }, ...
'InputName', { 'aileron','rudder' }, ...
'OutputName', { 'beta','p','r','phi','psi' });
disp('eig(A_lat) = '); disp(eig(A_lat))

