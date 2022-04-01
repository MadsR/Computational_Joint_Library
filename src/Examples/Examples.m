%% Dynamic System
% clc, clear, close all
Q_A.Gravity(1) = QA(1,[0;-9.81*10;0],'Gravity');

% drag = @(q,q_dot) [0;-1/2*0.5*1*1*q_dot(2)*abs(q_dot(2));0];
% Q_A.Drag(1) = QA(1,drag,'Drag');

J(1) = 10;
J(2) = 1 * 1^2 / 12;

Joints.Abs(1) = Abs(1,[-1;0],[0;0]);
% Joints.Abs(2) = Abs_x(1,[0;1],[0;0]);


Pendulum = System(Joints,[],Q_A,J);

q0 = [0; 0; 0];
q0_dot = [0; 0; 0];
t = [0 10];

[q,q_dot,q_ddot,t] = D_Solver(Pendulum,q0,q0_dot,t,0,1);

%% Plots
K.L = 1;
figure()
% Trajectory
subplot(3,3,[1 4 7])
plot(0,0, 'k+'), hold on
plot(q(1,:),q(2,:)), grid, title('Trajectory'), xlabel('x_1'), ylabel('y_1'), axis equal, 
% xlim([-K.L K.L]), ylim([-K.L-0.1 0.1]), hold off

% Position
subplot(3,3,2)
plot(t,q(1,:)), grid, hold on, title('Position'), xlabel('t [s]'), ylabel('Displacement [m]') 
plot(t,q(2,:))
legend('x_1','y_1')

% Velocity
subplot(3,3,5)
plot(t,q_dot(1,:)), grid, hold on, title('Velocity'), xlabel('t [s]'), ylabel('Velocity [m/s]') 
plot(t,q_dot(2,:))
legend('x_{1dot}','y_{1dot}')

% Acceleration
subplot(3,3,8)
plot(t,q_ddot(1,:)), grid, hold on, title('Acceleration'), xlabel('t [s]'), ylabel('Acceleration [m/s^2]') 
plot(t,q_ddot(2,:))
legend('x_{1ddot}','y_{1ddot}')

% Angular position
subplot(3,3,3)
plot(t,q(3,:)), grid, hold on, title('Angular position'), xlabel('t [s]'), ylabel('Angle [rad]') 
legend('\theta_{1}')

% Angular velocity
subplot(3,3,6)
plot(t,q_dot(3,:)), grid, hold on, title('Angular velocity'), xlabel('t [s]'), ylabel('Angular velocity [rad/s]') 
legend('\theta_{1dot}')

% Angular acceleration
subplot(3,3,9)
plot(t,q_ddot(3,:)), grid, hold on, title('Angular acceleration'), xlabel('t [s]'), ylabel('Angular acceleration [rad/s^2]') 
legend('\theta_{1ddot}')

%% Kinematic System
clc, clear
s_prime_i = [-1;0];
s_prime_j = [-1;2];
s_prime_k = [1;2];
 

% 
Joints.Abs(1) = Abs(1,[-1.5;0],[0;0]);
Joints.Abs(2) = Abs(2,[2;0],[2;0]);
Joints.Abs(3) = Abs(3,[3;0],[1;0]);

Joint.Trans(1) = Trans(1,2,s_prime_i,s_prime_j,[1,0],[1,1]);
 
Joints.Rev(1) = Rev(3,2,s_prime_i,s_prime_j);
Joints.Rev(2) = Rev(2,3,s_prime_j,s_prime_k);

Phi = @(q,t) q(3)-(3/2*pi+pi/4*cos(2*pi*t));
Jac = @(q,t) [0, 0, 1];
Phi_t = @(q,t) pi^2/2 * sin(2*pi*t);
Phi_qt = @(q,t) [0,0,0];
Phi_tt = @(q,t) pi^3*cos(2*pi*t);
Phi_qq_dotq = @(q,t) [0,0,0];

Pend_driver = Driver(1,Phi,Jac,Phi_t,Phi_qt,Phi_tt,Phi_qq_dotq);

Pendulum = System(Joints,Pend_driver);


q0 = [0; 0; -pi/4; 0; 3; 4; 5; 6; 7];

% Pendulum.Phi(q,0)
% Pendulum.Jac(q0,0)
% Pendulum.Nu(q,0.1)
% Pendulum.Gamma(q,q_dot,0)

kmax = 50; tol = 1E-8;
% q0 = [0.2; -0.2; -pi/4];
t = 0:0.01:4;

[q,q_dot,q_ddot,t] = K_Solver(Pendulum,q0,t,kmax,tol);


%% Plots
K.L = 1;
figure()
% Trajectory
subplot(3,3,[1 4 7])
plot(0,0, 'k+'), hold on
plot(q(1,:),q(2,:)), grid, title('Trajectory'), xlabel('x_1'), ylabel('y_1'), axis equal, 
% xlim([-K.L K.L]), ylim([-K.L-0.1 0.1]), hold off

% Position
subplot(3,3,2)
plot(t,q(1,:)), grid, hold on, title('Position'), xlabel('t [s]'), ylabel('Displacement [m]') 
plot(t,q(2,:))
legend('x_1','y_1')

% Velocity
subplot(3,3,5)
plot(t,q_dot(1,:)), grid, hold on, title('Velocity'), xlabel('t [s]'), ylabel('Velocity [m/s]') 
plot(t,q_dot(2,:))
legend('x_{1dot}','y_{1dot}')

% Acceleration
subplot(3,3,8)
plot(t,q_ddot(1,:)), grid, hold on, title('Acceleration'), xlabel('t [s]'), ylabel('Acceleration [m/s^2]') 
plot(t,q_ddot(2,:))
legend('x_{1ddot}','y_{1ddot}')

% Angular position
subplot(3,3,3)
plot(t,q(3,:)), grid, hold on, title('Angular position'), xlabel('t [s]'), ylabel('Angle [rad]') 
legend('\theta_{1}')

% Angular velocity
subplot(3,3,6)
plot(t,q_dot(3,:)), grid, hold on, title('Angular velocity'), xlabel('t [s]'), ylabel('Angular velocity [rad/s]') 
legend('\theta_{1dot}')

% Angular acceleration
subplot(3,3,9)
plot(t,q_ddot(3,:)), grid, hold on, title('Angular acceleration'), xlabel('t [s]'), ylabel('Angular acceleration [rad/s^2]') 
legend('\theta_{1ddot}')