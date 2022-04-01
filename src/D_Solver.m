function [q,q_dot,q_ddot,t,lambda] = D_Solver(System,q0,q0_dot,t,alpha,beta)
% Solves dynamic systems from given initial guess values, maximum iteration
% and termination tolerance.
%
% by Mads Rosenh?j, Aarhus University, 2019.
%
% Input:
%   System  = System object, containing all joints and drivers
%   q0      = Initial position vector
%   q_dot0  = Initial velocities vector
%   t       = Time vector for results
%   kmax    = Maximum number of iterations
%   tol     = Termination tolerance
%
% Output:
%   q       = Position matrix
%   q_dot   = Velocity matrix
%   q_ddot  = Acceleration matrix
%   t       = time vector
%   lambda  = ...

try
    alpha; beta;
catch
    alpha = []; beta = [];
end
m = length(q0);

y0 = [q0;q0_dot];
opts = odeset('MaxStep',1e-1,'InitialStep',1e-5);
[t,y] = ode45(@(t,y) diff(t,y,System,alpha,beta), t, y0,opts);
t = t';     y = y';
q = y(1:m,:);
q_dot = y(m+1:m*2,:);

[m,~] = size(q_dot);
q_ddot = zeros(m,length(t));
lambda = zeros(length(System.Phi(q0,0)),length(t));



for i = 1:length(t)
    res = acc(t(i),y(:,i),System,alpha,beta,1);
    q_ddot(:,i) = res(1:m);
    lambda(:,i) = res(m+1:length(res));
end


    function dydt = diff(t,y,System,alpha,beta)
        n = length(y)/2;
        y_dot = y(n+1:n*2);
        y_ddot = acc(t,y,System,alpha,beta,0);
        dydt = [y_dot; y_ddot];
    end

    function y_ddot = acc(t,y,System,alpha,beta,type)
        [n,~] = size(System.M);
        y_dot = y(n+1:n*2);
        
        A = [System.M, System.Jac(y,t)';
            System.Jac(y,t), zeros(n-1,n-1)];
        
        try
            Gamma = System.Gamma_hat(y,y_dot,t,alpha,beta);
        catch
            Gamma = System.Gamma(y,y_dot,t);
        end
        B = [System.getQA(y,y_dot);
            Gamma];
        
        C = A\B;
        
        if type == 0
            y_ddot = C(1:n);
        else
            y_ddot = C;
        end
    end
end







