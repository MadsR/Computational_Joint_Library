function [q,q_dot,q_ddot,t] = K_Solver(System,q0,t,kmax,tol)
% Calculates roots to given function (func) from the jacobian function, and
% from given initial guess values, maximum iteration and termination
% tolerance.
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

m = length(q0);

q = zeros(m,length(t));
q_dot = zeros(m,length(t));
q_ddot = zeros(m,length(t));



for i = 1:length(t)
    q(:,i) = NewRaph(@System.Phi,@System.Jac,q0,t(i),kmax,tol);
    q0 = q(:,i);
    
    q_dot(:,i) = System.Jac(q0,t(i))\System.Nu(q0,t(i));
    q0_dot = q_dot(:,i);
    
    q_ddot(:,i) = System.Jac(q0,t(i))\System.Gamma(q0,q0_dot,t(i));
end
end



