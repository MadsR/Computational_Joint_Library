% Revolute joint class
classdef Rev < Joint
    properties
        body_i
        body_j
        s_prime_i
        s_prime_j
    end
    methods
        function Joint = Rev(body_i,body_j,s_prime_i,s_prime_j)
            Joint.body_i = body_i;
            Joint.body_j = body_j;
            Joint.s_prime_i = s_prime_i;
            Joint.s_prime_j = s_prime_j;
            Joint.Type = 'Revolute joint';
            Joint.nb = 2;
        end
        function Phi = getPhi(obj,q)
            q_i = q(1:3);
            q_j = q(4:6);
            
            Phi = q_i(1:2) + A_mat(q_i(3))*obj.s_prime_i - ...
                q_j(1:2) - A_mat(q_j(3))*obj.s_prime_j;
        end
        function Jac_i = getJac_i(obj,q)
            I = [1,0;0,1];
            q_i = q(1:3);
            Jac_i = [I,B_mat(q_i(3))*obj.s_prime_i];
        end
        function Jac_j = getJac_j(obj,q)
            I = [1,0;0,1];
            q_j = q(4:6);
            
            Jac_j = [-I,-B_mat(q_j(3))*obj.s_prime_j];
        end
        function Nu = getNu(obj,q)
            Nu = [0;0];
        end
        function Gamma = getGamma(obj,q,q_dot)
            q_i = q(1:3);
            q_i_dot = q_dot(1:3);
            q_j = q(4:6);
            q_j_dot = q_dot(4:6);
            Gamma = A_mat(q_i(3)) * obj.s_prime_i * q_i_dot(3)^2 - ...
                A_mat(q_j(3)) * obj.s_prime_j * q_j_dot(3)^2;
        end
    end
end