% Absolute distance joint class
classdef Abs_d < Joint
    properties
        body
        s_prime
        constant
        constant_dist
    end
    methods
        function Joint = Abs_d(body,s_prime,constant,constant_dist)
            Joint.body = body;
            Joint.s_prime = s_prime;
            Joint.constant = constant;
            Joint.constant_dist = constant_dist;
            Joint.Type = 'Absolute joint';
            Joint.nb = 1;
        end
        function Phi = getPhi(obj,q)
            r_p = q(1:2) + A_mat(q(3))*obj.s_prime;
            Phi = (r_p - obj.constant)' * (r_p - obj.constant) - obj.constant_dist^2;
        end
        function Jac = getJac(obj,q)
            r_p = q(1:2) + A_mat(q(3))*obj.s_prime;
            Jac = 2 * [(r_p - obj.constant)', obj.s_prime' * B_mat(q(3))' * (r_p - obj.constant)];
        end
        function Nu = getNu(obj,q)
            Nu = 0;
        end
        function Gamma = getGamma(obj,q,q_dot)
            r_p = q(1:2) + A_mat(q(3))*obj.s_prime;
            r_p_dot = q_dot(1:2) + q_dot(3) * B_mat(q(3)) * obj.s_prime;
            Gamma = -2 * (r_p_dot' * r_p_dot - obj.s_prime' * A_mat(q(3))' * ...
                (r_p - obj.constant) * q_dot(3)^2);
        end
    end
end