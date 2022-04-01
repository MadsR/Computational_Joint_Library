% Absolute joint class
classdef Abs < Joint
    properties
        body
        s_prime
        constant
    end
    methods
        function Joint = Abs(body,s_prime,constant)
            Joint.body = body;
            Joint.s_prime = s_prime;
            Joint.constant = constant;
            Joint.Type = 'Absolute joint';
            Joint.nb = 1;
        end
        function Phi = getPhi(obj,q)
            Phi = q(1:2) + A_mat(q(3))*obj.s_prime - obj.constant;
        end
        function Jac = getJac(obj,q)
            I = [1,0;0,1];
            Jac = [I,B_mat(q(3))*obj.s_prime];
        end
        function Nu = getNu(obj,q)
            Nu = [0;0];
        end
        function Gamma = getGamma(obj,q,q_dot)
            Gamma = A_mat(q(3)) * obj.s_prime * q_dot(3)^2;
        end
    end
end