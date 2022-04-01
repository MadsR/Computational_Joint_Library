% Absolute joint class
classdef Abs_x < Joint
    properties
        body
        s_prime
        constant
    end
    methods
        function Joint = Abs_x(body,s_prime,constant)
            Joint.body = body;
            Joint.s_prime = s_prime;
            Joint.constant = constant;
            Joint.Type = 'Absolute joint';
            Joint.nb = 1;
        end
        function Phi = getPhi(obj,q)
            Phi = q(1) + obj.s_prime(1)*cos(q(3)) ...
                - obj.s_prime(2)*sin(q(3)) - obj.constant(1);
        end
        function Jac = getJac(obj,q)
            Jac = [1,0,-obj.s_prime(1)*sin(q(3)) - obj.s_prime(2)*cos(q(3))];
        end
        function Nu = getNu(obj,q)
            Nu = 0;
        end
        function Gamma = getGamma(obj,q,q_dot)
            Gamma = (obj.s_prime(1)*cos(q(3)) - obj.s_prime(2)*sin(q(3)))*q_dot(3)^2;
        end
    end
end