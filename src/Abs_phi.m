% Absolute joint class
classdef Abs_phi < Joint
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
            Phi = q(3) - obj.constant(3);
        end
        function Jac = getJac(obj,q)
            Jac = [0,0,1];
        end
        function Nu = getNu(obj,q)
            Nu = 0;
        end
        function Gamma = getGamma(obj,q,q_dot)
            Gamma = 0;
        end
    end
end