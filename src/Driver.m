% Driver class
classdef Driver
    properties
        body
        Jac
        Phi
        Phi_t
        Phi_qt
        Phi_tt
        Phi_qq_dotq
    end
    methods
        function System = Driver(body,Phi,Jac,Phi_t,Phi_qt,Phi_tt,Phi_qq_dotq)
            System.body = body;
            System.Jac = Jac;
            System.Phi = Phi;
            
            try
                System.Phi_t = Phi_t;
            catch
                System.Phi_t = @(q,t) 0;
            end
            
            try
                System.Phi_tt = Phi_tt;
            catch
                System.Phi_tt = @(q,t) 0;
            end
            
            try
                System.Phi_qt = Phi_qt;
            catch
                System.Phi_qt = @(q,t) [0,0,0];
            end
            
            try
                System.Phi_qq_dotq = Phi_qq_dotq;
            catch
                System.Phi_qq_dotq = @(q,t) [0,0,0];
            end
        end
        function Phi = getPhi(obj,q,t)
            Phi = obj.Phi(q,t);
        end
        function Jac = getJac(obj,q,t)
            Jac = obj.Jac(q,t);
        end
        function Nu = getNu(obj,q,t)
            Nu = - obj.Phi_t(q,t);
        end
        function Gamma = getGamma(obj,q,q_dot,t)
            Gamma = - obj.Phi_qq_dotq(q,t) * q_dot + 2 * obj.Phi_qt(q,t) ...
                * q_dot - obj.Phi_tt(q,t);
        end
    end
end