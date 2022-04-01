% Applied force
classdef QA
    properties
        body            % Body force is applied to
        Q_A             % Applied force vector
        Desc            % Applied force description
        nb              % Number of bodies used in force
    end
    methods
        function System = QA(body,Q_A,description)
            System.body = body;
            System.Q_A = Q_A;
            System.nb = 1;
            
            try
                System.Desc = description;
            catch
            end
        end
        function Q_A = getQ_A(obj,q,q_dot)
            try
                q;q_dot;
                Q_A = obj.Q_A(q,q_dot);
            catch
                try
                    q;
                    Q_A = obj.Q_A(q);
                catch
                    Q_A = obj.Q_A;
                end
            end
        end
    end
end