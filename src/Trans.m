% Translational joint class
classdef Trans < Joint
    properties
        body_i
        body_j
        s_prime_i
        s_prime_j
        v_prime_i
        v_prime_j
    end
    methods
        function Joint = Trans(body_i,body_j,s_prime_i,s_prime_j,v_prime_i,v_prime_j)
            Joint.body_i = body_i;
            Joint.body_j = body_j;
            Joint.s_prime_i = s_prime_i;
            Joint.s_prime_j = s_prime_j;
            Joint.v_prime_i = v_prime_i;
            Joint.v_prime_j = v_prime_j;
            Joint.Type = 'Translational joint';
            Joint.nb = 2;
        end
        function Phi = getPhi(obj,q)
            q_i = q(1:3);
            q_j = q(4:6);
            rp_i = q_i(1:2) + A_mat(q_i(3)) * obj.s_prime_i;
            rp_j = q_j(1:2) + A_mat(q_j(3)) * obj.s_prime_j;
            d_ij = rp_j - rp_i;
            B_ij = B_mat(q_j(3) - q_i(3));
            
            R = [0,-1;1,0];
            
            Phi = [obj.v_prime_i' * B_mat(q_i(3))' * d_ij - obj.v_prime_i' * ...
                B_ij * obj.s_prime_j - obj.v_prime_i' * R' * obj.v_prime_i;
                -obj.v_prime_i' * B_mat(q_j(3) - q_i(3)) * obj.v_prime_j];
        end
        
        function Jac_i = getJac_i(obj,q)
            q_i = q(1:3);
            q_j = q(4:6);
            rp_i = q_i(1:2) + A_mat(q_i(3)) * obj.s_prime_i;
            rp_j = q_j(1:2) + A_mat(q_j(3)) * obj.s_prime_j;
            d_ij = rp_j - rp_i;
            A_ij = A_mat(q_j(3) - q_i(3));
            
            Jac_i = [-obj.v_prime_i' * B_mat(q_i(3))', -obj.v_prime_i' * A_mat(q_i(3))' * ...
                d_ij - obj.v_prime_i' * A_ij * obj.s_prime_j;
                0, 0, -obj.v_prime_i' * A_ij * obj.v_prime_j];
        end
        
        function Jac_j = getJac_j(obj,q)
            q_i = q(1:3);
            q_j = q(4:6);
            A_ij = A_mat(q_j(3) - q_i(3));
            
            
            Jac_j = [obj.v_prime_i' * B_mat(q_i(3))' , obj.v_prime_i' * ...
                A_ij * obj.s_prime_j;
                0, 0, obj.v_prime_i' * A_ij * obj.v_prime_j];
        end
        function Nu = getNu(obj,q)
            Nu = [0;0];
        end
        function Gamma = getGamma(obj,q,q_dot)
            R = [0,-1;1,0];
            q_i = q(1:3);
            q_j = q(4:6);
            qd_i = q_dot(1:3);
            qd_j = q_dot(4:6);
            rp_i = q_i(1:2) + A_mat(q_i(3)) * obj.s_prime_i;
            rp_j = q_j(1:2) + A_mat(q_j(3)) * obj.s_prime_j;
            d_ij = rp_j - rp_i;
            rp_i_d = qd_i(1:2) + qd_i(3) * A_mat(q_i(3))*R*obj.s_prime_i;
            rp_j_d = qd_j(1:2) + qd_j(3) * A_mat(q_j(3))*R*obj.s_prime_j;
            d_ij_d = rp_j_d - rp_i_d;
            B_ij = B_mat(q_j(3) - q_i(3));
            
            Gamma = -[obj.v_prime_i' * (B_ij * obj.s_prime_j * (qd_j(3) - qd_i(3))^2 ...
                - B_mat(q_i(3))' * d_ij * qd_i(3)^2 - 2 * A_mat(q_i(3))' *...
                d_ij_d * qd_i(3));0];
        end
    end
end