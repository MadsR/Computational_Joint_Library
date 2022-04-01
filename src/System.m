% System class
classdef System
    properties
        Joints          % Structure containing all body joints
        Joint_names     % Names of Joints in Joints structure
        
        Drivers         % Structure of all drivers
        Drivers_names   % Names of drivers in driver class
        
        Q_A             % Structure containing classes with applied forces
        Q_A_names       % Names of Applied forces in structure
        
        M               % Inertia vector for bodies
        
        % Inputs:
        %       Joints  =   Structure containing all joint objectives
        %       driver  =   Optional structure containing all system drivers
        %       QA      =   Optional structure containing all system
        %                   applied forces
        %       M_J_vec =   Optional vector containing mass and intertia for 
        %                   system bodies in format:
        %                       M_J_vec(1) = mass_body_1
        %                       M_J_vec(2) = inertia_body_1
        %                       M_J_vec(3) = mass_body_2
        %                       M_J_vec(4) = inertia_body_2
    end
    methods
        function Setup = System(Joints,driver,QA,M_J_vec)
            Setup.Joints = Joints;
            Setup.Joint_names = string(fieldnames(Joints));
            
            try
                Setup.Drivers = driver;
                Setup.Drivers_names = string(fieldnames(driver));
            catch
            end
            
            try
                Setup.Q_A = QA;
                Setup.Q_A_names = string(fieldnames(QA));
            catch
            end
            
            try
                Setup.M = eye(length(M_J_vec)*3/2);
                index = 1;
                for x = 1:length(M_J_vec)/2
                    for i = 1:2
                        Setup.M(index,index) = M_J_vec(x*2 - 1);
                        index = index + 1;
                    end
                    Setup.M(index,index) = M_J_vec(x*2);
                    index = index + 1;
                end
            catch
            end
        end
        
        function getPhi = Phi(obj,q,t)
            
            i = 1;          % Joint number
            for n = 1:length(obj.Joint_names)       % n = Joint name number
                name = obj.Joint_names(n);
                Joint_Type = getfield(obj.Joints,name);
                
                if Joint_Type(1).nb == 1
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index = (Joint_i.body-1) * 3 + 1;           % Body index
                        q_i = q(index:index+2);                         % Body variables
                        
                        Phi_i = Joint_i.getPhi(q_i);
                        
                        getPhi(i:i+length(Phi_i)-1,1) = Phi_i;
                        i = i + length(Phi_i);
                    end
                else
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index_i = (Joint_i.body_i-1) * 3 + 1;           % Body index
                        q_i = q(index_i:index_i+2);                         % Body variables
                        
                        index_j = (Joint_i.body_j-1) * 3 + 1;           % Body index
                        q_j = q(index_j:index_j+2);                         % Body variables
                        
                        q_joint = [q_i;q_j];
                        
                        Phi_i = Joint_i.getPhi(q_joint);
                        
                        getPhi(i:i+length(Phi_i)-1,1) = Phi_i;
                        i = i + length(Phi_i);
                    end
                end
            end
            try
                for dr = 1:length(obj.Drivers)
                    index = (obj.Drivers.body-1) * 3 + 1;           % Body index
                    q_i = q(index:index+2);
                    Phi_i = obj.Drivers.getPhi(q_i,t);
                    getPhi(i:i+length(Phi_i)-1,1) = Phi_i;
                end
            catch
            end
        end
        function getJac = Jac(obj,q,t)
            i = 1;          % Joint number
            for n = 1:length(obj.Joint_names)       % n = Joint name number
                name = obj.Joint_names(n);
                Joint_Type = getfield(obj.Joints,name);
                
                if Joint_Type(1).nb == 1
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index = (Joint_i.body-1) * 3 + 1;           % Body index
                        q_i = q(index:index+2);                         % Body variables
                        
                        Jac_i = Joint_i.getJac(q_i);
                        [m,~] = size(Jac_i);
                        
                        getJac(i:i+m-1,index:index+2) = Jac_i;
                        i = i + m;
                    end
                else
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index_i = (Joint_i.body_i-1) * 3 + 1;           % Body index
                        q_i = q(index_i:index_i+2);                         % Body variables
                        
                        index_j = (Joint_i.body_j-1) * 3 + 1;           % Body index
                        q_j = q(index_j:index_j+2);                         % Body variables
                        
                        q_joint = [q_i;q_j];
                        
                        Jac_i = Joint_i.getJac_i(q_joint);
                        Jac_j = Joint_i.getJac_j(q_joint);
                        [m,~] = size(Jac_i);
                        
                        getJac(i:i+m-1,index_i:index_i+2) = Jac_i;
                        getJac(i:i+m-1,index_j:index_j+2) = Jac_j;
                        i = i + m;
                    end
                end
            end
            for dr = 1:length(obj.Drivers)
                index = (obj.Drivers.body-1) * 3 + 1;           % Body index
                q_i = q(index:index+2);
                Jac_i = obj.Drivers.getJac(q_i,t);
                [m,~] = size(Jac_i);
                getJac(i:i+m-1,index:index+2) = Jac_i;
            end
        end
        function getNu = Nu(obj,q,t)
            i = 1;          % Joint number
            for n = 1:length(obj.Joint_names)       % n = Joint name number
                name = obj.Joint_names(n);
                Joint_Type = getfield(obj.Joints,name);
                
                if Joint_Type(1).nb == 1
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index = (Joint_i.body-1) * 3 + 1;           % Body index
                        q_i = q(index:index+2);                         % Body variables
                        
                        Nu_i = Joint_i.getNu(q_i);
                        
                        getNu(i:i+length(Nu_i)-1,1) = Nu_i;
                        i = i + length(Nu_i);
                    end
                else
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index_i = (Joint_i.body_i-1) * 3 + 1;           % Body index
                        q_i = q(index_i:index_i+2);                         % Body variables
                        
                        index_j = (Joint_i.body_j-1) * 3 + 1;           % Body index
                        q_j = q(index_j:index_j+2);                         % Body variables
                        
                        q_joint = [q_i;q_j];
                        
                        Nu_i = Joint_i.getNu(q_joint);
                        
                        getNu(i:i+length(Nu_i)-1,1) = Nu_i;
                        i = i + length(Nu_i);
                    end
                end
            end
            try
                for dr = 1:length(obj.Drivers)
                    index = (obj.Drivers.body-1) * 3 + 1;           % Body index
                    q_i = q(index:index+2);
                    Nu_i = obj.Drivers.getNu(q_i,t);
                    getNu(i:i+length(Nu_i)-1,1) = Nu_i;
                end
            catch
            end
        end
        function getGamma = Gamma(obj,q,q_dot,t)
            i = 1;          % Joint number
            for n = 1:length(obj.Joint_names)       % n = Joint name number
                name = obj.Joint_names(n);
                Joint_Type = getfield(obj.Joints,name);
                
                if Joint_Type(1).nb == 1
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;    % Current joint
                        
                        index = (Joint_i.body-1) * 3 + 1;           % Body index
                        q_i = q(index:index+2);                     % Body variables
                        q_i_dot = q_dot(index:index+2);
                        
                        Gamma_i = Joint_i.getGamma(q_i,q_i_dot);
                        
                        getGamma(i:i+length(Gamma_i)-1,1) = Gamma_i;
                        i = i + length(Gamma_i);
                    end
                else
                    for j = 1:length(getfield(obj.Joints,name)) % Current Joint
                        Joint_i = Joint_Type(j) ;               % Current joint
                        
                        index_i = (Joint_i.body_i-1) * 3 + 1;           % Body index
                        q_i = q(index_i:index_i+2);                         % Body variables
                        q_i_dot = q_dot(index_i:index_i+2);
                        
                        index_j = (Joint_i.body_j-1) * 3 + 1;           % Body index
                        q_j = q(index_j:index_j+2);                         % Body variables
                        q_j_dot = q_dot(index_j:index_j+2);
                        
                        q_joint = [q_i;q_j];
                        q_joint_dot = [q_i_dot;q_j_dot];
                        
                        Gamma_i = Joint_i.getGamma(q_joint,q_joint_dot);
                        
                        getGamma(i:i+length(Gamma_i)-1,1) = Gamma_i;
                        i = i + length(Gamma_i);
                    end
                end
            end
            for dr = 1:length(obj.Drivers)
                index = (obj.Drivers.body - 1) * 3 + 1;           % Body index
                q_i = q(index:index+2);                      % Body variables
                q_i_dot = q_dot(index:index+2);
                Gamma_i = obj.Drivers.getGamma(q_i,q_i_dot,t);
                getGamma(i:i+length(Gamma_i)-1,1) = Gamma_i;
            end
        end
        
        function getGamma_hat = Gamma_hat(obj,q,q_dot,t,alpha,beta)
            Gamma = obj.Gamma(q,q_dot,t);
            Phi = obj.Phi(q,t);
            Phi_dot = obj.Phi_dot(q,q_dot,t);
            getGamma_hat = Gamma - 2 * alpha * Phi_dot - beta^2 * Phi;
        end
        function getPhi_dot = Phi_dot(obj,q,q_dot,t)
            getPhi_dot = obj.Jac(q,t) * q_dot - obj.Nu(q,t);
        end
        function getPhi_ddot = Phi_ddot(obj,q,q_dot,q_ddot,t)
            getPhi_ddot = obj.Jac(q,t) * q_ddot - obj.Gamma(q,q_dot,t);
        end
        
        function Q_A = getQA(obj,q,q_dot)
            [m,~] = size(obj.M);
            Q_A = zeros(m,1);
            for n = 1:length(obj.Q_A_names)       % n = Applied force name number
                name = obj.Q_A_names(n);
                Applied_Force = getfield(obj.Q_A,name);
                
                if Applied_Force(1).nb == 1
                    for j = 1:length(getfield(obj.Q_A,name)) % Current Force
                        Force_j = Applied_Force(j) ;         % Current Force
                        
                        index = (Force_j.body-1) * 3 + 1;           % Body index
                        q_i = q(index:index+2);                     % Body variables
                        q_i_dot = q_dot(index:index+2);
                        
                        Q_A_j = Force_j.getQ_A(q_i,q_i_dot);
                        
                        Q_A(index:index+length(Q_A_j)-1,1) = Q_A(index:index+length(Q_A_j)-1,1) + Q_A_j;
                    end
                else
                    % Used for several body forces
                end
            end
        end
        
        function res = Reaction_force(obj,q,lambda)
            
            
        end
    end
end


























