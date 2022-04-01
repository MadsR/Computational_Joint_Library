% This function calculates the A matrix
% Inputs: 
% b = Angle in radians
function a = A_mat(b)
a = [cos(b),-sin(b);sin(b),cos(b)];
end