% This function calculates the B matrix
% Inputs: 
% b = Angle in radians
function a = B_mat(b)
a = [-sin(b),-cos(b);
      cos(b),-sin(b)];
end