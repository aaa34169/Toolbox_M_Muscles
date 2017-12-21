function J = Criterion_Lagrange_Multipliers(X)

global W;

% Objective function
% Sum of squared forces
J = 1/2*X'*W*X;