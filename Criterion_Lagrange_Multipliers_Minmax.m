function J = Criterion_Lagrange_Multipliers_Minmax(X,m)

% J1: musculo-tendon forces
% J2: joint contact forces
% J3: ligament forces

x = 1:m;
J(1) = 1/2*(1/(length(x)^2))*((X(x,1))'*X(x,1));

x = m+[1;2;3;6;7;11;12;13;15;16;17];
J(2) = 1/2*(1/(length(x)^2))*((X(x,1))'*X(x,1));

x = m+[4;5;8;9;10;14];
J(3) = 1/2*(1/(length(x)^2))*((X(x,1))'*X(x,1));