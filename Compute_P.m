% FUNCTION
% Compute_P.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the generalised weight
%
% SYNOPSIS
% Model = Compute_P(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% OUTPUT
% Model (cf. data structure in user guide)
%
% DESCRIPTION
%
% REFERENCES
% R Dumas, L Cheze. 3D inverse dynamics in non-orthonormal segment 
% coordinate system. Medical & Biological Engineering & Computing 2007;
% 45(3):315-22
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
%__________________________________________________________________________


function Model = Compute_P(Segment,Model)

% Number of frames
n = size(Segment(2).rM,3);

% Initialisation
E33 = eye(3,3);
P = [];

for i = 2:5 % From foot (i = 2) to thigh (i = 5) 
P = [P; ...
    [Segment(i).nC(1,1)*E33,...
    (1 + Segment(i).nC(2,1))*E33,...
    - Segment(i).nC(2,1)*E33,...
    Segment(i).nC(3,1)*E33]'*Segment(i).m*[0;-9.81;0]];
end

% Generalised weight
Model.P = repmat(P,[1,1,n]);
