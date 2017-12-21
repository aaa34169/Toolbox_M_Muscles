% FUNCTION
% Compute_P.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the generalised ground reaction force
%
% SYNOPSIS
% Model = Compute_R(Segment,Model)
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
% Mprod_array3.m
% Minv_array3.m
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
%__________________________________________________________________________


function Model = Compute_R(Segment,Joint,Model)

% Number of frames
n = size(Segment(2).rM,3);
% Initialisation
Model.R = zeros(48,1,n);
M33 = zeros(3,3,n);

% Imput data
F1 = Joint(1).F; % Force of the foot on ground
M1 = Joint(1).M; % Moment of the foot on ground at the CoP
rP1 = Segment(1).Q(4:6,:,:); % rP1 = CoP

% Transpose of the interpolation matrix of rP2
NP2t = [zeros(3,3,n);...
    repmat(eye(3,3),[1,1,n]);...
    zeros(3,3,n);...
    zeros(3,3,n)];

% Segment parameters
u2 = Segment(2).Q(1:3,1,:);
rP2 = Segment(2).Q(4:6,1,:);
rD2 = Segment(2).Q(7:9,1,:);
w2 = Segment(2).Q(10:12,1,:);

% Composition of interpolation matrices and lever arms
M33(1:3,2,:) = rP2 - rD2;
M33(4:6,3,:) = -w2;
M33(7:9,3,:) = w2;
M33(10:12,1,:) = u2;

% Bstar
Bstar2(1:3,1,:) = cross(w2,u2);
Bstar2(1:3,2,:) = cross(u2,(rP2 - rD2));
Bstar2(1:3,3,:) = cross(-(rP2 - rD2),w2);

% Nstar2t (pseudo interpolation matrix)
Nstar2t = Mprod_array3(M33,Minv_array3(Bstar2));

% Generalised ground reaction forces
Model.R = [Mprod_array3(NP2t,-F1) + ...
    Mprod_array3(Nstar2t,(-M1 + cross((rP1 - rP2),-F1)));...
    zeros(12,1,n);...
    zeros(12,1,n);...
    zeros(12,1,n)];
          