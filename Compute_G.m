% FUNCTION
% Compute_G.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the generalised mass matrix
%
% SYNOPSIS
% [Segment,Model] = Compute_G(Segment,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% OUTPUT
% Segment (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Computation of the segment and model mass matrices 
%
% REFERENCES
% R Dumas, L Cheze. 3D inverse dynamics in non-orthonormal segment 
% coordinate system. Medical & Biological Engineering & Computing 2007;
% 45(3):315-22
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX)
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas, Florent Moissent, Edouard Jouan
% September 2012
%__________________________________________________________________________


function [Segment,Model] = Compute_G(Segment,Model)

% Number of frames
n = size(Segment(2).rM,3);

% Initialmisation
E33 = eye(3,3);
G = zeros(48,48);

for i = 2:5 % From foot (i = 2) to thigh (i = 5)
    
    % CoM Coordinates in generalized coordinates (u,rP-rD,w)
    % ---------------------------------------------------------------------
    Segment(i).nC = inv(Segment(i).B)*Segment(i).rCs;
    
    % 1st line
    % ---------------------------------------------------------------------
    Segment(i).G(1:3,1:3) = Segment(i).J(1,1)*E33;
    Segment(i).G(1:3,4:6) = (Segment(i).m*Segment(i).nC(1,1) + Segment(i).J(1,2))*E33;
    Segment(i).G(1:3,7:9) = - Segment(i).J(1,2)*E33;
    Segment(i).G(1:3,10:12) = Segment(i).J(1,3)*E33;
    
    % 2nd line
    % ---------------------------------------------------------------------
    Segment(i).G(4:6,1:3) = (Segment(i).m*Segment(i).nC(1,1) + Segment(i).J(1,2))*E33;
    Segment(i).G(4:6,4:6) = (Segment(i).m + 2*Segment(i).m*Segment(i).nC(2,1) + Segment(i).J(2,2))*E33;
    Segment(i).G(4:6,7:9) = (- Segment(i).m*Segment(i).nC(2,1) - Segment(i).J(2,2))*E33;
    Segment(i).G(4:6,10:12) = (Segment(i).m*Segment(i).nC(3,1) + Segment(i).J(2,3))*E33;
    
    % 3rd line
    % ---------------------------------------------------------------------
    Segment(i).G(7:9,1:3) = - Segment(i).J(1,2)*E33;
    Segment(i).G(7:9,4:6) = (- Segment(i).m*Segment(i).nC(2,1) - Segment(i).J(2,2))*E33;
    Segment(i).G(7:9,7:9) = Segment(i).J(2,2)*E33;
    Segment(i).G(7:9,10:12) = - Segment(i).J(2,3)*E33;
    
    % 4th line
    % ---------------------------------------------------------------------
    Segment(i).G(10:12,1:3) = Segment(i).J(1,3)*E33;
    Segment(i).G(10:12,4:6) = (Segment(i).m*Segment(i).nC(3,1) + Segment(i).J(2,3))*E33;
    Segment(i).G(10:12,7:9) = - Segment(i).J(2,3)*E33;
    Segment(i).G(10:12,10:12) = Segment(i).J(3,3)*E33;
    
    % Creation of the generalized mass matrix (48x48)
    % ---------------------------------------------------------------------
    s = 12*(i-2);
    for j = 1:12
        for k = 1:12
            G(s+j,s+k) = Segment(i).G(j,k);
        end
    end
    
end

% Generalized mass matrix 
Model.G = repmat(G,[1,1,n]);