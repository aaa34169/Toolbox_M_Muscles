% FUNCTION
% Compute_J.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of the pseudo-inertia matrix from the segment mass, 
% inertia and CoM position
%
% SYNOPSIS
% Segment = Compute_J(Segment)
%
% INPUT
% Segment (cf. data structure in user guide)
%
% OUTPUT
% Segment (cf. data structure in user guide)
%
% DESCRIPTION
% Computation of the inertia matrix at proximal endpoint in (u,rP-rD,w)
% using generalised paralell axis theorem 
%
% REFERENCES
% 
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
% Created by Rapha�l Dumas, Florent Moissent, Edouard Jouan
% September 2012
%__________________________________________________________________________


function Segment = Compute_J(Segment)

% Number of frames
n = size(Segment(2).rM,3);

% Initialisation
E33 = eye(3,3);

for i = 2:5 % From foot (i = 2) to thigh (i = 5)
    
        % J = invB*(Is + m*((rCs'*rCs)*E33 - rCs*rCs')*invB'
        Segment(i).J = ...
            Mprod_array3(Mprod_array3(Minv_array3(repmat(Segment(i).B,[1,1,n])),...
            repmat(Segment(i).Is + ...
            Segment(i).m*((Segment(i).rCs'*Segment(i).rCs)*E33 - ...
            Segment(i).rCs*Segment(i).rCs'),[1,1,n])),...
            permute(Minv_array3(repmat(Segment(i).B,[1,1,n])),[2,1,3]));
        % with transpose = permute( ,[2,1,3])
  
end