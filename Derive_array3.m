% FUNCTION
% Derive_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of time derivative of matrix or vector 
%
% SYNOPSIS
% dVdt = Derive_array3(V,dt)
% dMdt = Derive_array3(M,dt)
%
% INPUT
% V (i.e., vector) or M (i.e., matrix) 
% dt (i.e., sampling time that is to say dt = 1/f)
%
% OUTPUT
% dVdt or dMdt (i.e., vector or matrix)
%
% DESCRIPTION
% Gradient approximation in the 3rd dimension of matrix or vector 
% (i.e., all frames, cf. data structure in user guide)
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% None
% 
% MATLAB VERSION
% Matlab R2007b
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas
% March 2010
%__________________________________________________________________________

function dVdt = Derive_array3(V,dt)

% Trick for passing both versions R14 and R2007
if size(V,2) == 1
    [gr1,gr2,dVdt] = gradient(repmat(V,[1,size(V,1),1]),dt);
    dVdt = dVdt(1:size(V,1),1,:);
else
    % 3rd direction gradient
    % with time sampling (dt)
    [gr1,gr2,dVdt] = gradient(V,dt);
end
