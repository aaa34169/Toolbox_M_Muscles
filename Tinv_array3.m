% FUNCTION
% Tinv_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of homogenous matrix inverse 
%
% SYNOPSIS
% B = Tinv_array3(A)
%
% INPUT
% A (i.e., homogenous matrix) 
%
% OUTPUT
% B (i.e., homogenous matrix)
%
% DESCRIPTION
% Computation, for all frames (i.e., in 3rd dimension, cf. data structure
% in user guide), of the inverse of an homogenous matrix (T)
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% Mprod_array3.m
% 
% MATLAB VERSION
% Matlab R2007b
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphaël Dumas
% March 2010
%__________________________________________________________________________

function B = Tinv_array3(A)

% Initialization
A1 = [];
A2 = [];

% First matrix
% Transposition of rotational part
% with transpose = permute( ,[2,1,3])
A1(1:3,1:3,:) = permute(A(1:3,1:3,:),[2,1,3]);
A1(4,4,:) = 1;

% Second matrix
% Opposite sign of the translation part
A2(1:4,4,:) = - A(1:4,4,:);
A2(1,1,:) = 1;
A2(2,2,:) = 1;
A2(3,3,:) = 1;
A2(4,4,:) = 1;

% Matrix product
B = Mprod_array3(A1,A2);
