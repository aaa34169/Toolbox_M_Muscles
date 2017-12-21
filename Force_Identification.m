% FUNCTION
% Force_Identification.m
%__________________________________________________________________________
%
% PURPOSE
% Identify musculo-tendon, contact, ligament and bone forces
%
% SYNOPSIS
% Model = Force_Identification(Segment,Joint,Model)
%
% INPUT
% Model (cf. data structure in user guide)
%
% OUTPUT
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Identify musculo-tendon, contact, ligament and bone forces and save them
% in different substructures
%
% REFERENCE
% 
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% 
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Rapha�l Dumas, Florent Moissent, Edouard Jouan
% March 2015
%__________________________________________________________________________

function Model = Force_Identification(Model)

% Number of frames
n = size(Model.X,3);

% Number of muscles
m = size(Model.Lever,2);

% Musculo-tendon forces (tensile forces > 0)
Model.Fm = Model.X(1:m,1,1:n); % About lines of action

% Contact forces (reaction forces > 0 i.e., acting on the proximal segment)
Model.Fc = [Model.X(m+1:m+3,1,1:n); ... % 3D ankle contact in foot SCS
    Model.X(m+6,1,1:n);Model.X(m+7,1,1:n); ... % Tibio-femoral medial and lateral contact in shank SCS
    Model.X(m+11:m+13,1,1:n); ... % 3D patello-femoral contact in patella SCS
    Model.X(m+15:m+17,1,1:n)]; % 3D hip contact in thigh SCS

% Ligament forces (tensile forces > 0)
% About lines of action
Model.Fl = [Model.X(m+4:m+5,1,1:n); ... % TiCaL, CaFil
    Model.X(m+8:m+10,1,1:n); ... % ACL, PCL, MCL
    Model.X(m+14,1,1:n)]; % PT

% Bone forces (compression forces > 0)
% About segment Y axis
Model.Fb = Model.X(m+18:m+21,1,1:n); % Foot, tibia, patella, femur axial