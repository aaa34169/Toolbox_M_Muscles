% MAIN PROGRAM
% Main_Musculoskeletal_Model.m
%__________________________________________________________________________
%
% PURPOSE
% Building of lower limb model and computation of musculo-tendon, contact,
% ligament and bone forces
%
% SYNOPSIS
% N/A (i.e., main program)
%
% DESCRIPTION
%
% REFERENCE
% F Moissenet, L Cheze, R Dumas. A 3D lower limb musculoskeletal model for
% simultaneous estimation of musculo-tendon, joint contact, ligament and
% bone forces during gait. Journal of Biomechanics 2014;47(1):50-8.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Modify_Segment.m
% Multibody_Optimisation.m
% Compute_J.m
% Compute_G.m
% Compute_P.m
% Compute_R.m
% Compute_L.m
% Static_Optimisation_Lagrange_Multipliers.m
% 
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Raphael Dumas, Florent Moissenet, Edouard Jouan
% September 2012
%
% Modified by Raphael Dumas
% April 2013
% Computation of the non-optimized Lagrange multipiers
% Updated figure
%
% Modified by Raphael Dumas
% February 2014
% Jacobian of the kinematic constraints for spherical joints in SCS and
% full selection of lambda1
%
%__________________________________________________________________________

% -------------------------------------------------------------------------
% LOAD RECORDED DATA
% -------------------------------------------------------------------------
uiopen;

% -------------------------------------------------------------------------
% COMPUTE POSITIONS, ACCELERATIONS AND JACOBIAN MATRIX
% -------------------------------------------------------------------------
% Insert patella as segment 4
Segment = Modify_Segment(Segment);
% Optimisation
[Segment,Joint] = Multibody_Optimisation(Segment,Joint,f);

% % Figure
% Main_Joint_Kinematics

% Model velocities
Model.dQdt = [Segment(2).dQdt; ... % Foot (12*1*n)
    Segment(3).dQdt; ... % Shank (12*1*n)
    Segment(4).dQdt; ... % Patella (12*1*n)
    Segment(5).dQdt]; % Thigh (12*1*n)

% Model accelerations
Model.d2Qdt2 = [Segment(2).d2Qdt2; ... % Foot (12*1*n)
    Segment(3).d2Qdt2; ... % Shank (12*1*n)
    Segment(4).d2Qdt2; ... % Patella (12*1*n)
    Segment(5).d2Qdt2]; % Thigh (12*1*n)

% Model Jacobian
Model.K = [Joint(2).Kk(:,1:48,:); ... % Ankle (5*48*n)
    Joint(3).Kk(:,1:48,:); ... % Tibio-femoral (5*48*n)
    Joint(4).Kk(:,1:48,:); ... % Patello-femoral (6*48*n)
    Joint(5).Kk(:,1:48,:); ... % Hip (3*48*n)
    Segment(2).Kr(:,1:48,:); ... % Foot (6*48*n)
    Segment(3).Kr(:,1:48,:); ... % Shank (6*48*n)
    Segment(4).Kr(:,1:48,:); ... % Patella (6*48*n)
    Segment(5).Kr(:,1:48,:)]; % Thigh (6*48*n)

% -------------------------------------------------------------------------
% PREPARE SEGMENT KINETICS
% -------------------------------------------------------------------------
Segment = Compute_J(Segment); % Pseudo inertia matrix
[Segment,Model] = Compute_G(Segment,Model); % Mass matrix

% -------------------------------------------------------------------------
% COMPUTE EXTERNAL FORCES
% -------------------------------------------------------------------------
Model = Compute_P(Segment,Model); % Weight
Model = Compute_R(Segment,Joint,Model); % Ground reaction forces

% -------------------------------------------------------------------------
% COMPUTE MUSCULAR LEVER ARMS
% -------------------------------------------------------------------------
[Segment,Model] = Compute_L_KleinHorsman_modTLEM(Segment,Model,f,weight);   

% -------------------------------------------------------------------------
% ESTIMATE MUSCULO-TENDON, CONTACT, LIGAMENT AND BONE FORCES USING A 
% CONSTRAINED STATIC OPTIMISATION METHOD
% -------------------------------------------------------------------------
Model = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,weights,'wsm');
Model = Force_Identification(Model);
Model = Dickerson_new(Model,Emg,'KH',n);