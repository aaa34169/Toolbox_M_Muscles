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

% grandChallenge = 1;
% filename = {'jw_ngait_2' 'jw_ngait_3' 'jw_ngait_4' 'jw_ngait_5' 'jw_ngait_6'};
% weights = [1e0;1e0;1e0;1e-6;1e-6;1;5;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6];
grandChallenge = 2;
filename = {'dm_ngait4' 'dm_ngait10' 'dm_ngait11' 'dm_ngait12' 'dm_ngait13'};
weights = [1e0;1e0;1e0;1e-6;1e-6;2;1;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6];
% grandChallenge = 3;
% filename = {'SC_ngait_og5' 'SC_ngait_og6' 'SC_ngait_og7' 'SC_ngait_og8' 'SC_ngait_og9'};
% weights = [1e0;1e0;1e0;1e-6;1e-6;10;30;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6];
% grandChallenge = 4; % not used for abstract
% filename = {'jw_ngait_og1' 'jw_ngait_og2' 'jw_ngait_og3' 'jw_ngait_og4' 'jw_ngait_og7'};
% weights = [1e0;1e0;1e0;1e-6;1e-6;15;15;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6];
% grandChallenge = 5;
% filename = {'PS_ngait_og_ss1' 'PS_ngait_og_ss7' 'PS_ngait_og_ss8' 'PS_ngait_og_ss9' 'PS_ngait_og_ss11'};
% weights = [1e0;1e0;1e0;1e-6;1e-6;4;4;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6];
% grandChallenge = 6; % not used for abstract
% filename = {'DM_ngait_og1' 'DM_ngait_og3' 'DM_ngait_og4' 'DM_ngait_og5' 'DM_ngait_og6'};
% weights = [1e0;1e0;1e0;1e-6;1e-6;10;0.1;1e-6;1e-6;1e-6;1e0;1e0;1e0;1e-6;1e0;1e0;1e0;1e-6;1e-6;1e-6;1e-6];

% for i = 1:length(filename)
% % -------------------------------------------------------------------------
% % LOAD RECORDED DATA
% % -------------------------------------------------------------------------
% clearvars -except grandChallenge filename i weights;
% cd('C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2018\IMSD\2014-2015');
% addpath('C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2018\IMSD\data\');
% load(['C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2018\IMSD\data\grand_challenge_',num2str(grandChallenge),'\',filename{i},'_pPS.mat']);
% 
% % -------------------------------------------------------------------------
% % COMPUTE POSITIONS, ACCELERATIONS AND JACOBIAN MATRIX
% % -------------------------------------------------------------------------
% 
% % Insert patella as segment 4
% Segment = Modify_Segment(Segment);
% % Optimisation
% [Segment,Joint] = Multibody_Optimisation(Segment,Joint,f);
% 
% % % Figure
% % Main_Joint_Kinematics
% 
% % Model velocities
% Model.dQdt = [Segment(2).dQdt; ... % Foot (12*1*n)
%     Segment(3).dQdt; ... % Shank (12*1*n)
%     Segment(4).dQdt; ... % Patella (12*1*n)
%     Segment(5).dQdt]; % Thigh (12*1*n)
% 
% % Model accelerations
% Model.d2Qdt2 = [Segment(2).d2Qdt2; ... % Foot (12*1*n)
%     Segment(3).d2Qdt2; ... % Shank (12*1*n)
%     Segment(4).d2Qdt2; ... % Patella (12*1*n)
%     Segment(5).d2Qdt2]; % Thigh (12*1*n)
% 
% % Model Jacobian
% Model.K = [Joint(2).Kk(:,1:48,:); ... % Ankle (5*48*n)
%     Joint(3).Kk(:,1:48,:); ... % Tibio-femoral (5*48*n)
%     Joint(4).Kk(:,1:48,:); ... % Patello-femoral (6*48*n)
%     Joint(5).Kk(:,1:48,:); ... % Hip (3*48*n)
%     Segment(2).Kr(:,1:48,:); ... % Foot (6*48*n)
%     Segment(3).Kr(:,1:48,:); ... % Shank (6*48*n)
%     Segment(4).Kr(:,1:48,:); ... % Patella (6*48*n)
%     Segment(5).Kr(:,1:48,:)]; % Thigh (6*48*n)
% 
% % -------------------------------------------------------------------------
% % PREPARE SEGMENT KINETICS
% % -------------------------------------------------------------------------
% Segment = Compute_J(Segment); % Pseudo inertia matrix
% [Segment,Model] = Compute_G(Segment,Model); % Mass matrix
% 
% % -------------------------------------------------------------------------
% % COMPUTE EXTERNAL FORCES
% % -------------------------------------------------------------------------
% Model = Compute_P(Segment,Model); % Weight
% Model = Compute_R(Segment,Joint,Model); % Ground reaction forces
% 
% % -------------------------------------------------------------------------
% % COMPUTE MUSCULAR LEVER ARMS
% % -------------------------------------------------------------------------
% [Segment,Model] = Compute_L_KleinHorsman_modTLEM(Segment,Model,f,weight);   
% 
% % -------------------------------------------------------------------------
% % ESTIMATE MUSCULO-TENDON, CONTACT, LIGAMENT AND BONE FORCES USING A 
% % CONSTRAINED STATIC OPTIMISATION METHOD
% % -------------------------------------------------------------------------
% Model1 = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,weights,'wsm');
% Model1 = Force_Identification(Model1);
% Model1 = Dickerson_new(Model1,Emg,'KH',n);
% Model2 = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,weights,'mmm');
% Model2 = Force_Identification(Model2);
% Model2 = Dickerson_new(Model2,Emg,'KH',n);
% 
% save(['C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2018\IMSD\data\grand_challenge_',num2str(grandChallenge),'\',filename{i},'_results.mat']);
% end

outcomes = processing(grandChallenge,filename);
save(['C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2018\IMSD\data\grand_challenge_',num2str(grandChallenge),'\goodnessOfFit_results.mat']);
saveas(gca, ['C:\Users\florent.moissenet\Documents\Professionnel\publications\communications\2018\IMSD\results\grand_challenge_',num2str(grandChallenge),'_plot'], 'png');