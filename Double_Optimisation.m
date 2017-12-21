% FUNCTION
% Double_Optimisation.m
%__________________________________________________________________________
%
% PURPOSE
% Computation of optimised weights based on RMSE difference on tibiofemoral
% forces
%
% SYNOPSIS
% w = Double_Optimisation(Segment,Joint,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Joint (cf. data structure in user guide)
% Model (cf. data structure in user guide)
%
% OUTPUT
% w (diagonal items of the weight matrix - of the static optimisation - 
% corresponding to the weights of the contact/ligament/bone forces)
%
% DESCRIPTION
% Find the minimum of the sum of squared forces suject to equality 
% (i.e., dynamic equlibrium) and inequality constraints (i.e., positive
% forces
%
% REFERENCES
% 
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
%
% MATLAB VERSION
% Matlab R2011b Unix
%__________________________________________________________________________
%
% CHANGELOG
% Created by Florent Moissenet
% February 2013
%__________________________________________________________________________

function w = Double_Optimisation(Segment,Joint,Contact,Emg,Geometry,f,weight)

F_med = permute(Contact.KneeMedial,[3,1,2]);
F_lat = permute(Contact.KneeLateral,[3,1,2]);

figure(1); hold on;
plot(F_med+F_lat,'red');
% figure(2); hold on;
% plot(F_lat,'red');

% Number of frames
n = size(Segment(2).rM,3);

h = waitbar(1,'Optimization is running...');

Segment = Modify_Segment(Segment);
[Segment,Joint] = Multibody_Optimisation(Segment,Joint,f);
Model.dQdt = [Segment(2).dQdt; ... % Foot (12*1*n)
    Segment(3).dQdt; ... % Shank (12*1*n)
    Segment(4).dQdt; ... % Patella (12*1*n)
    Segment(5).dQdt]; % Thigh (12*1*n)
Model.d2Qdt2 = [Segment(2).d2Qdt2; ... % Foot (12*1*n)
    Segment(3).d2Qdt2; ... % Shank (12*1*n)
    Segment(4).d2Qdt2; ... % Patella (12*1*n)
    Segment(5).d2Qdt2]; % Thigh (12*1*n)
Model.K = [Joint(2).Kk(:,1:48,:); ... % Ankle (5*48*n)
    Joint(3).Kk(:,1:48,:); ... % Tibio-femoral (5*48*n)
    Joint(4).Kk(:,1:48,:); ... % Patello-femoral (6*48*n)
    Joint(5).Kk(:,1:48,:); ... % Hip (3*48*n)
    Segment(2).Kr(:,1:48,:); ... % Foot (6*48*n)
    Segment(3).Kr(:,1:48,:); ... % Shank (6*48*n)
    Segment(4).Kr(:,1:48,:); ... % Patella (6*48*n)
    Segment(5).Kr(:,1:48,:)]; % Thigh (6*48*n)
Segment = Compute_J(Segment); % Pseudo inertia matrix
[Segment,Model] = Compute_G(Segment,Model); % Mass matrix
Model = Compute_P(Segment,Model); % Weight
Model = Compute_R(Segment,Joint,Model); % Ground reaction forces
if strcmp(Geometry,'Delp')
    [Segment,Model] = Compute_L_Delp(Segment,Model,f,weight);
elseif strcmp(Geometry,'KH')
    [Segment,Model] = Compute_L_KleinHorsman(Segment,Model,f,weight); 
end

% Optimisation 2 (weighted sum method - condition C2)

% Included contact, ligament and bone forces
Model.K1(1,:,:) = Model.K(1,1:48,:) ... % Ankle contact about X axis
    *(-1); % Coefficient for contact force;
Model.K1(2,:,:) = Model.K(2,1:48,:) ... % Ankle contact about Y axis
    *(-1); % Coefficient for contact force
Model.K1(3,:,:) = Model.K(3,1:48,:) ... % Ankle contact about Z axis
    *(-1); % Coefficient for contact force;
Model.K1(4,:,:) = Model.K(4,1:48,:) ... % TiCaL
    *(1/(2*Joint(2).d(1,1))); % Coefficient for ligament force
Model.K1(5,:,:) = Model.K(5,1:48,:) ... % CaFiL
    *(1/(2*Joint(2).d(1,2))); % Coefficient for ligament force
Model.K1(6,:,:) = Model.K(6,1:48,:) ... % Knee medial contact
    *(-1); % Coefficient for contact force
Model.K1(7,:,:) = Model.K(7,1:48,:) ... % Knee lateral contact
    *(-1); % Coefficient for contact force
Model.K1(8,:,:) = Model.K(8,1:48,:) ... % ACL
    *(1/(2*Joint(3).d(1,3))); % Coefficient for ligament force
Model.K1(9,:,:) = Model.K(9,1:48,:) ... % PCL
    *(1/(2*Joint(3).d(1,4))); % Coefficient for ligament force
Model.K1(10,:,:) = Model.K(10,1:48,:) ... % MCL
    *(1/(2*Joint(3).d(1,4))); % Coefficient for ligament force;
Model.K1(11,:,:) = Model.K(11,1:48,:) ... % Patellar contact about X axis
    *(-1); % Coefficient for contact force
Model.K1(12,:,:) = Model.K(12,1:48,:) ... % Patellar contact about Y axis
    *(-1); % Coefficient for contact force
Model.K1(13,:,:) = Model.K(13,1:48,:) ... % Patellar contact about Z axis
    *(-1); % Coefficient for contact force
Model.K1(14,:,:) = Model.K(16,1:48,:) ... % PT
    *(1/(2*Joint(4).d(1,1))); % Coefficient for ligament force
Model.K1(15,:,:) = Model.K(17,1:48,:) ... % Hip contact about X axis
    *(-1); % Coefficient for contact force
Model.K1(16,:,:) = Model.K(18,1:48,:) ... % Hip contact about Y axis
    *(-1); % Coefficient for contact force
Model.K1(17,:,:) = Model.K(19,1:48,:) ...% Hip contact about Z axis
    *(-1); % Coefficient for contact force
Model.K1(18,:,:) = Model.K(23,1:48,:) ... % Foot axial
    *(-1/(2*Segment(2).L)); % Coefficient for bone force
Model.K1(19,:,:) = Model.K(29,1:48,:) ... % Tibia axial
    *(-1/(2*Segment(3).L)); % Coefficient for bone force
Model.K1(20,:,:) = Model.K(35,1:48,:) ... % Patella axial
    *(-1/(2*Segment(4).L)); % Coefficient for bone force
Model.K1(21,:,:) = Model.K(41,1:48,:) ... % Femur axial
    *(-1/(2*Segment(5).L)); % Coefficient for bone force

% Discarded forces
% Coefficients are introduced in order to associate Lagrange multiplieurs
% with external forces
Model.K2(1,:,:) = Model.K(14,1:48,:);
Model.K2(2,:,:) = Model.K(15,1:48,:);
Model.K2(3,:,:) = Model.K(20,1:48,:);
Model.K2(4,:,:) = Model.K(21,1:48,:);
Model.K2(5,:,:) = Model.K(22,1:48,:);
Model.K2(6,:,:) = Model.K(24,1:48,:);
Model.K2(7,:,:) = Model.K(25,1:48,:);
Model.K2(8,:,:) = Model.K(26,1:48,:);
Model.K2(9,:,:) = Model.K(27,1:48,:);
Model.K2(10,:,:) = Model.K(28,1:48,:);
Model.K2(11,:,:) = Model.K(30,1:48,:);
Model.K2(12,:,:) = Model.K(31,1:48,:);
Model.K2(13,:,:) = Model.K(32,1:48,:);
Model.K2(14,:,:) = Model.K(33,1:48,:);
Model.K2(15,:,:) = Model.K(34,1:48,:);
Model.K2(16,:,:) = Model.K(36,1:48,:);
Model.K2(17,:,:) = Model.K(37,1:48,:);
Model.K2(18,:,:) = Model.K(38,1:48,:);
Model.K2(19,:,:) = Model.K(39,1:48,:);
Model.K2(20,:,:) = Model.K(40,1:48,:);
Model.K2(21,:,:) = Model.K(42,1:48,:);
Model.K2(22,:,:) = Model.K(43,1:48,:);

% -------------------------------------------------------------------------
% RUN DOUBLE OPTIMISATION TO GET AN OPTIMIZED WEIGHT MATRIX
% -------------------------------------------------------------------------

% Initial guess, min and max boundaries
% -------------------------------------------------------------------------
w_max = 1e3*ones(21,1);
w_min = zeros(21,1);
w_ini = 10*ones(21,1);

% Optimisation (upper level)
% -------------------------------------------------------------------------
options = optimset('Display','off','TolFun',0.1,'TolX',1e-10,'TolCon',1e-10,...
    'algorithm','active-set','DiffMinChange',1e0);
w = fmincon(@(w) Criterion_upper_level(w,Joint,Model,Contact,Emg,Geometry,n,h),w_ini,[],[],[],[],w_min,w_max,[],options);

% -------------------------------------------------------------------------
% UPPER LEVEL CRITERION
% -------------------------------------------------------------------------
function J = Criterion_upper_level(w,Joint,Model,Contact,Emg,Geometry,n,h)

global W

% Number of muscles
if strcmp(Geometry,'Delp')
    m = 43;
elseif strcmp(Geometry,'KH')
    m = 129;
end

W = eye(m+21); % Number of muscles + number of Lagrange multipliers
weights = [w(1);w(2);w(3);w(4);w(5);w(6);w(7);w(8);w(9);w(10);w(11);w(12);w(13);w(14);w(15);w(16);w(17);w(18);w(19);w(20);w(21)];
for i = 1:21
    W(m+i,m+i) = weights(i);
end

for i = 1:n
    
    waitbar(i/n,h);
    
    % ---------------------------------------------------------------------
    % Partial parameter reduction and constraint equations
    % ---------------------------------------------------------------------
    [eigvector,~] = eig(Model.K2(:,:,i)'*Model.K2(:,:,i));
    Model.ZK2(:,:,i) = eigvector(:,1:26); % 26 first eigenvalues are 0
    Aeq = Model.ZK2(:,:,i)'*[Model.Lever(:,:,i), - Model.K1(:,:,i)'];
    Beq = Model.ZK2(:,:,i)'*...
        (Model.G(:,:,i)*Model.d2Qdt2(:,:,i) - Model.P(:,:,i) - Model.R(:,:,i));
    
    % Initial guess
    if i == 1
        Xini = zeros(m+21,1);
    else
        Xini = Model.X(:,1,i-1);
    end
    
    % Lower abounds
    Xmin = zeros(m+21,1);
    
    Xmin(m+1,1) = -Inf; % Ankle contact about X axis
    % Ankle contact about Y axis
    Xmin(m+3,1) = -Inf; % Ankle contact about Z axis
    % TiCaL
    % CaFiL
    % Knee medial contact
    % Knee lateral contact
    % ACL
    % PCL
    Xmin(m+10,1) = -Inf; % MCL
    Xmin(m+11,1) = -Inf; % Patellar contact about X axis
    Xmin(m+12,1) = -Inf; % Patellar contact about Y axis
    Xmin(m+13,1) = -Inf; % Patellar contact about Z axis
    % PT
    Xmin(m+15,1) = -Inf; % Hip contact about X axis
    % Hip contact about Y axis
    Xmin(m+17,1) = -Inf; % Hip contact about Z axis
    Xmin(m+18,1) = -Inf; % Foot axial
    % Tibia axial
    Xmin(m+20,1) = -Inf; % Patella axial
    % Femur axial
    
    % Upper bounds
    Xmax = Inf(m+21,1);
%    Xmax = [Model.Fmax(1:m,1);Inf(21,1)];
%    Xmax = [Model.PCSA(1:m,1)*61;Inf(21,1)]; % Muscle stress of 61 N/cm^2
    
    % ---------------------------------------------------------------------
    % Run optimization (fmincon)
    % ---------------------------------------------------------------------
    options = optimset('Display','off','MaxIter',100,'LargeScale','off','TolFun',0.1,'TolX',1e-10,'TolCon',1e-10,...
        'algorithm','active-set','GradObj','off');
    [X,~,exitflag] = fmincon(@Criterion_lower_level,Xini,[],[],Aeq,Beq,Xmin,Xmax,[],options);
    Model.X(:,1,i) = X;
%     flag(i) = exitflag;
    Model = Force_Identification(Model);
    
end

% Objective function (RMSE between predicted tibiofemoral forces and
% Grand Challenge implant measurements)
% -------------------------------------------------------------------------
F_med = permute(Contact.KneeMedial,[3,1,2]);
F_lat = permute(Contact.KneeLateral,[3,1,2]);

Fe_med = permute(Model.Fc(4,:,:),[3,1,2]);
Fe_lat = permute(Model.Fc(5,:,:),[3,1,2]);

% Model = Dickerson_new(Model,Emg,Geometry,n);

% J = (1/Model.concordance_EMG)*1e5*1e1 + ...
%     (1/Model.concordance_Perry)*1e5 + ...
%     sqrt(sum((Fe_med-F_med).^2)/(numel(F_med))) + ...
%     sqrt(sum((Fe_lat-F_lat).^2)/(numel(F_lat)))
J = sqrt(sum(((Fe_med+Fe_lat)-(F_med+F_lat)).^2)/(numel((F_med+F_lat))))
figure(1);
plot(Fe_med+Fe_lat);

% -------------------------------------------------------------------------
% LOWER LEVEL CRITERION
% -------------------------------------------------------------------------
function J = Criterion_lower_level(X)

global W;

% Objective function
% Sum of squared forces
J = 1/2*X'*W*X;
