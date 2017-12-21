% FUNCTION
% Static_Optimisation_Lagrange_Multipliers.m
%__________________________________________________________________________
%
% PURPOSE
% Computation musculo-tendon, contact, ligament and bone forces
%
% SYNOPSIS
% Model = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model)
%
% INPUT
% Segment (cf. data structure in user guide)
% Joint (cf. data structure in user guide)
% Model (cf. data structure in user guide)
% weights (cf. data structure in user guide)
% method (cf. data structure in user guide)
%
% OUTPUT
% Model (cf. data structure in user guide)
%
% DESCRIPTION
% Find the minimum of the sum of squared forces suject to equality 
% (i.e., dynamic equlibrium) and inequality constraints (i.e., positive
% forces
%
% REFERENCE
% F Moissenet, L Cheze, R Dumas. A 3D lower limb musculoskeletal model for
% simultaneous estimation of musculo-tendon, joint contact, ligament and
% bone forces during gait. Journal of Biomechanics 2014;47(1):50-8.
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Criterion_Lagrange_Multipliers.m
%
% MATLAB VERSION
% Matlab R2012a
%__________________________________________________________________________
%
% CHANGELOG
% Created by Rapha�l Dumas, Florent Moissent, Edouard Jouan
% September 2012
%
% Modified by Raphael Dumas
% February 2014
% Full selection of lambda1
%
% Modified by Florent Moissenet
% March 2015
% Possibility to choose between weighted sum method or min-max method
%__________________________________________________________________________

function Model = Static_Optimisation_Lagrange_Multipliers(Segment,Joint,Model,weights,method)

global W;

% Number of frames
n = size(Model.K,3);

% Number of muscles
m = size(Model.Lever,2);

% Waitbar
h = waitbar(1,'Optimization is running...');

% -------------------------------------------------------------------------
% Prepare reduced Jacobian matrices
% K1 is kept, K2 is removed
%
% Coefficients are introduced in order to associate Lagrange multiplieurs
% with external forces acting on the proximal segment in order to have 
% positive force about the SCS axes (or axes of the structure)
% Contact forces: reaction forces
% Ligament forces: tensile forces
% Bone forces: compression forces
%
% -------------------------------------------------------------------------

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
% Weight matrix
% -------------------------------------------------------------------------
if strcmp(method,'wsm')
    W = eye(m+21); % Number of muscles + number of Lagrange multipliers
    for i = 1:m
        W(i,i) = (1/m);
    end
    for i = 1:21
        W(m+i,m+i) = weights(i);
    end
end

% -------------------------------------------------------------------------
% RUN OPTIMISATION FOR EACH FRAME
% -------------------------------------------------------------------------
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
    
    % Lower bounds
    Xmin = zeros(m+21,1);    
    Xmin(m+1,1) = -Inf; % Ankle contact about X axis
    % Ankle contact about Y axis
    Xmin(m+3,1) = -Inf; % Ankle contact about Z axis
%     Xmin(m+4,1) = -Inf; % TiCaL
%     Xmin(m+5,1) = -Inf; % CaFiL
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
    Xmax = [Model.PCSA(1:m,1)*61;Inf(21,1)]; % Muscle stress of 61 N/cm^2 (Arnold et al. 2011)
    
    % ---------------------------------------------------------------------
    % Run optimization (fmincon)
    % ---------------------------------------------------------------------
    if strcmp(method,'wsm')
        options = optimset('Display','off','MaxIter',100,'LargeScale','off','TolFun',0.1,'TolX',1e-10,'TolCon',1e-10,...
            'algorithm','active-set','GradObj','off');
        [X,~,exitflag] = fmincon(@Criterion_Lagrange_Multipliers,Xini,[],[],Aeq,Beq,Xmin,Xmax,[],options);
    elseif strcmp(method,'mmm')
        options = optimset('Display','off','MaxIter',100,'LargeScale','off','TolFun',0.1,'TolX',1e-10,'TolCon',1e-10,...
            'GradObj','off');
        [X,~,~,exitflag] = fgoalattain(@(X) Criterion_Lagrange_Multipliers_Minmax(X,m),Xini,zeros(3,1),ones(3,1),[],[],Aeq,Beq,Xmin,Xmax,[],options);
    end
    Model.X(:,1,i) = X;
    flag(i) = exitflag;
    
end

close(h);
figure
plot(flag);