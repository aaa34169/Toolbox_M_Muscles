function Segment = Modify_Segment(Segment)

% Number of frames
n = size(Segment(2).rM,3);

% -------------------------------------------------------------------------
% Insert patella as Segment(4) 
% -------------------------------------------------------------------------
Segment(6).Q = Segment(5).Q;
Segment(6).rM = Segment(5).rM;
Segment(5).Q = Segment(4).Q;
Segment(5).m = Segment(4).m;
Segment(5).rCs = Segment(4).rCs;
Segment(5).Is = Segment(4).Is;
Segment(5).rM = Segment(4).rM;
Segment(4).Q = [];                  % to be defined
Segment(4).m = 0.0975;              % according to the OpenSim model
Segment(4).rCs = [0;0;0];           % rCs = rP4;
Segment(4).Is = zeros(3,3);         % patella is considered as a ponctual mass     
Segment(4).rM = [];                 % no associated marker

% -------------------------------------------------------------------------
% Model parameters to define Segment(4)
% -------------------------------------------------------------------------
% Sancisi - 2009 - A new kinematic model of the passive motion of the knee 
% inclusive of the patella
delta_femur = 0.09; % Azimuth of the hinge axis, expressed in the femur SCS
nu_femur = 0.09; % Altitude of the hinge axis, expressed in the femur SCS
toffset = 20;
Roffset = [cosd(toffset) 0 sind(toffset); ...
    0 1 0; ...
    -sind(toffset) 0 cosd(toffset)];
n_axis_femur = Roffset*[sin(nu_femur); ... % Orientation of the hinge axis, expressed in the femur SCS
                cos(nu_femur)*sin(delta_femur); ...
                cos(nu_femur)*cos(delta_femur)];  
P_axis_femur = [4.77; 10.47; -2.68+15]./1000;%[4.77; 10.47; 0]./1000 + ... % Point of the hinge axis, expressed in the femur SCS
               %(-2.68/1000)*n_axis_femur; % Distance between femur and patella points of the hinge axis
L(4) = 21/1000; % Position of patellar tendon origin in patella SCS

% -------------------------------------------------------------------------
% Define Segment(4).Q
% -------------------------------------------------------------------------
% Segment parameter of femur
L(5) = mean(sqrt(sum(Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:)).^2),3);
a(5) = mean(acosd(dot(Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:), Segment(5).Q(10:12,1,:))./...
            sqrt(sum((Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:)).^2))),3);
b(5) = mean(acosd(dot(Segment(5).Q(10:12,1,:), ...
            Segment(5).Q(1:3,1,:))),3);
c(5) = mean(acosd(dot(Segment(5).Q(1:3,1,:), ...
            Segment(5).Q(4:6,1,:) - Segment(5).Q(7:9,1,:))./...
            sqrt(sum((Segment(5).Q(4:6,1,:) - ...
            Segment(5).Q(7:9,1,:)).^2))),3);
B5 = [1, L(5)*cosd(c(5)), cosd(b(5)); ...
        0, L(5)*sind(c(5)), (cosd(a(5)) - cosd(b(5))*cosd(c(5)))/sind(c(5)); ...
        0, 0, sqrt(1 - cosd(b(5))^2 - ((cosd(a(5)) - cosd(b(5))*cosd(c(5)))/sind(c(5)))^2)];
invB5 = inv(B5);

% Interpolation matrices
% Orientation of the Z axis of femur SCS in the NSCS
nZ5 = invB5(:,3); 
NVZ5 =[nZ5(1,1)*eye(3),...
       (nZ5(2,1))*eye(3), ...
       - nZ5(2,1)*eye(3), ...
       nZ5(3,1)*eye(3)]; 
% V65: virtual marker 6 of segment 5
% Point of the hinge axis
Segment(5).nV(:,6) = [0;-1;0] + invB5*P_axis_femur;
NV65 = [Segment(5).nV(1,6)*eye(3),...
        (1 + Segment(5).nV(2,6))*eye(3),...
        - Segment(5).nV(2,6)*eye(3), ...
        Segment(5).nV(3,6)*eye(3)];
    
% Segment parameters
u4 = Vnorm_array3((Segment(3).Q(1:3,1,:) + Segment(5).Q(1:3,1,:))/2);           % Mean of thigh and shank u axes
v4 = Vnorm_array3(cross(Mprod_array3(repmat(NVZ5,[1,1,n]),Segment(5).Q),u4));   % Z axis of femur SCS in the NSCS
w4 = Vnorm_array3(cross(u4,v4));
rP4 = Mprod_array3(repmat(NV65,[1,1,n]),Segment(5).Q) + (42.48/1000)*u4;        % Distance from the point of the hinge axis to the origin of patella SCS 
rD4 = rP4 - L(4)*v4;
Segment(4).Q = [u4;rP4;rD4;w4];