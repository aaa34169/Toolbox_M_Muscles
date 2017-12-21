
% FUNCTION
%Validation_model.m
%__________________________________________________________________________
%
% PURPOSE
% validation of optimization results with EMG
%
% SYNOPSIS
% Model=Validation_model(Model,Emg,Force,n,weight)
%
% INPUT
% Model (cf. data structure in user guide)
% Emg (electromyograms of 14 muscles)
% Force (force of the knee contact mesured with a instrumented prothesis)
% n (number of frames)
% weight (mass of the subject)
%
% OUTPUT
% Model (cf. data structure in user guide)
%
% DESCRIPTIONEmg
% Validation based on the Dikerson's method (for muscles) couple with a RMS
% (for knee contacts)
% 
% REFERENCES
% Dikerson 08
% Perry 92
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM MUSCULO-SKELETAL TOOLBOX)
% Criterion_Pseudoinverse.m
%
% MATLAB VERSION
% Matlab R2007b
%__________________________________________________________________________
%
% CHANGELOG
% Created by Rapha�l Dumas, Florent Moissent, Edouard Jouan ,Matthieu
% Giroux
% January 2013
%__________________________________________________________________________


function Model = Dickerson(Joint,Model,Emg,n)

% -------------------------------------------------------------------------
% PHASE definition
% -------------------------------------------------------------------------
first = find(Joint(1,1).F(2,1,:),1,'first'); % frame initial contact
last = find(Joint(1,1).F(2,1,:),1,'last');   % frame between pre-swing and initial swing
swing = n-last;

Phase1 = first:ceil(0.17*last);                           % loading responce
Phase2 = ceil(0.17*last+1):ceil(0.5*last);                % midstance
Phase3 = ceil(0.5*last+1):ceil(0.83*last);                % terminal stance
Phase4 = ceil(0.83*last+1):last;                          % pre-Swing
Phase5 = ceil(last+1):ceil(0.32*swing+last);              % initial Swing
Phase6 = ceil(0.32*swing+last+1):ceil(0.67*swing+last);   % midswing
Phase7 = ceil(0.67*swing+last+1):n;                       % terminal Swing

% -------------------------------------------------------------------------
% EMG  Data
% -------------------------------------------------------------------------
EMG = zeros(14,n);
MAX_EMG = zeros(14,1);
MIN_EMG = zeros(14,1);

EMG(1,:) = Emg.gmax(:)';        % gluteus_maximus
EMG(2,:) = Emg.gmed(:)';        % gluteus_medius
EMG(3,:) = Emg.addmagnus(:)';   % adductus_magnus
EMG(4,:) = Emg.tfl(:)';         % tensor_fascia_late
EMG(5,:) = Emg.semimem(:)';     % semimembranosus
EMG(6,:) = Emg.bifem(:)';       % biceps_femoris_long_head
EMG(7,:) = Emg.rf(:)';          % rectus_femoris
EMG(8,:) = Emg.vasmed(:)';      % vastus_medialis
EMG(9,:) = Emg.vaslat(:)';      % vastus_lateralis
EMG(10,:) = Emg.medgas(:)';     % gastrocnemius_medialis
EMG(11,:) = Emg.latgas(:)';     % gastrocnemius_lateralis
EMG(12,:) = Emg.soleus(:)';     % soleus
EMG(13,:) = Emg.tibant(:)';     % tibialis_anterior
EMG(14,:) = Emg.peronl(:)';     % peroneus_longus

MAX_EMG(1) = max(Emg.gmax);     % gluteus_maximus
MAX_EMG(2) = max(Emg.gmed);     % gluteus_medius
MAX_EMG(3) = max(Emg.addmagnus);% adductus_magnus
MAX_EMG(4) = max(Emg.tfl);      % tensor_fascia_late
MAX_EMG(5) = max(Emg.semimem);  % semimembranosus
MAX_EMG(6) = max(Emg.bifem);    % biceps_femoris_long_head
MAX_EMG(7) = max(Emg.rf);       % rectus_femoris
MAX_EMG(8) = max(Emg.vasmed);   % vastus_medialis
MAX_EMG(9) = max(Emg.vaslat);   % vastus_lateralis
MAX_EMG(10) = max(Emg.medgas);  % gastrocnemius_medialis
MAX_EMG(11) = max(Emg.latgas);  % gastrocnemius_lateralis
MAX_EMG(12) = max(Emg.soleus);  % soleus
MAX_EMG(13) = max(Emg.tibant);  % tibialis_anterior
MAX_EMG(14) = max(Emg.peronl);  % peroneus_longus

MIN_EMG(1) = min(Emg.gmax);     % gluteus_maximus
MIN_EMG(2) = min(Emg.gmed);     % gluteus_medius
MIN_EMG(3) = min(Emg.addmagnus);% adductus_magnus
MIN_EMG(4) = min(Emg.tfl);      % tensor_fascia_late
MIN_EMG(5) = min(Emg.semimem);  % semimembranosus
MIN_EMG(6) = min(Emg.bifem);    % biceps_femoris_long_head
MIN_EMG(7) = min(Emg.rf);       % rectus_femoris
MIN_EMG(8) = min(Emg.vasmed);   % vastus_medialis
MIN_EMG(9) = min(Emg.vaslat);   % vastus_lateralis
MIN_EMG(10) = min(Emg.medgas);  % gastrocnemius_medialis
MIN_EMG(11) = min(Emg.latgas);  % gastrocnemius_lateralis
MIN_EMG(12) = min(Emg.soleus);  % soleus
MIN_EMG(13) = min(Emg.tibant);  % tibialis_anterior
MIN_EMG(14) = min(Emg.peronl);  % peroneus_longus

% -------------------------------------------------------------------------
% MFP  Data
% -------------------------------------------------------------------------
MFP = zeros(14,n);
MAX_MFP = zeros(14,1);
Model.X = squeeze(Model.X);

MFP(1,:) = Model.X(1,:)+Model.X(2,:)+Model.X(3,:);      % gluteus_maximus
MFP(2,:) = Model.X(4,:)+Model.X(5,:)+Model.X(6,:);      % gluteus_medius
MFP(3,:) = Model.X(12,:)+Model.X(13,:)+Model.X(14,:);   % adductus_magnus
MFP(4,:) = Model.X(21,:);                               % tensor_fascia_late
MFP(5,:) = Model.X(24,:);                               % semimembranosus
MFP(6,:) = Model.X(26,:);                               % biceps_femoris_long_head
MFP(7,:) = Model.X(28,:);                               % rectus_femoris
MFP(8,:) = Model.X(29,:);                               % vastus_medialis
MFP(9,:) = Model.X(31,:);                               % vastus_lateralis
MFP(10,:) = Model.X(32,:);                              % gastrocnemius_medialis
MFP(11,:) = Model.X(33,:);                              % gastrocnemius_lateralis
MFP(12,:) = Model.X(34,:);                              % soleus
MFP(13,:) = Model.X(36,:);                              % tibialis_anterior
MFP(14,:) = Model.X(39,:);                              % peroneus_longus

MAX_MFP(1) = max(Model.X(1,:))+max(Model.X(2,:))+max(Model.X(3,:));   % gluteus_maximus
MAX_MFP(2) = max(Model.X(4,:))+max(Model.X(5,:))+max(Model.X(6,:));   % gluteus_medius
MAX_MFP(3) = max(Model.X(12,:))+max(Model.X(13,:))+max(Model.X(14,:));% adductus_magnus
MAX_MFP(4) = max(Model.X(21,:));                                      % tensor_fascia_late
MAX_MFP(5) = max(Model.X(24,:));                                      % semimembranosus
MAX_MFP(6) = max(Model.X(26,:));                                      % biceps_femoris_long_head
MAX_MFP(7) = max(Model.X(28,:));                                      % rectus_femoris
MAX_MFP(8) = max(Model.X(29,:));                                      % vastus_medialis
MAX_MFP(9) = max(Model.X(31,:));                                      % vastus_lateralis
MAX_MFP(10) = max(Model.X(32,:));                                     % gastrocnemius_medialis
MAX_MFP(11) = max(Model.X(33,:));                                     % gastrocnemius_lateralis
MAX_MFP(12) = max(Model.X(34,:));                                     % soleus
MAX_MFP(13) = max(Model.X(36,:));                                     % tibialis_anterior
MAX_MFP(14) = max(Model.X(39,:));                                     % peroneus_longus
   
% -------------------------------------------------------------------------
% PERRY's data
% -------------------------------------------------------------------------
onoff_Perry = zeros(14,7);
onoff_Perry(1,:) = [1 0 0 0 0 0 1] ;    % gluteus_maximus
onoff_Perry(2,:) = [1 1 0 0 0 0 0] ;    % gluteus_medius
onoff_Perry(3,:) = [1 0 0 0 0 0 1] ;    % adductus_magnus
onoff_Perry(4,:) = [0 0 1 0 0 0 0] ;    % tensor_fascia_late
onoff_Perry(5,:) = [1 0 0 0 0 0 1] ;    % semimembranosus
onoff_Perry(6,:) = [1 0 0 0 0 0 1] ;    % biceps_femoris_long_head
onoff_Perry(7,:) = [0 0 0 1 0 0 0] ;    % rectus_femoris
onoff_Perry(8,:) = [1 0 0 0 0 0 1] ;    % vastus_medialis
onoff_Perry(9,:) = [1 0 0 0 0 0 1] ;    % vastus_lateralis
onoff_Perry(10,:) = [0 1 1 0 0 0 0] ;   % gastrocnemius_medialis
onoff_Perry(11,:) = [0 1 1 0 0 0 0] ;   % gastrocnemius_lateralis
onoff_Perry(12,:) = [0 1 1 0 0 0 0] ;   % soleus
onoff_Perry(13,:) = [1 0 0 1 1 1 1] ;   % tibialis_anterior
onoff_Perry(14,:) = [0 1 1 0 0 0 0] ;   % peroneus_longus

% -------------------------------------------------------------------------
% Construction of on/off Matrix 
% -------------------------------------------------------------------------
%initialization
onoff_EMG = zeros(14,7);
onoff_MFP = zeros(14,7);
Ratio_EMG = zeros(14,1);
Ratio_Perry = zeros(14,1);

thresholdEMG = 0.20;
thresholdMFP = 0.10;
for m = 1:14
    % EMG matrix
    if mean(EMG(m,Phase1)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
       onoff_EMG(m,1) = 1;
    end
    if mean(EMG(m,Phase2)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
       onoff_EMG(m,2) = 1;
    end
    if mean(EMG(m,Phase3)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
       onoff_EMG(m,3) = 1;
    end
    if mean(EMG(m,Phase4)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
       onoff_EMG(m,4) = 1;
    end
    if mean(EMG(m,Phase5)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
        onoff_EMG(m,5) = 1;
    end
    if mean(EMG(m,Phase6)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
          onoff_EMG(m,6) = 1;
    end
    if mean(EMG(m,Phase7)) > thresholdEMG*(MAX_EMG(m)-MIN_EMG(m))+MIN_EMG(m);
       onoff_EMG(m,7) = 1;
    end
    % MFP Matrix
    if mean(MFP(m,Phase1)) > thresholdMFP*MAX_MFP(m,:);
       onoff_MFP(m,1) = 1;
    end
    if mean(MFP(m,Phase2)) > thresholdMFP*MAX_MFP(m);
       onoff_MFP(m,2) = 1;
    end
    if mean(MFP(m,Phase3)) > thresholdMFP*MAX_MFP(m);
       onoff_MFP(m,3) = 1;
    end
    if mean(MFP(m,Phase4)) > thresholdMFP*MAX_MFP(m);
       onoff_MFP(m,4) = 1;
    end
    if mean(MFP(m,Phase5)) > thresholdMFP*MAX_MFP(m);
        onoff_MFP(m,5) = 1;
    end
    if mean(MFP(m,Phase6)) > thresholdMFP*MAX_MFP(m);
          onoff_MFP(m,6) = 1;
    end
    if mean(MFP(m,Phase7)) > thresholdMFP*MAX_MFP(m);
       onoff_MFP(m,7) = 1;
    end
end

% Concordance Matrix
onoff_EMGMFP = onoff_MFP-onoff_EMG;     % 0 = concordance ; 1 = EMG:off MFP:on; -1 = EMF:on MFP:off
onoff_PerryMFP = onoff_MFP-onoff_Perry; % 0 = concordance ; 1 = EMG:off MFP:on; -1 = EMF:on MFP:off

% Percentil of validation
length(find(onoff_EMGMFP==0))*100/(7*14)
Model.concordance_EMG = length(find(onoff_EMGMFP==0))*100/(7*14)     % '%' concordance EMG/model
Model.concordance_Perry = length(find(onoff_PerryMFP==0))*100/(7*14) % '%' concordance Perry/model

% Dickerson's ratio
for m = 1:14
    Ratio_EMG(m) = length(find(onoff_EMGMFP(m,:)==0))/(length(find(onoff_EMGMFP(m,:)==1))+length(find(onoff_EMGMFP(m,:)==-1)));
    Ratio_Perry(m) = length(find(onoff_PerryMFP(m,:)==0))/(length(find(onoff_PerryMFP(m,:)==1))+length(find(onoff_PerryMFP(m,:)==-1)));
end
