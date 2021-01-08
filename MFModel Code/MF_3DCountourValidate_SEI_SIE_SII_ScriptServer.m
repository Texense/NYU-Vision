%% A Rough Estimation Contour for S_EI and S_IE
% first, setup connctivity map
CurrentFolder = pwd
FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])

N_HC = 3; 
% Number of E and I neurons
n_E_HC = 54; n_I_HC = 31; % per side of HC
N_E = n_E_HC^2 * N_HC^2; % neuron numbers In all
N_I = n_I_HC^2 * N_HC^2;
% Grid sizes of E and I neurons; 
Size_HC = 0.500; % in mm;
Size_E = Size_HC/n_E_HC; Size_I = Size_HC/n_I_HC; 
% Projection: SD of distances
SD_E = 0.2/sqrt(2); SD_I = 0.125/sqrt(2);
Dist_LB = 0.36; % ignore the connection probability of dist>0.3mm
% Peak probability of projection
Peak_EE = 0.15; Peak_I = 0.6; 

% spatial indexes of E and I neurons
[NnE.X,NnE.Y] = V1Field_Generation(N_HC,1:N_E,'e');
[NnI.X,NnI.Y] = V1Field_Generation(N_HC,1:N_I,'i');

% determine connections between E, I
% sparse metrices containing 0 or 1. 
% Row_i Column_j means neuron j projects to neuron i
% add periodic boundary
C_EE = ConnectionMat(N_E,NnE,Size_E,...
                     N_E,NnE,Size_E,...
                     Peak_EE,SD_E,Dist_LB,1);

C_EI = ConnectionMat(N_E,NnE,Size_E,...
                     N_I,NnI,Size_I,...
                     Peak_I,SD_I,Dist_LB,0);

C_IE = ConnectionMat(N_I,NnI,Size_I,...
                     N_E,NnE,Size_E,...
                     Peak_I,SD_E,Dist_LB,0);

C_II = ConnectionMat(N_I,NnI,Size_I,...
                     N_I,NnI,Size_I,...
                     Peak_I,SD_I,Dist_LB,1);

%% Variables and Parameters
%RefTimeE = zeros(N_E,1); VE = 0.5*rand(N_E,1)-0.5; SpE = sparse(N_E,1); GE_ampa_R = zeros(N_E,1); GE_nmda_R = zeros(N_E,1); GE_gaba_R = zeros(N_E,1); GE_ampa_D = zeros(N_E,1); GE_nmda_D = zeros(N_E,1); GE_gaba_D = zeros(N_E,1);
%RefTimeI = zeros(N_I,1); VI = 1.5*rand(N_I,1)-0.5; SpI = sparse(N_I,1); GI_ampa_R = zeros(N_I,1); GI_nmda_R = zeros(N_I,1); GI_gaba_R = zeros(N_I,1); GI_ampa_D = zeros(N_I,1); GI_nmda_D = zeros(N_I,1); GI_gaba_D = zeros(N_I,1); 
InSs = load('Initials.mat');
%parameters
S_EE = 0.024; S_EI = 0.055*0.024/0.029; %S_IE = 0.0081; S_II = 0.048; %Tuned
%S_EE = 0.033; S_EI = 0.061; S_IE = 0.0087; S_II = 0.048;
%S_EE = 0.033; S_EI = 0.061; S_IE = 0.0087; S_II = 0.048;
%S_EE = 0.029; S_EI = 0.055; S_IE = 0.0081; S_II = 0.048; %Tuned
%S_EE = 0.028; S_EI = 0.056; S_IE = 0.0095; S_II = 0.042; % original
p_EEFail = 0.2; S_amb = 0.01;

tau_ampa_R = 0.5; tau_ampa_D = 3;
tau_nmda_R = 2; tau_nmda_D = 80;
tau_gaba_R = 0.5; tau_gaba_D = 5;
tau_ref = 2; % time unit is ms
dt = 0.2;
T = 3000; sampleT = 2; SimulationT = 2000;
gL_E = 1/20;  Ve = 14/3; S_Elgn = 0.059; rhoE_ampa = 0.8; rhoE_nmda = 0.2;
gL_I = 1/15;  Vi = -2/3; S_Ilgn = 0.084; rhoI_ampa = 0.67;rhoI_nmda = 0.33;

lambda_E = 0.08; % ~16 LGN spike can excite a E neurons. 0.25 spike/ms makes 64 ms for such period. 
lambda_I = 0.08; 
%rE_amb = 0.72; rI_amb = 0.36;
rE_amb = 0.54; rI_amb = 0.54;
% Replace S_EI by testing values
GridNum1 = 32;
GridNum2 = 32;
GridNum3 = 15;
S_EI_Mtp = [0.4, 4]; % of S_EE
S_IE_Mtp = [0.05,0.5]; % of S_II
S_II_Mtp = [0.5 3.5]; % of current S_EI
S_EItest = linspace(S_EI_Mtp(1),S_EI_Mtp(2),GridNum1)*S_EE;
S_IEtestMtp = linspace(S_IE_Mtp(1),S_IE_Mtp(2),GridNum2);%*S_EE; I only specify a vecter length here
S_IItest = linspace(S_II_Mtp(1),S_II_Mtp(2),GridNum3)*S_EI;
% Add lines boundaries
LineL1 = polyfit([0.06 0.14],[2.5 0.4],1); % S_IEMtp first, second S_EIMtp. Those numbers are multipliers of S_II and S_EE
LineL2 = polyfit([0.06 0.28],[1.5 0.4],1);
LineU1 = polyfit([0.1  0.5 ],[4   1.5],1);

% defining cluster. 
cluster = gcp('nocreate');
if isempty(cluster)
    cluster = parpool([4 64]);
    cluster.IdleTimeout = 1200;
end

%% Network Simulation to check if meeting MF+LIF
NWtestNum = 64;
Fr_NW_Valid = zeros(2,NWtestNum,length(S_IItest));
mV_NW_Valid = zeros(2,NWtestNum,length(S_IItest));

NWtestSeq3D = zeros(NWtestNum,2,length(S_IItest)); % 

E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
[E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);
    
I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
[I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);

N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));


for S_IIInd = 1:length(S_IItest)
    S_II = S_IItest(S_IIInd);
       

% Randomly select testing points
NWtestSeq = [];
Indtest = 0;
while (Indtest<NWtestNum)
    S_EIind = randi(GridNum1);   S_IEind = randi(GridNum2);
    S_EI_NW = S_EItest(S_EIind); S_IE_NW = S_IEtestMtp(S_IEind)*S_II;
    if (S_EI_NW/S_EE<=S_IE_NW/S_II*LineL1(1)+LineL1(2) || S_EI_NW/S_EE<=S_IE_NW/S_II*LineL2(1)+LineL2(2) )
        continue
    end
    
    if (S_EI_NW/S_EE>=S_IE_NW/S_II*LineU1(1)+LineU1(2) )
        continue
    end
    
    NWtestSeq = [NWtestSeq;[S_EIind,S_IEind]];
    Indtest = Indtest+1;
end

Fr_NW_Valid_FixSII = zeros(2,NWtestNum);
mV_NW_Valid_FixSII = zeros(2,NWtestNum);

parfor ValiInd = 1:length(S_EItest)
    % In case NW explode for too long:     
    ExplodeCount = 0;
    JumpCollectData = false;
    
    tic
    RefTimeE = InSs.RefTimeE; VE = InSs.VE; SpE = InSs.SpE; GE_ampa_R = InSs.GE_ampa_R; GE_nmda_R = InSs.GE_nmda_R; GE_gaba_R = InSs.GE_gaba_R; 
                                                            GE_ampa_D = InSs.GE_ampa_D; GE_nmda_D = InSs.GE_nmda_D; GE_gaba_D = InSs.GE_gaba_D;
    RefTimeI = InSs.RefTimeI; VI = InSs.VI; SpI = InSs.SpI; GI_ampa_R = InSs.GI_ampa_R; GI_nmda_R = InSs.GI_nmda_R; GI_gaba_R = InSs.GI_gaba_R; 
                                                            GI_ampa_D = InSs.GI_ampa_D; GI_nmda_D = InSs.GI_nmda_D; GI_gaba_D = InSs.GI_gaba_D; 
  
    S_EI = S_EItest(NWtestSeq(ValiInd,1));
    S_IE = S_IEtestMtp(NWtestSeq(ValiInd,2))*S_II;
    %E_SpHis= []; I_SpHis=[];
    E_Sp = []; I_Sp = [];
    VE_T = []; VI_T = [];
    sampleN = floor(sampleT/dt); % sample each 2 ms
    SimulationN = floor(SimulationT/dt); % show and check every 200ms
    for TimeN = 1:floor(T/dt)
    [oRefTimeE,oVE,oSpE,oGE_ampa_R,oGE_nmda_R,oGE_gaba_R,...
          oGE_ampa_D,oGE_nmda_D,oGE_gaba_D,...
     oRefTimeI,oVI,oSpI,oGI_ampa_R,oGI_nmda_R,oGI_gaba_R,...
          oGI_ampa_D,oGI_nmda_D,oGI_gaba_D] = ...
          V1NetworkUpdate(RefTimeE,VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R,...
                                          GE_ampa_D,GE_nmda_D,GE_gaba_D,...
                          RefTimeI,VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R,...
                                          GI_ampa_D,GI_nmda_D,GI_gaba_D,...
                          C_EE,C_EI,C_IE,C_II,...
                          S_EE,S_EI,S_IE,S_II,...
                          tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... % time unit is ms
                          dt,p_EEFail,...
                          gL_E,Ve,S_Elgn,rhoE_ampa,rhoE_nmda,...
                          gL_I,Vi,S_Ilgn,rhoI_ampa,rhoI_nmda,...
                          S_amb,lambda_E,lambda_I,rE_amb,rI_amb);
                      
    RefTimeE = oRefTimeE; VE = oVE;SpE = oSpE;GE_ampa_R = oGE_ampa_R; GE_nmda_R = oGE_nmda_R; GE_gaba_R = oGE_gaba_R;
                                              GE_ampa_D = oGE_ampa_D; GE_nmda_D = oGE_nmda_D; GE_gaba_D = oGE_gaba_D;
    RefTimeI = oRefTimeI; VI = oVI;SpI = oSpI;GI_ampa_R = oGI_ampa_R; GI_nmda_R = oGI_nmda_R; GI_gaba_R = oGI_gaba_R;
                                              GI_ampa_D = oGI_ampa_D; GI_nmda_D = oGI_nmda_D; GI_gaba_D = oGI_gaba_D;


    E_Sp = [E_Sp;[find(oSpE),ones(size(find(oSpE)))*TimeN*dt]];
    I_Sp = [I_Sp;[find(oSpI),ones(size(find(oSpI)))*TimeN*dt]];
    VE_T = [VE_T;nanmean(oVE(E_Ind))];
    VI_T = [VI_T;nanmean(oVI(I_Ind))];
    
    % Count explode step. When exceed certain number, break
    if sum(isnan(oVE))>0.7*N_E
        ExplodeCount = ExplodeCount+1;  
    else
        ExplodeCount = 0;
    end
    
    if ExplodeCount>floor(T/dt)/10
        JumpCollectData = true;
        disp('NW Exploded too long, break...')
        break
    end
    
    end
    
    ParameterDisp = ['S_EE = ' num2str(S_EE) ', S_EI = ' num2str(S_EI) ...
                         ', S_IE = ' num2str(S_IE) ', S_II = ' num2str(S_II) '\n' ...
                         'S_Elgn = ' num2str(S_Elgn) ', lambda_E = ' num2str(lambda_E) ...
                         ', S_Ilgn = ' num2str(S_Ilgn) ', lambda_I = ' num2str(lambda_I) '\n' ...
                         'S_amb = ' num2str(S_amb) ', rE_amb = ' num2str(rE_amb) ', rI_amb = ' num2str(rI_amb)];
    
    fprintf(ParameterDisp)
     
    % if network exploded too long, give up collecting statistics
    if JumpCollectData 
    Fr_NW_Valid_FixSII(:,ValiInd) = [nan;nan];
    mV_NW_Valid_FixSII(:,ValiInd) = [nan,nan];
    else
    %E_Ind = 10000:13000; I_Ind = 3300:4300;
    scatterE = find(ismember(E_Sp(:,1),E_Ind));
    scatterI = find(ismember(I_Sp(:,1),I_Ind));
    [~,E_Fire_Ind] = ismember(E_Sp(scatterE,1),E_Ind);
    [~,I_Fire_Ind] = ismember(I_Sp(scatterI,1),I_Ind);

    WinSize = SimulationT;
    T_RateWindow = [TimeN*dt-WinSize TimeN*dt];
    E_SpInd = find(E_Sp(:,2)>=T_RateWindow(1) & E_Sp(:,2)<=T_RateWindow(2) & ismember(E_Sp(:,1),E_Ind));
    I_SpInd = find(I_Sp(:,2)>=T_RateWindow(1) & I_Sp(:,2)<=T_RateWindow(2) & ismember(I_Sp(:,1),I_Ind));
    E_Rate = length(E_SpInd)/(WinSize/1000)/length(E_Ind);
    I_Rate = length(I_SpInd)/(WinSize/1000)/length(I_Ind);
    
    Fr_NW_Valid_FixSII(:,ValiInd) = [E_Rate;I_Rate];
    mV_NW_Valid_FixSII(:,ValiInd) = [nanmean(VE_T),nanmean(VI_T)];
    end
    toc
end
Fr_NW_Valid(:,:,S_IIInd) = Fr_NW_Valid_FixSII;
mV_NW_Valid(:,:,S_IIInd) = mV_NW_Valid_FixSII;
NWtestSeq3D(:,:,S_IIInd) = NWtestSeq;
end

%% save data
CommentString = sprintf('NWValidate_S_II_rEamb%.2f_rIamb%.2f',rE_amb,rI_amb);
ContourValidate = struct('Fr_NW_Valid',Fr_NW_Valid,'mV_NW_Valid',mV_NW_Valid,'NWtestSeq3D',NWtestSeq3D,...
                         'S_EE',S_EE, 'S_EItest',S_EItest, 'S_IItest',S_IItest, 'S_IEtestMtp',S_IEtestMtp,...
                         'rE_amb',rE_amb, 'rI_amb',rI_amb);
save(['ContourData_S_EE=' num2str(S_EE) CommentString '.mat'],'ContourValidate')