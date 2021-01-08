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
load('Initials.mat')
%parameters
S_EE = 0.024; S_II = 0.120; %S_IE = 0.0081; S_II = 0.048; %Tuned
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
gL_E = 1/20;  Ve = 14/3; S_Elgn = 2*S_EE; rhoE_ampa = 0.8; rhoE_nmda = 0.2;
gL_I = 1/15;  Vi = -2/3; %S_Ilgn = 0.084; 
rhoI_ampa = 0.67;rhoI_nmda = 0.33;

lambda_E = 0.08; % ~16 LGN spike can excite a E neurons. 0.25 spike/ms makes 64 ms for such period. 
lambda_I = 0.08; 
%rE_amb = 0.72; rI_amb = 0.36;
rE_amb = 0.50; rI_amb = 0.50;

%L6 input
S_IEOneTime = 0.20*S_II;
S_EL6 = 1/3*S_EE; S_IL6 = 1/3*S_IEOneTime;
rE_L6 = 0.25; % rI_L6 to be determined

% Replace S_EI by testing values
GridNum1 = 16;
GridNum2 = 16;
S_EI_Mtp = [0.4, 4]; % of S_EE
S_IE_Mtp = [0.05,0.5]; % of S_II
S_EItest = linspace(S_EI_Mtp(1),S_EI_Mtp(2),GridNum1)*S_EE;
S_IEtest = linspace(S_IE_Mtp(1),S_IE_Mtp(2),GridNum2)*S_II;%*S_EE; I only specify a vecter length here
% Panel: Two proportions
PanelNum1 = 5;
PanelNum2 = 6;
S_Ilgn_Mtp = [1, 3]; % of S_Elgn
rI_L6_Mtp  = [1, 6]; % of rE_L6
S_Ilgntest = linspace(S_Ilgn_Mtp(1),S_Ilgn_Mtp(2),PanelNum1)*S_Elgn;
rI_L6test = linspace(rI_L6_Mtp(1),rI_L6_Mtp(2),PanelNum2)*rE_L6;
%% Calculate input current
mVE = 0:0.1:1;
EInput = (Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb + rE_L6*S_EL6)*1000;
LeakE = -gL_E*mVE*1000;

figure('Name','net current input')
hold on
plot(mVE, EInput+LeakE,'r','DisplayName','Net input to E')

mVI = 0:0.1:1;
for S_Ilgn = S_Ilgntest
    for rI_L6 = rI_L6test

IInput = (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb + rI_L6*S_IL6)*1000;
LeakI = -gL_I*mVI*1000;
plot(mVI, IInput+LeakI, 'b','DisplayName',sprintf('Net input to I, S_Ilgn = %.3f; rI_{L6} = %d',S_Ilgn,floor(rI_L6*1000)) )
    end
end
plot(mVE, zeros(size(mVE)), 'DisplayName', 'zero')
legend
hold off
xlabel('mV'); ylabel('net current Input')

