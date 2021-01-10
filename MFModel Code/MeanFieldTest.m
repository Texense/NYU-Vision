%% Test BkGd Input
S_EE = 0.029; S_EI = 0.055; S_IE = 0.0081; S_II = 0.048;

rE_amb = 0.72; rI_amb = 0.32; S_amb = 0.01;

 lambda_Eall = [0.01:0.01:0.3];
 f_EnI = [zeros(size(lambda_Eall));zeros(size(lambda_Eall))];
 meanVs = zeros(size(f_EnI));
 loop = zeros(size(lambda_Eall));
 for lamInd = 1:length(lambda_Eall)
tic
lambda_E = lambda_Eall(lamInd);
lambda_I = 0.8*lambda_E;
[f_EnI(:,lamInd),meanVs(:,lamInd),loop(lamInd)] = MeanFieldEst_BkGd_Indep(C_EE,C_EI,C_IE,C_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                   rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                   gL_E,gL_I,Ve,Vi);
toc
 end
 
 figure(5)
 hold on
 subplot 211
 plot(lambda_Eall(1:end),f_EnI(1,1:end))
 subplot 212
 plot(lambda_Eall(1:end),f_EnI(2,1:end))