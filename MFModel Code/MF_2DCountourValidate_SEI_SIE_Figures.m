load('ContourData_S_EE=0.029-40s.mat','ContourData');
load('ContourData_S_EE=0.029NWValidate.mat','ContourValidate');

Fr_MFComp = zeros(size(ContourValidate.Fr_NW_Valid));

for FrInd = 1:size(Fr_MFComp,2)
    Fr_MFComp(:,FrInd) = ContourData.Fr_NoFix(:,ContourValidate.NWtestSeq(FrInd,1),ContourValidate.NWtestSeq(FrInd,2));
    
end

Accuracy = abs(Fr_MFComp - ContourValidate.Fr_NW_Valid)./ContourValidate.Fr_NW_Valid;
c = max(Accuracy);
c(isnan(c)) = 1;
c(isinf(c)) = 1;
hold on
scatter(ContourData.S_IEtest(ContourValidate.NWtestSeq(:,2)),ContourData.S_EItest(ContourValidate.NWtestSeq(:,1)),25,c,'filled')
hold off
colorbar
caxis([0, 0.4]);
xlabel('S_{IE}');ylabel('S_{EI}')
title('S_{EE} = 0.029')