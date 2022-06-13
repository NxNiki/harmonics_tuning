close all; clear;clc
load('result/CommonGlobalHarmonics.mat')
load('result/CommonGlobalHarmonics_65.mat')
load('result/CommonLocalizedRegion.mat')

%% 
LocalizedRegion = CommonLocalizedRegion;
CommonHarmonics_65 = CommonGlobalHarmonics_65;
CommonHarmonics = CommonGlobalHarmonics;

[AD, CN, EMCI, LMCI, SMC] = Proprocess_original_data('Data/signal Data/Amyloid_SUVR.xlsx');
Combined_Amyloid_all = [AD.Data, CN.Data, EMCI.Data, SMC.Data, LMCI.Data];
[n, m] = size(Combined_Amyloid_all);

Total_C_Reconstruction_loss = zeros(m,1);
Total_C65_Reconstruction_loss = zeros(m,1);
Total_L_Reconstruction_loss = zeros(m,1);

for j=1:n

    L_nonzero_index = find(LocalizedRegion(j).Region_mask > 0);
    L_Combined_Amyloid_all = zeros(size(Combined_Amyloid_all));
    L_Combined_Amyloid_all(L_nonzero_index,:) = Combined_Amyloid_all(L_nonzero_index,:);
    
    % Basic (55) common harmonics reconstruction
    CommonHarmonics_Power = CommonHarmonics' * L_Combined_Amyloid_all;
    CommonHarmonics_Reconstruction = CommonHarmonics * CommonHarmonics_Power;
    C_Reconstruction_loss = zeros(m,1);
    
    for i=1:m
        C_Reconstruction_loss(i) = abs(CommonHarmonics_Reconstruction(j,i) - L_Combined_Amyloid_all(j,i));
    end
    Total_C_Reconstruction_loss = Total_C_Reconstruction_loss + C_Reconstruction_loss;
    
    % 70 common harmonics reconstruction
    CommonHarmonics_65_Power = CommonHarmonics_65' * L_Combined_Amyloid_all;
    CommonHarmonics_65_Reconstruction = CommonHarmonics_65 * CommonHarmonics_65_Power;
    C65_Reconstruction_loss = zeros(m, 1);
    
    for i=1:m
        C65_Reconstruction_loss(i) = abs(CommonHarmonics_65_Reconstruction(j,i) - L_Combined_Amyloid_all(j,i));
    end
    Total_C65_Reconstruction_loss = Total_C65_Reconstruction_loss + C65_Reconstruction_loss;

    % Localized harmonics reconstruction
    LocalizedHarmonics_Power = LocalizedRegion(j).Harmonics' * L_Combined_Amyloid_all;
    LocalizedHarmonics_Reconstruction = LocalizedRegion(j).Harmonics * LocalizedHarmonics_Power;
    CommonHarmonics_Power = CommonHarmonics' * L_Combined_Amyloid_all;
    CommonHarmonics_Reconstruction = CommonHarmonics * CommonHarmonics_Power;    
    All_Reconstruction = CommonHarmonics_Reconstruction + LocalizedHarmonics_Reconstruction;
    L_Reconstruction_loss = zeros(m,1);
    
    for i=1:m
       L_Reconstruction_loss(i) = abs(All_Reconstruction(j,i) - L_Combined_Amyloid_all(j,i));
    end
    
    Total_L_Reconstruction_loss = Total_L_Reconstruction_loss + L_Reconstruction_loss;
    
end

Total_C_Reconstruction_loss = Total_C_Reconstruction_loss/n;
Total_C65_Reconstruction_loss = Total_C65_Reconstruction_loss/n;
Total_L_Reconstruction_loss = Total_L_Reconstruction_loss/n;
mean(Total_C_Reconstruction_loss)
mean(Total_C65_Reconstruction_loss)
mean(Total_L_Reconstruction_loss)

figure
boxplot([Total_C_Reconstruction_loss, Total_C65_Reconstruction_loss, Total_L_Reconstruction_loss], 'Labels', {'55 CH','65 CH', '55+10 LCH'})
ylabel('Reconstruction loss')