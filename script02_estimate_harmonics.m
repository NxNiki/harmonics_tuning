%% load brain connectivity network
% this script is used to tune the parameters with signal reconstruction.
% Xin June 13 2022.

% NOTE: to prevent system crashing, do not process more than 100 files!

close all;
clear;clc


% add more files here:
input_dir = "/home/xin/Downloads/Harmonics_signal_reconstruction/hcp_harmonics/hcp_out01_corr_matrix_scan1/";
% work_dir = "/home/xin/Downloads/Harmonics_signal_reconstruction/wd_hcp/hcp_harmonic_wavelets";

% use training set to estimate harmonics wavelets:
files = dir(strcat(input_dir, '/corr_matrix*.csv'));
files = files(1:9:end);

output_dir = "hcp_out02_harmonics";
out_file_name = 'harmonics_all_65.mat';

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

%% construct graph:

Graph = struct;
GroupNum = 1;
% p=55;% p is the number of eignvectors %55
p = 65;

for i = 1:size(files, 1)
    
    f = strcat(files(i).folder, filesep, files(i).name);
    
    fprintf(f);
    fprintf('\n');
    
    temp = readmatrix(f);
    temp = temp(2:end,:); % remove header in correlation matrix.
    
    subject_id = f(strfind(f, 'sub-'): strfind(f, '.txt')-1);
    
    % check for NaNs:
    if any(isnan(temp), 'all')
        fprintf('NaN values on matrix for subject: %d', subject_id)
        break
    end
    
    if ~isempty(temp)
        temp(temp<0.007)=0; % all correlations less than .007 are removed.
        v=sum(temp,2);

        if sum(v==0)==0
            %disp('compute harmonics...')
            D=diag(v); 
            temp=D^-1*temp;
            Graph(GroupNum).W=(temp+temp')/2; %Adjacency matrix
            Graph(GroupNum).D=diag(sum(Graph(GroupNum).W,2)); % Degree matrix
            Graph(GroupNum).L=Graph(GroupNum).D-Graph(GroupNum).W;% Laplacian matrix
            [Phi_temp,value]=eig(Graph(GroupNum).L);

            if sum(Phi_temp(:,1))<0
                Phi_temp=-Phi_temp;
            end
            
            Graph(GroupNum).Phi{1}=Phi_temp(:,1:p);
            Graph(GroupNum).Eigenvalue=value(1:p,1:p);
            Graph(GroupNum).SubjectID=subject_id;
            %Graph(GroupNum).PTID=Data_profile{i,2};
            %Graph(GroupNum).VISCODE=Data_profile{i,4};
            %Graph(GroupNum).DX_bl=Data_profile{i,5};
            %Graph(GroupNum).Age=Data_profile{i,6};
            %Graph(GroupNum).Gender=Data_profile{i,7};

            GroupNum=GroupNum+1;
        end
    end
end

save(strcat(output_dir, '/', out_file_name), 'Graph')

%% Estimating global common harmonic waves
load(strcat(output_dir, '/', out_file_name), 'Graph')
fprintf('Estimating global common harmonic waves!\n');
CommonHarmonics = Hcp_Estimate_glob_com_harmonics(Graph, p);
save(strcat(output_dir, '/', out_file_name), 'CommonHarmonics', '-append')
fprintf('Done!\n')

%% Constructing region-adaptive individual harmonic wavelets
fprintf('Constructing region-adaptive individual harmonic wavelets!\n');
% load(strcat(output_dir, '/', out_file_name), 'CommonHarmonics', 'Graph')
Graph = Hcp_Construct_individual_harmonic_wavelets(CommonHarmonics, Graph);
save(strcat(output_dir, '/', out_file_name), 'Graph', '-append')
fprintf('Done!\n')

%% Identifying region-adaptive common harmonic wavelets
fprintf('Identifying region-adaptive common harmonic wavelets!\n');
% load(strcat(output_dir, '/', out_file_name), 'Graph')
CommonHarWavelets = Hcp_Identify_glob_com_har_wavelets(Graph, CommonHarmonics);
CommonHarmonics = CommonHarmonics(:, 1:p);
save(strcat(output_dir, '/', out_file_name), 'CommonHarWavelets', 'CommonHarmonics', '-append')
fprintf('Finished!\n')


%% save results to txt files:
% writematrix(CommonHarmonics, strcat(output_dir, '/CommonHarmonics.csv'))
% 
% for i = 1:length(CommonHarmonics)
%     
%     writematrix(CommonHarWavelets(i).Region_mask, strcat(output_dir, '/CommonHarWavelets_Region_Mask', sprintf('%03d', i), '.csv'))
%     writematrix(CommonHarWavelets(i).Harmonics, strcat(output_dir, '/CommonHarWavelets_Harmonics', sprintf('%03d',i), '.csv'))
%     
% end



