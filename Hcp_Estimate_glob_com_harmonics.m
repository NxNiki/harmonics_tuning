function [CommonNetwork,CommonHarmonics]=Hcp_Estimate_glob_com_harmonics(Graph,p)

%% Initialize common Harmonics as the average of all subject's networks
SubjectNum=size(Graph,2);
NodeNum=size(Graph(1).W,1);

% calculate CommonNetwork by the average of individual networks
CommonNetwork=zeros(NodeNum);
for i=1:SubjectNum
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/SubjectNum;
% CommonNetwork(CommonNetwork<0.005)=0;

% CommonNetwork should already be a symmetric matrix, what is the purpose of this step?
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

% calculate graph Laplacian:
temp_D=diag(sum(CommonNetwork,2));
LatentLaplacian=temp_D-CommonNetwork;
[Phi_temp,~]=eig(LatentLaplacian);
Phi_ave=Phi_temp(:,1:p);

%% Estimating global common harmonic waves
% Initializing parameters
% lambda_1=0.001;
lambda_1 = .01;

gama=0.001/(lambda_1);   
iter=1;
maxiter = 50; % add by Xin

% Phi=cell(5,1);
Phi=cell(maxiter,1);

Diff=zeros(maxiter,1);
Diff(1)=100;

% ObjectiveFuncValue=zeros(2,1);
ObjectiveFuncValue=zeros(maxiter,1);

Phi{1}=Phi_ave;

% Runing Algrithom   
%while Diff(iter)>0.1&&iter<50
while Diff(iter)>0.1&&iter<maxiter
    %Step1: Updating individual Phi
    for i=1:SubjectNum
        Graph(i).Phi{iter+1}=Calculate_IndividualPhi(Graph(i).Phi{iter},Graph(i).L,Phi{iter},lambda_1);  
    end

    %Step2: Updating Common Harmonics Phi    
    Phi_k=Graph(20).Phi{iter+1}; %% ** WHY USING THE 20TH GRAPH? XIN.**
    err=1;    
    iter2=1;
    
    %Step2ObjectFunction=zeros(2,1);
    Step2ObjectFunction=zeros(500,1); %% what is the purpose of the Step2ObjectFunction? Xin.
    for i=1:SubjectNum
        Step2ObjectFunction(iter2)=Step2ObjectFunction(iter2)+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi_k));
    end
    while err>0.0001&&iter2<500
        Phi_increment=zeros(NodeNum,p);
        for i=1:SubjectNum
            Phi_increment=Phi_increment+Phi_k*Graph(i).Phi{iter+1}'*Phi_k-Graph(i).Phi{iter+1};
        end
        Phi_increment=-gama*lambda_1*Phi_increment;
        [Q,R]=qr((eye(NodeNum)-Phi_k*Phi_k')*Phi_increment,0);
        A=Phi_k'*Phi_increment;
        BC=expm([A,-R';R,zeros(p)])*[eye(p);zeros(p)];        
        Phi{iter+1}=Phi_k*BC(1:p,:)+Q*BC(p+1:2*p,:);
        Phi_k=Phi{iter+1};        
        err=norm(Phi_increment,'fro');
        temp=0;
        for i=1:SubjectNum
            % temp=temp+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi_k));
            
            % check NaN value added by Xin Mar-31-2022.
            ttemp = lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi_k));
            if isnan(ttemp)
                sprintf('NaN in Step2ObjectFunction found for subject %d', i);
                continue
            end
            temp=temp+ttemp;
        end
        Step2ObjectFunction(iter2+1)=temp;
        iter2=iter2+1;
    end
    % Convergent condition
    temp=0;
    for i=1:SubjectNum
        %temp=temp+trace(Graph(i).Phi{iter+1}'*Graph(i).L*Graph(i).Phi{iter+1})+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi{iter+1}));
        
        % check NaN. added by Xin Mar-31-2022.
        ttemp = trace(Graph(i).Phi{iter+1}'*Graph(i).L*Graph(i).Phi{iter+1})+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi{iter+1}));
        if isnan(ttemp)
            sprintf('NaN in ObjectiveFuncValue found for subject %d', i);
            continue
        end
        temp = temp + ttemp;
    end
    ObjectiveFuncValue(iter+1)=temp;
    Diff(iter+1)=abs(ObjectiveFuncValue(iter+1)- ObjectiveFuncValue(iter));   
    fprintf('The iteration: %d, the error: %f ....\n', iter, Diff(iter+1));
   
    
    iter=iter+1;   
end
CommonHarmonics=Phi{iter};
end
