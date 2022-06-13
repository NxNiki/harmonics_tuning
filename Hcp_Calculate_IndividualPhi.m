function TransformedPhi=Hcp_Calculate_IndividualPhi(IndividualPhi,IndividualL,Phi,lambda)
NodeNum=size(IndividualL,1);
p=size(IndividualPhi,2);
Phi_temp=IndividualPhi;
L=IndividualL;
B=lambda*Phi;
E=eig(L);
alpha=E(NodeNum);
L_tilde=alpha*eye(NodeNum,NodeNum)-L;
iter1=1;
Step1Diff=1;

% Step1ObjectiveFuncValue=zeros(2,1);
% changed by Xin:
maxiter = 100;
Step1ObjectiveFuncValue=zeros(maxiter,1); 

Step1ObjectiveFuncValue(iter1)=trace(Phi_temp'*L*Phi_temp)+lambda*trace((Phi_temp-Phi)'*(Phi_temp-Phi));
while Step1Diff>0.00001&&iter1<maxiter          
    M=L_tilde*Phi_temp+B;
    
    % replace NaNs with 0. by Xin Mar-31-2022.
    % M(isnan(M)) = 0;
    
    [U,S,V]=svd(M);
    Phi_temp=U(:,1:p)*V';
    iter1=iter1+1;
    Step1ObjectiveFuncValue(iter1)=trace(Phi_temp'*L*Phi_temp)+lambda*trace((Phi_temp-Phi)'*(Phi_temp-Phi));
    Step1Diff=abs(Step1ObjectiveFuncValue(iter1)-Step1ObjectiveFuncValue(iter1-1));           
end
TransformedPhi=Phi_temp;
end