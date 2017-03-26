function OmegaVal=LMI_Aut16a_th2(A,B,K,h,r1,muM,alpha,sigma)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "Predictor-based networked control under uncertain transmission delays," Automatica, vol. 70, pp. 101–108, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B   - the parameters of the system (3); 
% K     - the controller gain from (23); 
% h     - the maximum sampling period; 
% r1    - the controller-to-actuators known constant delay; 
% muM   - a bound for the controller-to-actuators time-varying delay uncertainty; 
% alpha - a desired decay rate; 
% sigma - the event-triggering parameter from (24); 

% Output: 
% OmegaVal - the value of Omega from (24). If OmegaVal is empty, the LMIs are not feasible

n=size(A,1); 
m=size(B,2); 

%% Decision variables 
P=sdpvar(n); 
S=sdpvar(n); 
S0=sdpvar(n); 
S1=sdpvar(n); 
R0=sdpvar(n); 
R1=sdpvar(n); 
Omega=sdpvar(m); 
P2=sdpvar(n,n,'f'); 
P3=sdpvar(n,n,'f'); 
G0=sdpvar(n,n,'f'); 
G1=sdpvar(n,n,'f'); 

%% Notations
tauTilde=h+muM; 
rhoTilde=exp(-2*alpha*(r1+tauTilde)); 
rhoM=exp(-2*alpha*(r1+muM)); 

%% The LMI for Sigma 
Sigma=blkvar; 
Sigma(1,1)=2*alpha*P+S+P2'*(A+B*K)+(A+B*K)'*P2; 
Sigma(1,2)=P-P2'+(A+B*K)'*P3; 
Sigma(1,3)=-P2'*expm(A*r1)*B*K; 
Sigma(1,5)=P2'*expm(A*r1)*B*K;
Sigma(2,2)=muM^2*R0+h^2*R1-P3-P3'; 
Sigma(2,3)=-P3'*expm(A*r1)*B*K; 
Sigma(2,5)=P3'*expm(A*r1)*B*K; 
Sigma(3,3)=exp(-2*alpha*r1)*(S0-S)-rhoM*R0; 
Sigma(3,4)=rhoM*R0; 
Sigma(4,4)=-rhoM*(R0+S0-S1)-rhoTilde*R1; 
Sigma(4,5)=rhoTilde*(R1-G1); 
Sigma(4,6)=rhoTilde*G1; 
Sigma(5,5)=-rhoTilde*(2*R1-G1-G1'); 
Sigma(5,6)=rhoTilde*(R1-G1); 
Sigma(6,6)=-rhoTilde*(S1+R1); 
Sigma=sdpvar(Sigma); 

%% The LMI for Xi 
Xi=blkvar; 
Xi(1,1)=2*alpha*P+S+P2'*(A+B*K)+(A+B*K)'*P2; 
Xi(1,2)=P-P2'+(A+B*K)'*P3; 
Xi(1,3)=-P2'*expm(A*r1)*B*K; 
Xi(1,4)=P2'*expm(A*r1)*B*K; 
Xi(1,7)=P2'*expm(A*r1)*B; 
Xi(2,2)=muM^2*R0+h^2*R1-P3-P3'; 
Xi(2,3)=-P3'*expm(A*r1)*B*K; 
Xi(2,4)=P3'*expm(A*r1)*B*K; 
Xi(2,7)=P3'*expm(A*r1)*B; 
Xi(3,3)=exp(-2*alpha*r1)*(S0-S)-rhoM*R0; 
Xi(3,4)=rhoM*(R0-G0); 
Xi(3,5)=rhoM*G0; 
Xi(4,4)=-rhoM*(2*R0-G0-G0')+sigma*K'*Omega*K; 
Xi(4,5)=rhoM*(R0-G0); 
Xi(5,5)=rhoM*(S1-S0-R0)-rhoTilde*R1; 
Xi(5,6)=rhoTilde*R1; 
Xi(6,6)=-rhoTilde*(S1+R1); 
Xi(7,7)=-Omega; 
Xi=sdpvar(Xi); 

%% Park's conditions 
Park0=[R0 G0; G0' R0]; 
Park1=[R1 G1; G1' R1];

%% Solution of LMIs
LMIs=[P>=0, S>=0, S0>=0, S1>=0, R0>=0, R1>=0, Omega>=0, Sigma<=0, Xi<=0, Park0>=0, Park1>=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    if min(primal)>=0 && primal(1)>0 
        OmegaVal=value(Omega); 
    end
else
    yalmiperror(sol.problem) 
end
