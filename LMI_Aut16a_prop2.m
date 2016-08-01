function OmegaVal=LMI_Aut16a_prop2(A,B,K,h,alpha,sigma)
% This MATLAB program checks the feasibility of LMIs from Proposition 2 of the paper 
% A. Selivanov and E. Fridman, "Predictor-based networked control under uncertain transmission delays," Automatica, vol. 70, pp. 101–108, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B   - the parameters of the system (3); 
% K     - the controller gain from (23); 
% h     - the maximum sampling period; 
% alpha - a desired decay rate; 
% sigma - the event-triggering parameter from (24); 

% Output: 
% OmegaVal - the value of Omega from (24). If OmegaVal is empty, the LMIs are not feasible

n=size(A,1); 
m=size(B,2); 

%% Decision variables 
P=sdpvar(n); 
S=sdpvar(n); 
R=sdpvar(n); 
Omega=sdpvar(m); 
P2=sdpvar(n,n,'f'); 
P3=sdpvar(n,n,'f'); 
G=sdpvar(n,n,'f'); 

%% Notations
rhoh=exp(-2*alpha*h); 

%% The LMI for M 
M=blkvar; 
M(1,1)=2*alpha*P+S-rhoh*R+P2'*A+A'*P2; 
M(1,2)=P-P2'+A'*P3; 
M(1,3)=rhoh*(R-G)+P2'*B*K; 
M(1,4)=rhoh*G; 
M(2,2)=h^2*R-P3-P3'; 
M(2,3)=P3'*B*K; 
M(3,3)=-rhoh*(2*R-G-G'); 
M(3,4)=rhoh*(R-G); 
M(4,4)=-rhoh*(S+R); 
M=sdpvar(M); 

%% The LMI for N 
N=blkvar; 
N(1,1)=2*alpha*P+S-rhoh*R+sigma*K'*Omega*K+P2'*(A+B*K)+(A+B*K)'*P2; 
N(1,3)=rhoh*R; 
N(1,4)=P2'*B; 
N(2,2)=h^2*R-P3-P3'; 
N(2,4)=P3'*B; 
N(3,3)=-rhoh*(S+R); 
N(4,4)=-Omega; 
N=sdpvar(N); 

%% Park's condition 
Park=[R G; G' R]; 

%% Solution of LMIs
LMIs=[P>=0, S>=0, R>=0, Omega>=0, M<=0, N<=0, Park>=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    if min(primal)>=0 && primal(1)>0 
        OmegaVal=double(Omega); 
    end
else
    yalmiperror(sol.problem); 
end