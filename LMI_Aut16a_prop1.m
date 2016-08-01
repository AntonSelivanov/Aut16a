function OmegaVal=LMI_Aut16a_prop1(A,B,K,h,etaM,alpha,sigma)
% This MATLAB program checks the feasibility of LMIs from Proposition 1 of the paper 
% A. Selivanov and E. Fridman, "Predictor-based networked control under uncertain transmission delays," Automatica, vol. 70, pp. 101–108, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B   - the parameters of the system (15); 
% K     - the controller gain from (7); 
% h     - the maximum sampling period; 
% etaM  - a bound for the sensors-to-controller time-varying delay uncertainty; 
% alpha - a desired decay rate; 
% sigma - the event-triggering parameter from (13); 

% Output: 
% OmegaVal - the value of Omega from (13). If OmegaVal is empty, the LMIs are not feasible

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
tauBar=h+etaM; 
rhoBar=exp(-2*alpha*tauBar); 

%% The LMI for Psi
Psi=blkvar; 
Psi(1,1)=2*alpha*P+S-rhoBar*R+P2'*A+A'*P2; 
Psi(1,2)=P-P2'+A'*P3; 
Psi(1,3)=rhoBar*(R-G)+P2'*B*K; 
Psi(1,4)=rhoBar*G; 
Psi(1,5)=P2'*B; 
Psi(2,2)=tauBar^2*R-P3-P3'; 
Psi(2,3)=P3'*B*K; 
Psi(2,5)=P3'*B; 
Psi(3,3)=-rhoBar*(2*R-G-G')+sigma*K'*Omega*K; 
Psi(3,4)=rhoBar*(R-G); 
Psi(4,4)=-rhoBar*(S+R); 
Psi(5,5)=-Omega; 
Psi=sdpvar(Psi); 

%% Park's condition 
Park=[R G; G' R];

%% Solution of LMIs
LMIs=[Psi<=0, P>=0, S>=0, R>=0, Omega>=0, Park>=0]; 
options=sdpsettings('solver','sedumi','verbose',0); 
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    if min(primal)>=0 && primal(2)>0 
        OmegaVal=double(Omega); 
    end
else
    yalmiperror(sol.problem); 
end