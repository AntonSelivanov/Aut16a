% This MATLAB program checks the feasibility of LMIs from Theorems 1, 2 and Propositions 1, 2 of the paper 
% A. Selivanov and E. Fridman, "Predictor-based networked control under uncertain transmission delays," Automatica, vol. 70, pp. 101–108, 2016.

%% System parameters
M=10;   % the cart mass
m=1;    % the pendulum mass
l=3;    % the length of the pendulum arm
g=10;   % the gravitational acceleration

A=[0 1 0 0; 0 0 -m*g/M 0; 0 0 0 1; 0 0 g/l 0]; 
B=[0; 1/M; 0; -1/(M*l)]; 
K=[2 12 378 210]; 

%% LMIs of Theorem 1 and Proposition 1 
h=.0315; r0=.2; etaM=.01; r1=.2; muM=.01; alpha=.01; sigma=.01; 
display(['Theorem 1: Omega=' num2str(LMI_Aut16a_th1(A,B,K,h,r0,etaM,r1,muM,alpha,sigma))]); 
display(['Proposition 1: Omega=' num2str(LMI_Aut16a_prop1(A,B,K,h,etaM,alpha,sigma))]); 

%% LMIs of Theorem 2 and Proposition 2 
h=.104; r1=.2; muM=.01; alpha=.01; sigma=.13; 
display(['Theorem 2: Omega=' num2str(LMI_Aut16a_th2(A,B,K,h,r1,muM,alpha,sigma))]); 
display(['Proposition 2: Omega=' num2str(LMI_Aut16a_prop2(A,B,K,h,alpha,sigma))]); 