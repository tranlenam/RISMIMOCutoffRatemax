clear 
clc
rng(1)
Nt = 8;
Nr = 2;
Nris = 225;
lt = 20;
lr = 20;
D = 500;
no_mat = 1;             % number of channel realizations to be generated
K = 1;                  % Rician K factor 
f = 2000e6;             % frequency
dist_ris = 30;
N0 = db2pow(-110);      % noise power
M = 4;                  % modulatio order QPSK
alpha_dir = 3;          % path loss exponent of the direct link
Pow = Nr;               % preconding matrix power constrain 
no_iter = 40;           % number of iterations

Ns = (M^Nr);          % number of sympbols

% line search paramaters
delta = 1e-4;
rho = 1/2;
step_size_theta=1e5; % initial value for step size for phase shift
step_size_P=1e3; % initial value for step size for precoder

% initialization 
% channel matrices 
[Hdir,H1,H2] = generateChannels(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);

%% conventional PGM
% initial precoding matrix 
P0 = rand(Nt,Nr)+1i*rand(Nt,Nr);
P = P0/norm(P0,'fro')*sqrt(Nr);

% initial RIS phase shifts
mytheta0 = rand(Nris,1)+1i*rand(Nris,1);
mytheta = mytheta0./abs(mytheta0);

%% sequential projection


Hdirt = Hdir{1}; 
H1t = H1{1}; 
H2t = H2{1};
DIR = 1; % DIR = 1 if direct link present, DIR = 0 other wise 
% Remove the direct link if DIR = 0
if ~DIR
    Hdirt = zeros(size(Hdirt));
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = Hdirt + H2t*diag(mytheta)*H1t; 
x = symb_vec_set(M,size(H,1));


objSPGM = zeros(no_iter+1,1);
objSPGM_MI = zeros(no_iter+1,1);
objSPGM(1) = computeobjective(H*P*x,N0); % the initial objective
objSPGM_MI(1) = computeMI(H*P*x,N0);
for iIter =1:no_iter
    % line search over phase shift
    [mytheta,step_size_theta]=linesearch_RIS_Phaseshift(Hdirt,P,H1t,H2t,mytheta,x,N0,step_size_theta,delta,rho);
    Hnew = Hdirt + H2t*diag(mytheta)*H1t;
    
    % line search over precoder
    [P,step_size_P]=linesearch_Tx_Precoder(Hdirt,P,H1t,H2t,mytheta,x,N0,step_size_P,delta,rho,Pow);
    objSPGM(iIter+1) =  computeobjective(Hnew*P*x,N0);
    objSPGM_MI(iIter+1) =  computeMI(Hnew*P*x,N0);

end 
plot(1:no_iter+1,-log2(objSPGM/Ns^2),'b-','DisplayName','$R_0$, PGM');
hold on
plot(1:no_iter+1,objSPGM_MI,'r--','DisplayName','MI, PGM');
