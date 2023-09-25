clear
clc
rng(1)
Nt = 8;
Nr = 2;
Nris = 225;
lt = 20;
lr = 20;
D = 500;
no_mat = 1;             % number of channel realizations
K = 1;                  % Rician K factor 
f = 2000e6;             % frequency
dist_ris = 30;
N0 = db2pow(-110);      % noise power
M = 4;                  % modulation order QPSK
alpha_dir = 3;          % path loss exponent of the direct link
Pow = Nr;               % preconding matrix power constrain 
nSCAiter = 20;
Ns = M^Nr;

% initialization 
% channel matrices 
[Hdir,H1,H2] = generateChannels(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);

Hdirt = Hdir{1}; 
H1t = H1{1}; 
H2t = H2{1};

DIR = 1; % DIR = 1 if direct link present, DIR = 0 other wise 
% Remove the direct link if DIR = 0
if ~DIR
    Hdirt = zeros(size(Hdirt));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = symb_vec_set(M,size(Hdirt,1));

% precoding matrix 
P = rand(Nt,Nr)+1i*rand(Nt,Nr);
P = P/norm(P,'fro')*sqrt(Nr);

% RIS phase shifts
mytheta = rand(Nris,1) +1i*rand(Nris,1);
mytheta = mytheta./abs(mytheta);


mytheta_variable = sdpvar(Nris,1,'full','complex');
P_variable = sdpvar(Nt,Nr,'full','complex');
solversettings = sdpsettings('solver','mosek','verbose',0);

Hnew = Hdirt + H2t*diag(mytheta)*H1t; % the channel with current theta
objSCA = zeros(nSCAiter+1,1);
objSCA_MI = zeros(nSCAiter+1,1);
objSCA(1) = computeobjective(Hnew*P*x,N0); % the initial objective
objSCA_MI(1) = computeMI(Hnew*P*x,N0); % the initial objective

for iSCA=1:nSCAiter
    % update P_variable
    %{
    cvx_begin quiet
    variable P_variable(Nt,Nr) complex
    obj = 0;
    for ind1 = 1:Ns
        for ind2 = ind1+1:Ns
            eij = x(:,ind1)-x(:,ind2);
            myPhi = norm(Hnew*P*eij,'fro')^2+2*real(trace((Hnew'*Hnew*P*(eij*eij'))'*(P_variable-P)));
            obj  = obj+exp(-myPhi/(4*N0));
        end
    end
    minimize(obj);
    subject to
    norm(P_variable,'fro') <= sqrt(Nr);
    cvx_end
    P=P_variable;
    %}
    Hnew = Hdirt + H2t*diag(mytheta)*H1t; % update channel with newly obtained theta

    obj = 0;
    for ind1 = 1:Ns
        for ind2 = ind1+1:Ns
            eij = x(:,ind1)-x(:,ind2);
            %myPhi = norm(Hnew*P*eij,'fro')^2+2*real(trace((Hnew'*Hnew*P*(eij*eij'))'*(P_variable-P)));
            myPhi = norm(Hnew*P*eij,'fro')^2+2*real(eij'*P'*(Hnew')*Hnew*(P_variable-P)*eij);
            obj  = obj+exp(-myPhi/(4*N0));
        end
    end
    diagnotics = optimize(norm(P_variable,'fro') <= sqrt(Nr),obj,solversettings);
    if(diagnotics.problem==0)
        P = value(P_variable);
    end

    % objSCA = [objSCA double(obj)];

    % update theta
    %{
    cvx_begin quiet
    variable mytheta_variable(Nris,1) complex
    obj = 0;
    for ind1 = 1:Ns
        for ind2 = ind1+1:Ns
            eij = x(:,ind1)-x(:,ind2);
            gradtheta = diag(H2t'*Hnew*P*(eij*eij')*P'*H1t');
            myPhi = norm(Hnew*P*eij,'fro')^2+2*real((gradtheta)'*(mytheta_variable-mytheta));
            obj  = obj+exp(-myPhi/(4*N0));
        end
    end
    minimize(obj);
    subject to
    abs(mytheta_variable) <= 1;
    cvx_end
    %}

    obj = 0;
    for ind1 = 1:Ns
        for ind2 = ind1+1:Ns
            eij = x(:,ind1)-x(:,ind2);
            gradtheta = diag(H2t'*Hnew*P*(eij*eij')*P'*H1t');
            myPhi = norm(Hnew*P*eij,'fro')^2+2*real((gradtheta)'*(mytheta_variable-mytheta));
            obj  = obj+exp(-myPhi/(4*N0));
        end
    end
    diagnotics = optimize(abs(mytheta_variable)<=1,obj,solversettings);
    if(diagnotics.problem==0)
        mytheta = value(mytheta_variable);
    end
    objSCA(iSCA+1) = 2*double(obj)+Ns;
    objSCA_MI(iSCA+1) = computeMI((Hdirt + H2t*diag(mytheta)*H1t)*P*x,N0);
end
plot(1:length(objSCA),-log2(objSCA/Ns^2),'r','DisplayName','$R_0$, SCA');
hold on 
plot(1:length(objSCA_MI),objSCA_MI,'r--','DisplayName','MI, SCA');
xlabel('Iteration'); ylabel('Cut-off rate / MI [bpcu]');
legend('show','Location','SouthEast','Interpreter','latex');
saveas(gcf,'../results/Fig2a.png')

