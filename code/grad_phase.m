function y = grad_phase(Hdir,P,H1,H2,RIS_phase,x,N0)
H = Hdir + H2*diag(RIS_phase)*H1; 
no_symb = size(x,2);
s = 0; 
for i = 1:no_symb
    for j = i+1:no_symb
        eij = x(:,i)-x(:,j);
        s = s+exp(-norm(H*P*eij)^2/(4*N0))*diag(H2'*H*P*eij*eij'*P'*H1');     
    end
end
s = 2*s;
y = -1/(4*N0)*s;
end 