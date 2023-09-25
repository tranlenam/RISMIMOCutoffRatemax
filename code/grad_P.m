function y = grad_P(Hdir,P,H1,H2,RIS_phase,x,N0)
H = Hdir + H2*diag(RIS_phase)*H1;
no_symb = size(x,2);
s = 0; 
for iSym = 1:no_symb
    for jSym = iSym+1:no_symb
        eij = x(:,iSym)-x(:,jSym);
        s = s+exp(-norm(H*P*eij)^2/(4*N0))*(eij)*(eij');     
    end
end
s = 2*s;
y = -1/(4*N0)*(H')*H*P*s;
end 