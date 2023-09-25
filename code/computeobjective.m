% calculate the cut-off rate 
function obj = computeobjective(x,N0)
no_symb = size(x,2);
s = 0;
for iSym = 1:no_symb
    x1 = x(:,iSym)*ones(1,no_symb-iSym)-x(:,iSym+1:end);
    arg = sum(abs(x1).^2,1)/(4*N0);
    s = s + sum(exp(-arg));
end
obj = 2*s+no_symb;
end