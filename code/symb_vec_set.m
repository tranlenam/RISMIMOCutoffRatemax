function symb_vec = symb_vec_set(M,len)
s = qammod(0:M-1,M);
symb_ind_comb = combinator(M,len,'p','r')';
symb_vec = s(symb_ind_comb);
end