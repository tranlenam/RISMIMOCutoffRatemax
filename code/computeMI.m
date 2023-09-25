function MI = computeMI(x,N0)
no_symb = size(x,2);
Nr = size(x,1);
nNoiseSamples = 5000; % the number of noise samples
s = 0;
for iSymbol = 1:no_symb
    NoiseSamples = sqrt(N0/2)*(randn(Nr,nNoiseSamples)+1i*randn(Nr,nNoiseSamples));
    stemp = zeros(1,nNoiseSamples);
    x1 = x(:,iSymbol)*ones(1,no_symb)-x;
    for iNoiseSample=1:nNoiseSamples
        x1it = x1+NoiseSamples(:,iNoiseSample)*ones(1,no_symb);
        arg = (sum(abs(x1it).^2,1)-norm(NoiseSamples(:,iNoiseSample))^2)/N0;
        stemp(iNoiseSample) = log2(sum(exp(-arg)));
    end
    s = s+mean(stemp);
end
MI = log2(no_symb)-1/no_symb*s;
end
