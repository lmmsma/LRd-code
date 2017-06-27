bcls_tested = [1000:-50:600 590:-10:510 505 500:-10:70];%Vector of bcls
f=load('fp_found/compiled_fp');
fp= f.fp_found;

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');

for i = 1:17
    h(i) = figure;
    hold on;
end

for i = 1:17
    figure(h(i));
    xlabel('BCL (ms)');
    ylabel(statenames(i,:));
    title([statenames(i,:) 'vs BCL'])
    plot(bcls_tested, fp(i, :))
end