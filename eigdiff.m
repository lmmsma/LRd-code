file1= 'Eigenvalues/eigfile_red';;%select reduced def param eigfile
file2= 'Eigenvalues/eigfile_adj';%select adjusted parameter eigfile
load(file1, 'alleigsabs');
alleigs_def= alleigsabs;
load(file2, 'alleigsabs')
alleigs_adj= alleigsabs;
diff= zeros(17, length(alleigs_adj));
for i=1:length(alleigs_def)
   s1= sort(alleigs_def{i});
   s2= sort(alleigs_adj{i});
   diff(:, i)= s2-s1;
end

%plot the difference in Eigenvalues in a semi-log plot
for i=1:46
    for j=1:17
        scatter(bcls(i), log10(abs(diff(j, i))))
        hold on
    end
    xlabel('BCL (ms)');
    ylabel('log_1_0(difference)');
    title('Difference in Eigenvalues between parameter sets');
end