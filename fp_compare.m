f1=load('fp_found/compiled_fp');
fp1= f1.fp_found;
f1=load('fp_found2/compiled_fp');
fp2= f1.fp_found;
short_bcls= 400:-10:70;
shortfp1= fp1(:,30:end);
shortfp2= fp2(:, 13:end);
fp_diff= shortfp2-shortfp1;
statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');

%plot diff for V on Default parameters
V_fp_diff1= diff(shortfp1(1, :));
scatter(short_bcls(1:33), V_fp_diff1);
title('V Diff for Def Parameter');
figure
V_fp_diff2= diff(shortfp2(1, :));
scatter(short_bcls(1:33), V_fp_diff2);
title('V Diff for Alt Parameter');
figure
Ca_fp_diff1= diff(shortfp1(6, :));
scatter(short_bcls(1:33), Ca_fp_diff1);
title('Ca Diff for Def Parameter');
figure
Ca_fp_diff2= diff(shortfp2(6, :));
scatter(short_bcls(1:33), Ca_fp_diff2);
title('Ca Diff for Alt Parameter');



for i= [1 6]
    h(i) = figure;
    hold on;
end

for i= [1 6]
    figure(h(i));
    xlabel('BCL (ms)');
    ylabel(statenames(i,:));
    title([statenames(i,:) 'difference vs BCL'])
    scatter(short_bcls, fp_diff(i, :))
end