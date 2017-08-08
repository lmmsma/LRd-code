folder= 'Stored/';
load('b1000fsolem12variable_amplitudes.mat');
eval(['load ' folder 'eps=10_-7/jacfile alljacs']);
jacA_1000= alljacs{46};
eval(['load ' folder 'eps=10_-6/jacfile alljacs']);
jacB_1000= alljacs{46};
eval(['load ' folder 'eps=10_-5/jacfile alljacs']);
jacC_1000= alljacs{46};
eval(['load ' folder 'eps=10_-4/jacfile alljacs']);
jacD_1000= alljacs{46};
eval(['load ' folder 'eps=10_-3/jacfile alljacs']);
jacE_1000= alljacs{43};

diff1=jacB_1000-jacA_1000;
diff2=jacC_1000-jacB_1000; 
diff3=jacD_1000-jacC_1000;
diff4=jacE_1000-jacD_1000;
Smat = diag(1./varamp);
Smatinv = inv(Smat); % only need to compute once

df1= zeros(17,1); df2= zeros(17,1); df3= zeros(17,1); df4= zeros(17,1);

for i= 1:17
    df1(i)= sum(diff1(i , :));
    df2(i)= sum(diff2(i , :));
    df3(i)= sum(diff3(i , :));
    df4(i)= sum(diff4(i , :));
end

norm_df1= Smat*df1;
norm_df2= Smat*df2;
norm_df3= Smat*df3;
norm_df4= Smat*df4;
n1= norm(norm_df1)
n2= norm(norm_df2)
n3= norm(norm_df3)
n4= norm(norm_df4)