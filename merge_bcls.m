clear variables;
jacfolder= 'Jacobians/';
fullfilename= char(strcat(jacfolder, 'jacfile'));

eval(['load ' jacfolder 'jacfile_def *']);
jacfolder= 'Jacobians/';
eval(['load ' jacfolder 'pacedownsettings data']);
new_bcls= [1000:-50:400 390:-10:70];
new_alljacs = cell(1,length(new_bcls));

i=1; j=1;

while i<=length(bcls)
   if bcls(i)==new_bcls(j)
       new_alljacs{j} = alljacs{i};

       j = j+1;
   end
   i = i+1;
end
bcls= new_bcls;
alljacs= new_alljacs;
% allv= new_allv;
% alleigs= new_eigs;
% alleigsabs= new_eigsabs;
clear new_alljacs new_bcls new_allv new_eigs new_eigsabs Eigenvalues Jacobians L bcl eigfolder folder fp fname epsln eigv subdiv_per_;

% save(fullfilename, 'bcls', 'alljacs', 'allv', 'alleigs', 'alleigsabs', 'modelname', 'data')
save(fullfilename, '*')