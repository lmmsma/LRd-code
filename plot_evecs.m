clear variables;
oc_folder = 'OCValues/';

eval(['load ' oc_folder '/def/ocfile *']);
def.sum_ov= zeros(numstate, 1);
def.sum_nov= zeros(numstate, 1);
def.sum_cv= zeros(numstate, 1);
def.sum_ncv= zeros(numstate, 1);
def.sum_ov_sc= zeros(numstate, 1);
def.sum_nov_sc= zeros(numstate, 1);
def.sum_cv_sc= zeros(numstate, 1);
def.sum_ncv_sc= zeros(numstate, 1);
def.eigof= eigof;
def.eigcf= eigcf;

eval(['load ' oc_folder '/adj/ocfile *']);
adj.sum_ov= zeros(numstate, 1);
adj.sum_nov= zeros(numstate, 1);
adj.sum_cv= zeros(numstate, 1);
adj.sum_ncv= zeros(numstate, 1);
adj.sum_ov_sc= zeros(numstate, 1);
adj.sum_nov_sc= zeros(numstate, 1);
adj.sum_cv_sc= zeros(numstate, 1);
adj.sum_ncv_sc= zeros(numstate, 1);
adj.eigof= eigof;
adj.eigcf= eigcf;

def_adj = cell(1, 2);
def_adj{1} = def; def_adj{2} = adj; %Stores the def and adj parameter sets so they can be referenced within loop

for i= 1:nbcls
    for adj_yn= 1:2
        for j=1:numstate
            def_adj{adj_yn}.sum_ov = def_adj{adj_yn}.sum_ov + sum(abs(def_adj{adj_yn}.eigof.ovec{i, j}), 2);
            def_adj{adj_yn}.sum_nov= def_adj{adj_yn}.sum_nov - sum(abs(eigof.novec{i, j}), 2);
            def_adj{adj_yn}.sum_ov_sc = def_adj{adj_yn}.sum_ov_sc + sum(abs(def_adj{adj_yn}.eigof.ovec_sc{i, j}), 2);
            def_adj{adj_yn}.sum_nov_sc= def_adj{adj_yn}.sum_nov_sc - sum(abs(def_adj{adj_yn}.eigof.novec_sc{i, j}), 2);
%           def_adj{adj_yn}.sum_cv = def_adj{adj_yn}.sum_cv + sum(abs(def_adj{adj_yn}.eigcf.cvec{i, j}), 2);
%           def_adj{adj_yn}.sum_ncv= def_adj{adj_yn}.sum_ncv - sum(abs(def_adj{adj_yn}.eigcf.ncvec{i, j}), 2);
        end
    end
end


%Plot the sum of the maagnitude of each variable in the Observable and
%Unobservable Eigenvectors of the system

fo= figure('Name', 'Observability Vector mags', 'NumberTitle', 'off');
set(gcf, 'Position', [317 87 1251 849]);
for adj_yn= 1:2
    if adj_yn==1
        param = 'Default';
    else
        param = 'Adjusted';
    end
    subplot(2, 2, adj_yn)
    plot(def_adj{adj_yn}.sum_ov);
    set(gca,'XTick',1:17,'XTickLabel',statenames)
    title({[param, ': Scaled obs eigenvec mags'], ['rankcutoff = ' num2str(rankcutoff)]});
    xlabel('State Variable');
    ylabel('Sum of Magnitudes');
    grid
    hold on
    plot(def_adj{adj_yn}.sum_nov);
    legend('Observable', 'Unobservable');
    
    %plot scaled versions on same figure
    subplot(2, 2, adj_yn+2)
    plot(def_adj{adj_yn}.sum_ov_sc);
    set(gca,'XTick',1:17,'XTickLabel',statenames)
    title({[param, ': Scaled obs eigenvec mags'], ['rankcutoff = ' num2str(rankcutoff)]});
    xlabel('State Variable');
    ylabel('Sum of Magnitudes');
    grid
    hold on
    plot(def_adj{adj_yn}.sum_nov_sc);
    legend('Observable', 'Unobservable');
end

saveas(fo,[oc_folder 'obsvf_eig_vec' '_E' num2str(log10(rankcutoff))])
saveas(fo,[oc_folder 'obsvf_eig_vec' '_E' num2str(log10(rankcutoff)) '.jpeg'])



% %Plot the sum of the maagnitude of each variable in the Controallable and
% %Uncontrollable Eigenvectors of the system
% fv= figure('Name', 'Controllability Vector mags', 'NumberTitle', 'off');
% set(gcf, 'Position', [111 229 754 564]);
% plot(sum_cv);
% set(gca,'XTick',1:17,'XTickLabel',statenames)
% title({'Controllability Eigenvector Magnitudes', 'Rank Cutoff =' rankcutoff});
% xlabel('State Variable');
% ylabel('Sum of Magnitudes');
% hold on
% plot(sum_ncv);
% legend('Observable', 'Unobservable');
% saveas(fv,[oc_folder 'ctrbf_eig_vec' '_E' num2str(log10(rankcutoff))])
% saveas(fv,[oc_folder 'ctrbf_eig_vec' '_E' num2str(log10(rankcutoff)) '.jpeg'])
% 


