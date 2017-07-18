% This code is adapted from the original compute_obsv_ctrb code. It will
% plot the observability and/or controllability plots for each
% measurement for both the adjusted and default parameter sets, as well as
% the results for these matrices after the multiplication of a scaling
% factor. The four plots are then placed on a single plot showing all four
% results for comparison
%% Observability Plots
clear variables;
oc_folder = 'OCValues/';
eval(['load ' oc_folder '/def/ocfile *']);
def.evof = eigof.oval;
def.evnof = eigof.noval;
def.evof_scaled = eigof.oval_sc;
def.evnof_scaled = eigof.noval_sc;

eval(['load ' oc_folder '/adj/ocfile *']);
adj.evof = eigof.oval;
adj.evnof = eigof.noval;
adj.evof_scaled = eigof.oval_sc;
adj.evnof_scaled = eigof.noval_sc;

clearvars -except def adj rankcutoff bcls numstate statenames nbcls oc_folder;
def_adj = cell(1, 2);
def_adj{1} = def; def_adj{2} = adj; %Stores the def and adj parameter sets so they can be referenced within loop

% Cycle through possible outputs (measurements), unscaled system  
for kk = 1:numstate
    eval(['h' num2str(kk) '= figure;']);
    set(gcf, 'Position', [317 155 1102 781]);
    for adj_yn= 1:2
        if adj_yn==1
            param = 'Default';
        else
            param = 'Adjusted';
        end
        subplot(2, 2, adj_yn)
        title({[param, ': Observable eigenvalues'], ['Meas. = ' strtrim(statenames(kk,:)) ', rankcutoff = ' num2str(rankcutoff)]});
        ylabel('Eigenvalue magnitude');
        xlabel('BCL (ms)');
        grid on;
        hold on;
        for i=1:nbcls
            if ~isempty(def_adj{adj_yn}.evof{i,kk})
                p1 = plot(bcls(i)*ones(size(def_adj{adj_yn}.evof{i,kk})),abs(def_adj{adj_yn}.evof{i,kk}),'b*');
            end
            if ~isempty(def_adj{adj_yn}.evnof{kk})
                p2 = plot(bcls(i)*ones(size(def_adj{adj_yn}.evnof{i,kk})),abs(def_adj{adj_yn}.evnof{i,kk})','ro');
            end
            %legend([p1 p2],'observable','unobservable')
            %saveas(gcf,[oc_folder 'obsvf_eig_meas' num2str(kk)])
        end
        
        %Plots the Scaled versions on the same figure
        subplot(2, 2, adj_yn+2)
        title({[param, ': Scaled system, observable eig.'], ['Meas. = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]});
        ylabel('Eigenvalue magnitude');
        xlabel('BCL (ms)');
        grid on;
        hold on;
        for i=1:nbcls
            if ~isempty(def_adj{adj_yn}.evof_scaled{i,kk})
                p1 = plot(bcls(i)*ones(size(def_adj{adj_yn}.evof_scaled{i,kk})),abs(def_adj{adj_yn}.evof_scaled{i,kk}),'b*');
            end
            if ~isempty(def_adj{adj_yn}.evnof{kk})
                p2 = plot(bcls(i)*ones(size(def_adj{adj_yn}.evnof_scaled{i,kk})),abs(def_adj{adj_yn}.evnof_scaled{i,kk})','ro');
            end
            %legend([p1 p2],'observable','unobservable')
        end
    end
    saveas(gcf,[oc_folder 'obsvf_eig_meas_subplots' num2str(kk) '_E' num2str(log10(rankcutoff))])
    saveas(gcf,[oc_folder 'obsvf_eig_meas_subplots' num2str(kk) '_E' num2str(log10(rankcutoff)) '.jpeg'])
end


%% Controllability Plots (Omitted 

% % Having 17x4 = 68 figures open at once may cause problems, so close some
% % of them, if desired
% cfflag = input('Before generating controllability figures, enter 1 to close observability figures, or enter 0 to keep the figures: ');
% if cfflag
%     close all; 
% end
% 
% % Cycle through possible control inputs, unscaled system  
% for kk = 1:numstate
%     figure
%     title([param, ': Controllable eigenvalues, input = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]);
%     ylabel('Eigenvalue magnitude');
%     xlabel('BCL (ms)');
%     grid on;
%     hold on;
%     for i=1:nbcls
%         if ~isempty(evcf{i,kk})
%             p1 = plot(bcls(i)*ones(size(evcf{i,kk})),abs(evcf{i,kk}),'b*');
%         end
%         if ~isempty(evncf{kk})
%             p2 = plot(bcls(i)*ones(size(evncf{i,kk})),abs(evncf{i,kk})','r*');
%         end
%         legend([p1 p2],'controllable','uncontrollable')
%         saveas(gcf,[OCvalues 'ctrbf_eig_meas' num2str(kk)])
%     end
% end
% 
% % Cycle through possible control inputs, scaled system
% for kk = 1:numstate
%     figure
%     title([param, ': Scaled system, controllable eig., input = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]);
%     ylabel('Eigenvalue magnitude');
%     xlabel('BCL (ms)');
%     grid on;
%     hold on;
%     for i=1:nbcls
%         if ~isempty(evcf_scaled{i,kk})
%             p1 = plot(bcls(i)*ones(size(evcf_scaled{i,kk})),abs(evcf_scaled{i,kk}),'b*');
%         end
%         if ~isempty(evncf_scaled{kk})
%             p2 = plot(bcls(i)*ones(size(evncf_scaled{i,kk})),abs(evncf_scaled{i,kk})','r*');
%         end
%         legend([p1 p2],'controllable','uncontrollable')
%         saveas(gcf,[OCvalues 'ctrbf_scaled_eig_meas' num2str(kk)])
%     end
% end