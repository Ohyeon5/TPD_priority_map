function fig = fc_plot_pe_rt(anal_pe, anal_rt, sess_names, condNms, conds, fig_title)
    % Plot percent error and response time difference &/ themselves
    
    % parameters
    % ----------
    %  anal_pe : array [nSubj x nSess x nCond] or [nSess x nCond] when nSubj==1
    
    colorMaps = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
    if length(size(anal_pe))==2
        nSubj = 1; 
    else
        nSubj = size(anal_pe,1);
    end


    rr = 2;
    cc = 2;
    fig = figure('position',[100 100 cc*length(sess_names)*200 rr*400]);
    
    plot_pe = anal_pe;
    plot_rt = anal_rt;
    diff_pe = [];
    diff_rt = [];
    
    % compute difference from baseline 
    if length(conds)==3
        if nSubj == 1   % individual plot
            diff_pe(:,:)=anal_pe(:,1:2) - repmat(anal_pe(:,3),1,2);
            diff_rt(:,:)=anal_rt(:,1:2) - repmat(anal_rt(:,3),1,2);
            plot_pe = diff_pe;
            plot_rt = diff_rt;
        else            % group plot
            for sess = 1:length(sess_names)
                plot_pe(:,sess,:) = anal_pe(:,sess,:) - anal_pe(:,sess,3);
                plot_rt(:,sess,:) = anal_rt(:,sess,:) - anal_rt(:,sess,3);
            end
            plot_pe(:,:,3) = [];
            plot_rt(:,:,3) = [];
        end
    end
    
    subplot(rr,cc,1);
    if nSubj > 1  % if plotting group results
        mean_pe = squeeze(mean(plot_pe,1));
        errorbar((repmat([1:length(sess_names)],2,1)+[-0.14;0.14])',mean_pe,squeeze(std(plot_pe,1))/sqrt(nSubj),'.k')
        b = bar(mean_pe, 0.8);
    else          % For the individual result
        mean_pe = plot_pe;
        b = bar(mean_pe,0.8);
    end
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    % draw each value on the graph directly
    text([1:length(sess_names)]-0.2,mean_pe(:,1)+0.01,strsplit(num2str(mean_pe(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(sess_names)]+0.2,mean_pe(:,2)-0.01,strsplit(num2str(mean_pe(:,2)',2)),'HorizontalAlignment', 'center')
    title([ fig_title ' - Percent error difference'])
    legend(condNms(conds));
    if length(conds)~=3     
        ylim([0,1])   
        ylabel('percent error')
        xticklabels(sess_names)
    else
        ylim([-0.2,0.2])
        ylabel('percent error difference')
        % show base line performance as xticks
        for ii=1:length(sess_names)  
            if nSubj==1; base_pe = anal_pe(ii,3); else; base_pe = squeeze(mean(anal_pe(:,ii,3))); end
            temp_labels{ii} = [sess_names{ii} ': ' num2str(base_pe,2)];
        end
        xticklabels(temp_labels)
    end  
    % plot RT difference
    subplot(rr,cc,2);
    if nSubj > 1  % if plotting group results
        mean_rt = squeeze(mean(plot_rt,1));
        errorbar((repmat([1:length(sess_names)],2,1)+[-0.14;0.14])',mean_rt,squeeze(std(plot_rt,1))/sqrt(nSubj),'.k')
        b = bar(mean_rt, 0.8);
    else          % For the individual result
        mean_rt = plot_rt;
        b = bar(mean_rt,0.8);
    end   
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    text([1:length(sess_names)]-0.2,mean_rt(:,1)+0.01,strsplit(num2str(mean_rt(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(sess_names)]+0.2,mean_rt(:,2)-0.01,strsplit(num2str(mean_rt(:,2)',2)),'HorizontalAlignment', 'center')
    title(['Reaction time difference'])
    if length(conds)~=3     
        ylim([0,1])   
        ylabel('RT (s)')
        xticklabels(sess_names)
    else
        ylim([-0.2,0.2])
        ylabel('RT difference (s)')
        % show base line performance as xticks
        for ii=1:length(sess_names)
            if nSubj==1; base_rt = anal_rt(ii,3); else; base_rt = squeeze(mean(anal_rt(:,ii,3))); end
            temp_labels{ii} = [sess_names{ii} ': ' num2str(base_rt,2)];
        end
        xticklabels(temp_labels)
    end
    
    subplot(rr,cc,3);
    if nSubj > 1  % if plotting group results
        mean_pe = squeeze(mean(anal_pe,1));
        errorbar((repmat([1:length(sess_names)],3,1)+[-0.14;0;0.14])',mean_pe,squeeze(std(anal_pe,1)),'.k')
        b = bar(mean_pe, 0.8);
    else          % For the individual result
        mean_pe = anal_pe;
        b = bar(mean_pe,0.8);
    end
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    b(3).FaceColor = colorMaps(3,:);
    text([1:length(sess_names)]-0.2,mean_pe(:,1)+0.01,strsplit(num2str(mean_pe(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(sess_names)]    ,mean_pe(:,2)+0.01,strsplit(num2str(mean_pe(:,2)',2)),'HorizontalAlignment', 'center')    
    text([1:length(sess_names)]+0.2,mean_pe(:,3)+0.01,strsplit(num2str(mean_pe(:,3)',2)),'HorizontalAlignment', 'center')
    title(['Percent error'])
    legend(condNms);
    ylim([0,max(mean_pe(:))+0.3])   
    ylabel('percent error')
    xticklabels(sess_names)
        
    subplot(rr,cc,4);
    if nSubj > 1  % if plotting group results
        mean_rt = squeeze(mean(anal_rt,1));
        errorbar((repmat([1:length(sess_names)],3,1)+[-0.14;0;0.14])',mean_rt,squeeze(std(anal_rt,1)),'.k')
        b = bar(mean_rt, 0.8);
    else          % For the individual result
        mean_rt = anal_rt;
        b = bar(mean_rt,0.8);
    end
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    b(3).FaceColor = colorMaps(3,:);
    text([1:length(sess_names)]-0.2,mean_rt(:,1)+0.01,strsplit(num2str(mean_rt(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(sess_names)]    ,mean_rt(:,2)+0.01,strsplit(num2str(mean_rt(:,2)',2)),'HorizontalAlignment', 'center')    
    text([1:length(sess_names)]+0.2,mean_rt(:,3)+0.01,strsplit(num2str(mean_rt(:,3)',2)),'HorizontalAlignment', 'center')
    title(['Reaction time'])
    ylim([0.5,1.5])   
    ylabel('RT (s)')
    xticklabels(sess_names)

end % end of fc_plot_pc_rt()