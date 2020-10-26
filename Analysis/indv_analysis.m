function results = indv_analysis(subjID, isSaveAnal)
    % quick overview of the data, combine all the data from same session
    % and get pc and rt results depends on the conditions
    % usage: ex. indv_analysis('OO', 1) 
    %      : Analyze subject 'OO's data and save the bar plot for each session
    
    % get all tested conditions
    data_path = ['..' filesep 'Data' filesep 'current' filesep upper(subjID)];
    f_list    = dir([data_path filesep upper(subjID) '*.dv']);
    ses_nm    = {};
    % get unique session name
    ses_unm = {};
    ses_unm_lab = {};
    for ii=1:length(f_list)        
        tempSplit = strsplit(f_list(ii).name,'_');
        ses_nm{ii,1}  = strjoin(tempSplit(3:end-2),'_');
        if ~ismember(ses_nm{ii},ses_unm)&&~isempty(ses_nm{ii})&&any(strfind(ses_nm{ii},'main'))
            ses_unm{length(ses_unm)+1} = ses_nm{ii};
            ses_unm_lab{length(ses_unm_lab)+1} = ses_nm{ii}(12:end);
        end
    end
    
    results = [];
    anal_pe = [];
    anal_rt = [];
    
    % analysis each session
    for ii = 1:length(ses_unm)
        ses_fn = dir([data_path filesep upper(subjID) '*' ses_unm{ii} '*.dv']);
        ses_dt = [];
        for jj=1:length(ses_fn)     
            ses_dv = lpsy.readDvFile([data_path filesep ses_fn(jj).name]);   % load dv file
            ses_dt = vertcat(ses_dt,struct2table(ses_dv.pool0));              % read data table only
        end       
        
        conds  = unique(ses_dt.condition);  % condition types
        condNms= {'Suppression', 'Enhancement', 'Baseline'};
        colorMaps = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
        
        for jj = 1:length(conds)
            anal_pe(ii,jj) = 1 - mean(ses_dt.iscorr(ses_dt.condition == conds(jj)...
                                                & ses_dt.keyidx~=0));
            anal_rt(ii,jj) = mean(ses_dt.rt(ses_dt.condition == conds(jj)...
                                            & ses_dt.keyidx~=0));
        end
        
        results(ii).sessNm = ses_unm{ii};
        results(ii).table  = ses_dt;
        
    end % end of session info collect
    plot_pe = anal_pe;
    plot_rt = anal_rt;
    diff_pe = [];
    diff_rt = [];
    % compute difference from baseline 
    if length(conds)==3
        diff_pe(:,:)=anal_pe(:,1:2) - repmat(anal_pe(:,3),1,2);
        diff_rt(:,:)=anal_rt(:,1:2) - repmat(anal_rt(:,3),1,2);
        plot_pe = diff_pe;
        plot_rt = diff_rt;
    end

    rr = 2;
    cc = 2;
    fig = figure('position',[100 100 cc*length(ses_unm)*200 rr*400]);
    
    subplot(rr,cc,1);
    b = bar(plot_pe,0.8);
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    text([1:length(ses_unm)]-0.2,plot_pe(:,1)+0.01,strsplit(num2str(plot_pe(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(ses_unm)]+0.2,plot_pe(:,2)-0.01,strsplit(num2str(plot_pe(:,2)',2)),'HorizontalAlignment', 'center')
    title(['Subject : ' subjID ' - Percent error difference'])
    legend(condNms(conds));
    if length(conds)~=3     
        ylim([0,1])   
        ylabel('percent error')
        xticklabels(ses_unm_lab)
    else
        ylim([-0.2,0.2])
        ylabel('percent error difference')
        for ii=1:length(ses_unm)
            temp_labels{ii} = [ses_unm_lab{ii} ': ' num2str(anal_pe(ii,3),2)];
        end
        xticklabels(temp_labels)
    end  
    subplot(rr,cc,2);
    b = bar(plot_rt,0.8);    
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    text([1:length(ses_unm)]-0.2,plot_rt(:,1)+0.01,strsplit(num2str(plot_rt(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(ses_unm)]+0.2,plot_rt(:,2)-0.01,strsplit(num2str(plot_rt(:,2)',2)),'HorizontalAlignment', 'center')
    title(['Reaction time difference'])
    if length(conds)~=3     
        ylim([0,1])   
        ylabel('RT (s)')
        xticklabels(ses_unm_lab)
    else
        ylim([-0.2,0.2])
        ylabel('RT difference (s)')
        for ii=1:length(ses_unm)
            temp_labels{ii} = [ses_unm_lab{ii} ': ' num2str(anal_rt(ii,3),2)];
        end
        xticklabels(temp_labels)
    end
    
    subplot(rr,cc,3);
    b = bar(anal_pe,0.8);
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    b(3).FaceColor = colorMaps(3,:);
    text([1:length(ses_unm)]-0.2,anal_pe(:,1)+0.01,strsplit(num2str(anal_pe(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(ses_unm)]    ,anal_pe(:,2)+0.01,strsplit(num2str(anal_pe(:,2)',2)),'HorizontalAlignment', 'center')    
    text([1:length(ses_unm)]+0.2,anal_pe(:,3)+0.01,strsplit(num2str(anal_pe(:,3)',2)),'HorizontalAlignment', 'center')
    title(['Percent error'])
    legend(condNms);
    ylim([0,max(anal_pe(:))+0.3])   
    ylabel('percent error')
    xticklabels(ses_unm_lab)
        
    subplot(rr,cc,4);
    b = bar(anal_rt,0.8);
    b(1).FaceColor = colorMaps(1,:);
    b(2).FaceColor = colorMaps(2,:);
    b(3).FaceColor = colorMaps(3,:);
    text([1:length(ses_unm)]-0.2,anal_rt(:,1)+0.01,strsplit(num2str(anal_rt(:,1)',2)),'HorizontalAlignment', 'center')
    text([1:length(ses_unm)]    ,anal_rt(:,2)+0.01,strsplit(num2str(anal_rt(:,2)',2)),'HorizontalAlignment', 'center')    
    text([1:length(ses_unm)]+0.2,anal_rt(:,3)+0.01,strsplit(num2str(anal_rt(:,3)',2)),'HorizontalAlignment', 'center')
    title(['Reaction time'])
    ylim([0.5,1.5])   
    ylabel('RT (s)')
    xticklabels(ses_unm_lab)
        
    if isSaveAnal
        saveas(fig, [data_path filesep subjID '_ind_anal.png'])
    end

    fprintf(['Done: ' subjID])
end % end of indv_analysis()