function indv_analysis(subjID, isSaveFig)
    % run short analysis on each participant and plot the percent correct results
    % usage: ex. indv_analysis('OO', 1) 
    %      : Analyze subject 'OO's data and save the bar plot for each session
    
    % get all tested conditions
    data_path = ['..' filesep 'Data' filesep upper(subjID)];
    f_list     = dir(data_path);
    data_fn    = {};
    for ii=1:length(f_list)        
        if (~f_list(ii).isdir) && any(strfind(f_list(ii).name,'.dv'))
            data_fn{ii,1} = [data_path filesep f_list(ii).name]; 
        end
    end
    data_fn = data_fn(~cellfun('isempty',data_fn));
    
    % analysis each session
    for ii = 1:length(data_fn)
        ses_fn = data_fn{ii,1};
        tempSplit = strsplit(ses_fn,'_');
        ses_nm = cell2mat(tempSplit(3:end-3));
        runNum = strsplit(tempSplit{end},'.');
        runNum = runNum{1};
        
        ses_dt = lpsy.readDvFile(ses_fn);   % load dv file
        ses_dt = ses_dt.pool0;              % read data table only
        
        conds  = unique(ses_dt.condition);  % condition types
        condNms= {'Suppression', 'Enhancement', 'Baseline'};
        anal_pc = [];
        anal_rt = [];
        for jj = 1:length(conds)
            anal_pc(jj) = mean(ses_dt.iscorr(ses_dt.condition == conds(jj)));
            anal_rt(jj) = mean(ses_dt.rt(ses_dt.condition == conds(jj)));
        end
        
        fig = figure('position',[100 100 length(conds)*250 400]);
        subplot(1,2,1);
        bar(anal_pc,0.8)
        text([1:length(conds)],anal_pc+0.05,strsplit(num2str(anal_pc)),'HorizontalAlignment', 'center')
        title(['PC session: ' ses_nm])
        ylabel('percent correct')
        ylim([0,1])
        xticklabels(condNms(conds))
        xtickangle(30)
        subplot(1,2,2);
        bar(anal_rt,0.8)
        text([1:length(conds)],anal_rt+0.2,strsplit(num2str(anal_rt)),'HorizontalAlignment', 'center')
        title(['RT'])
        ylabel('reaction time (s)')
        ylim([0,3])
        xticklabels(condNms(conds))
        xtickangle(30)
        
        if isSaveFig
            saveas(fig, [data_path filesep ses_nm '_RUN_' runNum '.png'])
        end
    end % end of session_analysis
    
end % end of indv_analysis()