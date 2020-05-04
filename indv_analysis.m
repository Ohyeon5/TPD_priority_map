function indv_analysis(subjID, isSaveFig)
    % quick overview of the data, combine all the data from same session
    % and get pc and rt results depends on the conditions
    % usage: ex. indv_analysis('OO', 1) 
    %      : Analyze subject 'OO's data and save the bar plot for each session
    
    % get all tested conditions
    data_path = ['..' filesep 'Data' filesep upper(subjID)];
    f_list    = dir([data_path filesep upper(subjID) '*.dv']);
    data_fn   = {};
    ses_nm    = {};
    % get unique session name
    ses_unm = {};
    for ii=1:length(f_list)        
        data_fn{ii,1} = [data_path filesep f_list(ii).name];
        tempSplit = strsplit(data_fn{ii,1},'_');
        ses_nm{ii,1}  = strjoin(tempSplit(3:end-2),'_');
        if ~ismember(ses_nm{ii},ses_unm)
            ses_unm{length(ses_unm)+1} = ses_nm{ii};
        end
    end
    
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
        anal_pc = [];
        anal_rt = [];
        for jj = 1:length(conds)
            anal_pc(jj) = mean(ses_dt.iscorr(ses_dt.condition == conds(jj)&ses_dt.keyidx~=0));
            anal_rt(jj) = mean(ses_dt.rt(ses_dt.condition == conds(jj)&ses_dt.keyidx~=0));
        end
        
        fig = figure('position',[100 100 length(conds)*250 400]);
        subplot(1,2,1);
        bar(anal_pc,0.8)
        text([1:length(conds)],anal_pc+0.05,strsplit(num2str(anal_pc)),'HorizontalAlignment', 'center')
        title('PC')
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
            saveas(fig, [data_path filesep ses_unm{ii} '.png'])
        end
    end % end of session_analysis
    
end % end of indv_analysis()