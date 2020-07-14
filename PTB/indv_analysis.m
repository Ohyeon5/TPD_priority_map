function indv_analysis(subjID, isSaveFig)
    % quick overview of the data, combine all the data from same session
    % and get pc and rt results depends on the conditions
    % usage: ex. indv_analysis('OO', 1) 
    %      : Analyze subject 'OO's data and save the bar plot for each session
    
    % get all tested conditions
    data_path = ['..' filesep 'Data' filesep upper(subjID)];
    f_list    = dir([data_path filesep upper(subjID) '*.dv']);
    ses_nm    = {};
    % get unique session name
    ses_unm = {};
    for ii=1:length(f_list)        
        tempSplit = strsplit(f_list(ii).name,'_');
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
        nEmptyBtw =  unique(ses_dt.nemptybetween);  % condition types
        anal_pc = [];
        anal_rt = [];
        for kk = 1:length(nEmptyBtw)
            for jj = 1:length(conds)
                anal_pc(kk,jj) = mean(ses_dt.iscorr(ses_dt.condition == conds(jj)...
                                                    & ses_dt.keyidx~=0 ... 
                                                    & ses_dt.nemptybetween == nEmptyBtw(kk)));
                anal_rt(kk,jj) = mean(ses_dt.rt(ses_dt.condition == conds(jj)...
                                                & ses_dt.keyidx~=0 ...
                                                & ses_dt.nemptybetween == nEmptyBtw(kk)));
            end
        end
        plot_pc = anal_pc;
        plot_rt = anal_rt;
        diff_pc = [];
        diff_rt = [];
        % compute difference from baseline 
        if length(conds)==3
            for kk = 1:length(nEmptyBtw)
                diff_pc(kk,:)=anal_pc(kk,1:2) - anal_pc(kk,3);
                diff_rt(kk,:)=anal_rt(kk,1:2) - anal_rt(kk,3);
            end
            plot_pc = diff_pc;
            plot_rt = diff_rt;
        end
        
        rr = floor(sqrt(length(nEmptyBtw)));
        cc = ceil(length(nEmptyBtw)/rr);
        fig = figure('position',[100 100 rr*length(conds)*250 cc*400]);
        
        for kk=1:length(nEmptyBtw)
            subplot(rr,cc*2,2*kk-1);
            bar(plot_pc(kk,:),0.8)
            text([1:length(plot_pc(kk,:))],plot_pc(kk,:)+0.05,strsplit(num2str(plot_pc(kk,:))),'HorizontalAlignment', 'center')
            if length(conds)~=3
                title([num2str(nEmptyBtw(kk)) ' empty between: PC'])
                ylabel('percent correct')
                ylim([0,1])
                xticklabels(condNms(conds))
            else
                title([num2str(nEmptyBtw(kk)) ' empty between: PC, Baseline: ' num2str(anal_pc(kk,3))])
                ylabel('percent correct difference')
                ylim([-0.2,0.2])
                xticklabels(condNms(conds(1:2)))
            end            
            xtickangle(30)
            subplot(rr,cc*2,2*kk);
            bar(plot_rt(kk,:),0.8)
            text([1:length(plot_rt(kk,:))],plot_rt(kk,:)+0.2,strsplit(num2str(plot_rt(kk,:))),'HorizontalAlignment', 'center')
            if length(conds)~=3
                title(['RT'])
                ylabel('reaction time (s)')
                ylim([0,3])
                xticklabels(condNms(conds))
            else
                title(['RT, Baseline: ' num2str(anal_rt(kk,3))])
                ylabel('reaction time (s) difference')
                ylim([-0.5,0.5])
                xticklabels(condNms(conds(1:2)))
            end
            
            xtickangle(30)
        end 
        
        if isSaveFig
            saveas(fig, [data_path filesep ses_unm{ii} '.png'])
            writetable(ses_dt,[data_path filesep upper(subjID) '_' ses_unm{ii} '.csv'])
        end
    end % end of session_analysis
    
end % end of indv_analysis()