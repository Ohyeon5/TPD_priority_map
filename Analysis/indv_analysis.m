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
    
    fig = fc_plot_pe_rt(anal_pe, anal_rt, ses_unm_lab, condNms, conds, subjID);
            
    if isSaveAnal
        saveas(fig, [data_path filesep subjID '_ind_anal.png'])
    end

    fprintf(['Done: ' subjID])
end % end of indv_analysis()