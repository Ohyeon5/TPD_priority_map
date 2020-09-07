function group_result = group_analysis(isSaveIndvAnal, isPlot)

    data_path = ['..' filesep 'Data' filesep 'current' filesep];
    subjIDs = readtable([pwd filesep '..' filesep 'logbook.xlsx']); 
    subjIDs = subjIDs.subjID;
    subjIDs = subjIDs(~cellfun(@isempty, subjIDs));

    % remove excluded participants
    foldList = dir(data_path);
    incIdx   = [];
    for ii=1:numel(foldList)
        incIdx = [incIdx; find(contains(subjIDs,foldList(ii).name))];
    end
    incIdx = unique(incIdx);
    subjIDs = subjIDs(incIdx);
        
    group_result = [];
    
    for subjNum = 1:length(subjIDs) % individual analysis
        indv_result = indv_analysis(subjIDs{subjNum}, isSaveIndvAnal);
        for sess = 1:length(indv_result)
            sess_name = indv_result(sess).sessNm;  % session name
            conds     = unique(indv_result(sess).table.condition);
            
            % per condition: 1)Suppression 2)Enhancement 3)Baseline
            for cd = 1:length(conds)
                ind_cd = indv_result(sess).table.condition==conds(cd) & indv_result(sess).table.keyidx~=0;
                group_result.(sess_name).cond.pc(subjNum,cd) = 1-mean(indv_result(sess).table.iscorr(ind_cd));
                group_result.(sess_name).cond.rt(subjNum,cd) = mean(indv_result(sess).table.rt(ind_cd));
            end
            
            % per congruency
            % 1) diff orient & cue target-disctractor incong 2) same orient & cue target-disctractor incong  
            % 3) diff orient & cue target-disctractor cong   4) same orient & cue target-disctractor cong  
            icong = 1;
            for m = -1:2:1
                for c = -1:2:1
                    ind = (indv_result(sess).table.cue_target_tilt .* indv_result(sess).table.prob_target_tilt)==c & indv_result(sess).table.keyidx~=0;
                    group_result.(sess_name).cong.pc(subjNum,icong) = 1-mean(indv_result(sess).table.iscorr(ind));
                    group_result.(sess_name).cong.rt(subjNum,icong) = mean(indv_result(sess).table.rt(ind));
                    icong = icong+1;
                end % end of task: cue-probe target congruent loop
            end % end of cue: target-distractor match loop
        end % end of sess loop        
    end % end of subjIDs loop
    
    % plots 
    if isPlot
        sess_names = fieldnames(group_result);
        % initialize parameters
        condNms= {'Suppression', 'Enhancement', 'Baseline'};
        nSubj  = length(subjIDs);
        nSess  = length(sess_names);
        nConds = length(condNms);
        nCongs = size(group_result.(sess_names{1}).cong.pc,2);
        cond_pe = zeros(nSubj,nSess,nConds);
        cond_rt = zeros(nSubj,nSess,nConds);
        cong_pe = zeros(nSubj,nSess,nCongs);
        cong_rt = zeros(nSubj,nSess,nCongs);    
        
        for sess = 1: nSess
            sessNm  = sess_names{sess};
                        
            % per condition plot
            cond_pe(:,sess,:) = reshape(group_result.(sessNm).cond.pc(:,:),nSubj,1,nConds);
            cond_rt(:,sess,:) = reshape(group_result.(sessNm).cond.rt(:,:),nSubj,1,nConds);
            cong_pe(:,sess,:) = reshape(group_result.(sessNm).cong.pc(:,:),nSubj,1,nCongs);
            cong_rt(:,sess,:) = reshape(group_result.(sessNm).cong.rt(:,:),nSubj,1,nCongs);
            
            cond_pe(:,sess,:) = cond_pe(:,sess,:) - cond_pe(:,sess,3);
            cond_rt(:,sess,:) = cond_rt(:,sess,:) - cond_rt(:,sess,3);
        end % end to collect values to plot
        
        cond_pe(:,:,3) = [];
        cond_rt(:,:,3) = [];
        mean_cond_pe = squeeze(mean(cond_pe,1));
        mean_cond_rt = squeeze(mean(cond_rt,1));

        fig1 = figure('position',[100 100 3*300 2*400]);
        subplot(2,1,1)
        hold on
        errorbar((repmat([1:nSess],2,1)+[-0.14;0.14])',mean_cond_pe,squeeze(std(cond_pe,1)),'.k')
        b=bar(mean_cond_pe, 0.8);
        text([1:nSess]-0.2,mean_cond_pe(:,1)+0.01,strsplit(num2str(mean_cond_pe(:,1)')),'HorizontalAlignment', 'center')
        text([1:nSess]+0.2,mean_cond_pe(:,2)-0.01,strsplit(num2str(mean_cond_pe(:,2)')),'HorizontalAlignment', 'center')
        xticks([1:nSess])
        xticklabels(sess_names)
        ylim([-0.1 0.1])
        ylabel('Percent Error difference')
        legend(b,condNms{conds})
        subplot(2,1,2)
        hold on
        errorbar((repmat([1:nSess],2,1)+[-0.14;0.14])',mean_cond_rt,squeeze(std(cond_rt,1)),'.k')
        b=bar(mean_cond_rt, 0.8)
        text([1:nSess]-0.2,mean_cond_rt(:,1)+0.01,strsplit(num2str(mean_cond_rt(:,1)')),'HorizontalAlignment', 'center')
        text([1:nSess]+0.2,mean_cond_rt(:,2)-0.01,strsplit(num2str(mean_cond_rt(:,2)')),'HorizontalAlignment', 'center')
        xticks([1:nSess])
        xticklabels(sess_names)
        ylim([-0.3 0.3])
        ylabel('Reaction time difference')
        legend(b,condNms{conds})
        saveas(fig1, [data_path filesep 'byCondition.png'])

        % per congruency plot            
        mean_cong_pe = squeeze(mean(cong_pe,1));
        mean_cong_rt = squeeze(mean(cong_rt,1));
        fig2 = figure('position',[100 100 3*300 2*400]);
        subplot(2,3,1)
        bar(mean_cong_pe(:,[1,3]), 0.8)
        xticks([1,2,3,4])
        xticklabels(sess_names)
        xtickangle(30)
        ylim([0 0.5])
        ylabel('Percent Error')
        title('Cue-probe target tilt in the same orientation')
        subplot(2,3,2)
        bar(mean_cong_pe(:,[2,4]), 0.8)
        xticks([1,2,3,4])
        xticklabels(sess_names)
        xtickangle(30)
        ylim([0 0.5])
        title('Different orientation')
        subplot(2,3,3)
        bar([mean(mean_cong_pe(:,1:2),2),mean(mean_cong_pe(:,3:4),2)], 0.8)
        xticks([1,2,3,4])
        xticklabels(sess_names)
        xtickangle(30)
        ylim([0 0.5])
        title('Combined')
        legend({'target-disctractor: congruent','target-disctractor: incongruent'})
        subplot(2,3,4)
        bar(mean_cong_rt(:,[1,3]), 0.8)
        xticks([1,2,3,4])
        xticklabels(sess_names)
        xtickangle(30)
        ylabel('reaction time')
        ylim([0.5 1.5])
        subplot(2,3,5)
        bar(mean_cong_rt(:,[2,4]), 0.8)
        xticks([1,2,3,4])
        xticklabels(sess_names)
        xtickangle(30)
        ylim([0.5 1.5])
        subplot(2,3,6)
        bar([mean(mean_cong_rt(:,1:2),2),mean(mean_cong_rt(:,3:4),2)], 0.8)
        xticks([1,2,3,4])
        xticklabels(sess_names)
        xtickangle(30)
        ylim([0.5 1.5])
        saveas(fig2, [data_path filesep 'byCongruency.png'])
    end % end of plot functions
end % end of the main function