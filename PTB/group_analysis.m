% Group analysis for pilot experiment

data_path = ['..' filesep 'Data' filesep 'current' filesep];
f_list    = dir([data_path filesep '*.csv']);
temp = readtable([data_path filesep f_list(1).name]);
nemptybetween = unique(temp.nemptybetween);
betweenNms = {'2disk between','4disk between','6disk between','8disk between'};
conds = unique(temp.condition);
condNms= {'Suppression', 'Enhancement', 'Baseline'};

cond_pe = [];
cond_rt = [];
cong_pe = [];
cong_rt = [];
for nn=1:length(f_list)
    indiv = readtable([data_path filesep f_list(nn).name]);
    for ff = 1:length(nemptybetween)
        ind_f = (indiv.nemptybetween == nemptybetween(ff)) & indiv.keyidx~=0;
        % by condition
        for cd = 1:length(conds)
            ind_cd = ind_f & indiv.condition==conds(cd);
            cond_pe(nn,ff,cd) = 1-mean(indiv.iscorr(ind_cd));
            cond_rt(nn,ff,cd) = mean(indiv.rt(ind_cd));
        end % end of condition loop
        cond_pe(nn,ff,:) = cond_pe(nn,ff,:)-cond_pe(nn,ff,3);
        cond_rt(nn,ff,:) = cond_rt(nn,ff,:)-cond_rt(nn,ff,3);
       
        % by congruency (match & congruent)
        icond = 1;
        for m = -1:2:1
            ind_m = ind_f & (indiv.cue_target_tilt .* indiv.cue_distractor_tilt)==m;
            for c = -1:2:1
                ind = ind_m & (indiv.cue_target_tilt .* indiv.prob_target_tilt)==c;
                cong_pe(nn,ff,icond)= 1-mean(indiv.iscorr(ind));
                cong_rt(nn,ff,icond)= mean(indiv.rt(ind));
                icond = icond+1;
            end % end of task: cue-probe target congruent loop
        end % end of cue: target-distractor match loop
    end % end of nemptybetween loop    
end % end of subject loop

cond_pe(:,:,3) = [];
cond_rt(:,:,3) = [];
mean_cond_pe = squeeze(mean(cond_pe,1));
mean_cond_rt = squeeze(mean(cond_rt,1));

% plots 
% per condition plot
fig1 = figure('position',[100 100 3*300 2*400]);
subplot(2,1,1)
hold on
errorbar((repmat([1:length(nemptybetween)],2,1)+[-0.14;0.14])',mean_cond_pe,squeeze(std(cond_pe,1)),'.k')
b=bar(mean_cond_pe, 0.8);
text([1:length(nemptybetween)]-0.2,mean_cond_pe(:,1)+0.01,strsplit(num2str(mean_cond_pe(:,1)')),'HorizontalAlignment', 'center')
text([1:length(nemptybetween)]+0.2,mean_cond_pe(:,2)-0.01,strsplit(num2str(mean_cond_pe(:,2)')),'HorizontalAlignment', 'center')
xticks([1,2,3,4])
xticklabels(betweenNms)
ylim([-0.1 0.1])
ylabel('Percent Error difference')
legend(b,condNms{conds})
subplot(2,1,2)
hold on
errorbar((repmat([1:length(nemptybetween)],2,1)+[-0.14;0.14])',mean_cond_rt,squeeze(std(cond_rt,1)),'.k')
b=bar(mean_cond_rt, 0.8)
text([1:length(nemptybetween)]-0.2,mean_cond_rt(:,1)+0.01,strsplit(num2str(mean_cond_rt(:,1)')),'HorizontalAlignment', 'center')
text([1:length(nemptybetween)]+0.2,mean_cond_rt(:,2)-0.01,strsplit(num2str(mean_cond_rt(:,2)')),'HorizontalAlignment', 'center')
xticks([1,2,3,4])
xticklabels(betweenNms)
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
xticklabels(betweenNms)
xtickangle(30)
ylim([0 0.1])
ylabel('Percent Error')
title('Cue-probe target tilt in the same orientation')
subplot(2,3,2)
bar(mean_cong_pe(:,[2,4]), 0.8)
xticks([1,2,3,4])
xticklabels(betweenNms)
xtickangle(30)
ylim([0 0.1])
title('Different orientation')
subplot(2,3,3)
bar([mean(mean_cong_pe(:,1:2),2),mean(mean_cong_pe(:,3:4),2)], 0.8)
xticks([1,2,3,4])
xticklabels(betweenNms)
xtickangle(30)
ylim([0 0.1])
ylim([0 0.1])
title('Combined')
legend({'target-disctractor: congruent','target-disctractor: incongruent'})
subplot(2,3,4)
bar(mean_cong_rt(:,[1,3]), 0.8)
xticks([1,2,3,4])
xticklabels(betweenNms)
xtickangle(30)
ylabel('reaction time')
ylim([0.5 0.9])
subplot(2,3,5)
bar(mean_cong_rt(:,[2,4]), 0.8)
xticks([1,2,3,4])
xticklabels(betweenNms)
xtickangle(30)
ylim([0.5 0.9])
subplot(2,3,6)
bar([mean(mean_cong_rt(:,1:2),2),mean(mean_cong_rt(:,3:4),2)], 0.8)
xticks([1,2,3,4])
xticklabels(betweenNms)
xtickangle(30)
ylim([0.5 0.9])
saveas(fig2, [data_path filesep 'byCongruency.png'])