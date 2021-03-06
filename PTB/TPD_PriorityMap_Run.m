function TPD_PriorityMap_Run
	global param
    rng Shuffle;
    
    param.numCal = 10;
    param.OBS_DIST = 65;                  % [cm] % Martin's setting
    param.stimElev = 0;                   % [arcmin] (-) upward, (+) downward
    param.backgrColor = [64, 64, 64];     % 8 cd/m2 background color
    param.respKeys = {'up','down'};       % keys to listen for responses 
    param.isIncorrFeedb = false;          % in correct feedback
    param.isISFix = logical(1);           % fixation dot in isi frame
    param.stimDursToBeUsedFrms = [4];     % [frames] sci.hz = 60Hz % default  4:  67ms
    param.isiDursToBeUsedFrms  = [11];    % [frames] sci.hz = 60Hz % default 11: 183ms
    
    %Eyetracker settings
    param.etDevice = 'none';          % Which eyetracker shall be used 'TET', 'SMI', 'GP3', 'eyelink', 'none'
    param.etOnline = true;
    
    % GUI: get user inputs
    param.subjID = char(upper(inputdlg({'Please enter the subject ID:'}, 'SubjID')));
    
    condNm = {'main_3disk_flicker','main_3disk_nonret', 'main_3disk_ret',...
         'train_3disk_flicker','train_3disk_nonret','train_3disk_ret',... 
         'training', 'gif_making'};
    condSelect = listdlg('PromptString', 'Which exptype would you like to run?', ...
                         'SelectionMode', 'single', 'ListString', ...
                        {condNm{:}});
    if isempty(condSelect) %if user canceled dialog window
	    return             %end program silently
	end 
	param.cond = condNm{1,condSelect};	
    param.condNmPostfix = '_None';     % to specify data filename % should be cleared in this script
    
    param.nx144 = 3; % total num of trials will be param.nx144 * 144 (default: 3) % should be cleared in this script
    
    switch param.cond
        case 'main_3disk_flicker' 
            param.nDisks = 3;
            param.nEmptyBetween = [2];   % the num of empty disks in between cue-probe stim (should be even number)
            param.isLeftAlways  = 1; % if the disks are always presented on the left (if 1 stationally 3 disk condition / flickerling condition)
            param.isNonRet      = 1; % the non-ret (1) or ret (0) condition  (transferred to vector in initialize_params())
            param.condNmPostfix = '_3diskflicker'; % to specify data filename % should be cleared in this script
        case 'main_3disk_nonret' 
            param.nDisks = 3;
            param.nEmptyBetween = [2];   % the num of empty disks in between cue-probe stim (should be even number)
            param.isLeftAlways  = 0;
            param.isNonRet      = 1; % the non-ret (1) or ret (0) condition  (transferred to vector in initialize_params())
            param.condNmPostfix = '_3disknonret'; % to specify data filename % should be cleared in this script
        case 'main_3disk_ret' 
            param.nDisks = 3;
            param.nEmptyBetween = [2];   % the num of empty disks in between cue-probe stim (should be even number)
            param.isLeftAlways  = 0;
            param.isNonRet      = 0; % the non-ret (1) or ret (0) condition  (transferred to vector in initialize_params())
            param.condNmPostfix = '_3diskret'; % to specify data filename % should be cleared in this script        
       
        case 'train_3disk_flicker'
            param.nDisks = 3;
            param.nEmptyBetween = [2];   % the num of empty disks in between cue-probe stim (should be even number)
            param.isLeftAlways  = 1; % if the disks are always presented on the left (if 1 stationally 3 disk condition / flickerling condition)
            param.isNonRet      = 1; % the non-ret (1) or ret (0) condition  (transferred to vector in initialize_params())
            init_train_params()
        case 'train_3disk_nonret'
            param.nDisks = 3;
            param.nEmptyBetween = [2];   % the num of empty disks in between cue-probe stim (should be even number)
            param.isLeftAlways  = 0; 
            param.isNonRet      = 1; % the non-ret (1) or ret (0) condition  (transferred to vector in initialize_params())
            init_train_params()
        case 'train_3disk_ret'
            param.nDisks = 3;
            param.nEmptyBetween = [2];   % the num of empty disks in between cue-probe stim (should be even number)
            param.isLeftAlways  = 0; 
            param.isNonRet      = 0; % the non-ret (1) or ret (0) condition  (transferred to vector in initialize_params())
            init_train_params()
        
        case 'training'
            param.numCal = 10;
            param.nDisks = 3;
            param.stimDursToBeUsedFrms = [12];    % [frames] sci.hz = 60Hz 
            param.isiDursToBeUsedFrms  = [12];    % [frames] sci.hz = 60Hz 
            param.isNonRet       = zeros(param.numCal,1);
            param.condition      = 2*ones(param.numCal,1);           % 1. Suppression(target position == d), 2. Enhancement(target position == p), 3. Baseline(target position == ~d & ~p) conditions
            param.distractor_pos = 1*ones(param.numCal,1);           % Cue frame's Distractor position 1/3
            param.cue_target_pos = 2*ones(param.numCal,1);           % Cue frame's target position 1/2 (except distractor position. It's relative position compared to distractor position. So, the exact position values are specified in the trial sessions.
            param.prob_target_pos= param.cue_target_pos;        % prob-target position: depends on the conditions
            param.cue_target_tilt    = (round(rand(param.numCal,1))-0.5)*2;    % Cue frame's target tilt orientation -1: left, 1: right
            param.cue_distractor_tilt= (round(rand(param.numCal,1))-0.5)*2;% Cue frame's distractor tilt orientation -1: left, 1: right
            param.prob_target_tilt   = (round(rand(param.numCal,1))-0.5)*2;   % Probe frame's target tilt orientation -1: left, 1: right
        case 'gif_making'
            param.numCal = 5;
            param.nDisks = 3;
            param.stimDursToBeUsedFrms = [12];    % [frames] sci.hz = 60Hz 
            param.isiDursToBeUsedFrms  = [12];    % [frames] sci.hz = 60Hz
            param.isISFix = logical(1);  % fixation dot in isi frame
            param.isLeftAlways  = 1; % if the disks are always presented on the left (if 1 stationally 3 disk condition / flickerling condition)
            param.isNonRet       = ones(param.numCal,1);
            param.condition      = ones(param.numCal,1);           % 1. Suppression(target position == d), 2. Enhancement(target position == p), 3. Baseline(target position == ~d & ~p) conditions
            param.distractor_pos = randi(2,param.numCal,1,1);           % Cue frame's Distractor position 1/3
            param.cue_target_pos = 3*ones(param.numCal,1);              % Cue frame's target position 1/2 (except distractor position. It's relative position compared to distractor position. So, the exact position values are specified in the trial sessions.
            param.prob_target_pos= randi(3,param.numCal,1,1); % prob-target position: depends on the conditions
            param.cue_target_tilt    = (round(rand(param.numCal,1))-0.5)*2;    % Cue frame's target tilt orientation -1: left, 1: right
            param.cue_distractor_tilt= (round(rand(param.numCal,1))-0.5)*2;% Cue frame's distractor tilt orientation -1: left, 1: right
            param.prob_target_tilt   = (round(rand(param.numCal,1))-0.5)*2;   % Probe frame's target tilt orientation -1: left, 1: right
    end    
    
    initialize_params();
    
    if param.numCal == 0
        sprintf(['All the trials in ' param.cond 'have been tested'])
        return;
    end
    
    % Call stimulus program
    run = struct2pvcell(param);
    % Use Root Stim Script
    if ~strcmp(param.cond, 'gif_making')
        TPD_PriorityMap_Stim(run{:});
    else
        TPD_PriorityMap_gifmaking_Stim(run{:});
    end
end % end of TPD_PriorityMap_Run()

function init_train_params()
    global param
    param.numCal = 18;
    param.isIncorrFeedb = true;
    param.isNonRet           = ones(param.numCal,1)*param.isNonRet;
    param.cue_target_tilt    = (round(rand(param.numCal,1))-0.5)*2;    % Cue frame's target tilt orientation -1: left, 1: right
    param.cue_distractor_tilt= (round(rand(param.numCal,1))-0.5)*2;% Cue frame's distractor tilt orientation -1: left, 1: right
    param.prob_target_tilt   = (round(rand(param.numCal,1))-0.5)*2;   % Probe frame's target tilt orientation -1: left, 1: right
        
    param.condition      = randi(3,param.numCal,1);     % 1. Suppression(target position == d), 2. Enhancement(target position == p), 3. Baseline(target position == ~d & ~p) conditions
    bar_comb = [1,2,3;1,3,2;2,1,3;2,3,1;3,1,2;3,2,1];   % combinations of 3 positions
    bar_positions    = randi(6,param.numCal,1);         % conbination type: 1-6
    bar_positions    = bar_comb(bar_positions,:);       % trial by trial bar positions
    param.distractor_pos = bar_positions(:,1);              % cue-distractor position
    param.cue_target_pos = bar_positions(:,2);              % cue-target position
    subind = sub2ind(size(bar_positions), [1:param.numCal], param.condition');
    param.prob_target_pos= bar_positions(subind');  % prob-target position: depends on the conditions

end % end of init_train_params()

function initialize_params
    global param
    % set-up file-location and the filenames
    [ProgDir,Fl,Ext] = fileparts(which(mfilename)); % get program directory
    dataPath = [ProgDir filesep '..' filesep 'Data' filesep param.subjID]; % path to data folder
    flieList = dir([ProgDir filesep '..' filesep 'Data']); % check if there is a already made subjID in Data folder (manually check same ID for diff subjects)
    % Check RunNum of current condition
    if ~any(strcmp({flieList.name}, param.subjID))
        runNum = 1; % count num of runs
        mkdir(dataPath)
    end
	fileListSub = dir(dataPath);    
    
    if any(strfind(lower(param.cond),lower('main')))
        nTot = 144*param.nx144; % total trials should be factor of 144
        param.numCal = 108;
        % Check num of succeeded main sessions and set current session num
        succeedNum = numel(strfind(strcat(fileListSub.name),['sessionSucceed_' param.cond param.condNmPostfix])); 

        if succeedNum == 0 % if first main session: 14 configurations
            %counterbalance condition(3), bar_positions (6), cue_target_tilt(2), cue_distractor_tilt(2), probe_target_tilt(2) conditions
            a = fullfact([3 6 2 2 2]);                      % condition(3), bar_positions (6),cue_target_tilt(2), cue_distractor_tilt(2), probe_target_tilt(2)
            b = repmat(a, nTot/size(a,1),1);             % replicate the matrix until it has as many rows as there are trials 
            idx = Shuffle([1:size(b,1)]);                   % create a random index for the rows in b
            c = b(idx,:);                                   % reorder rows in b according to the random index
            c(c(:,3) == 2,3) = -1;                          % In columns 3:4, replace the 2s by -1s so we only have -1s and 1s in c(:,3:4).  
            c(c(:,4) == 2,4) = -1;                          % In columns 3:4, replace the 2s by -1s so we only have -1s and 1s in c(:,3:4).  
            c(c(:,5) == 2,5) = -1;                          % In column 5, replace the 2s by -1s so we only have -1s and 1s in c(:,5).  
            param.condition          = c(:,1);              % 1. Suppression(target position == d), 2. Enhancement(target position == p), 3. Baseline(target position == ~d & ~p) conditions
            param.cue_target_tilt    = c(:,3);              % Cue frame's target tilt orientation -1: left, 1: right
            param.cue_distractor_tilt= c(:,4);              % Cue frame's distractor tilt orientation -1: left, 1: right
            param.prob_target_tilt   = c(:,5);              % Probe frame's target tilt orientation -1: left, 1: right
            param.isNonRet           = ones(nTot,1)*param.isNonRet;  % the non-ret (1) or ret (0) condition
            param.distractor_pos     = zeros(nTot,1);       % Cue frame's Distractor position 1/3
            param.cue_target_pos     = zeros(nTot,1);       % Cue frame's target position 1/2 (except distractor position. It's relative position compared to distractor position. So, the exact position values are specified in the trial sessions.
            param.prob_target_pos    = zeros(nTot,1);       % Prob frame's target position, it depends on condition and cue frame's distractor & target positions. 

            % set-up positions relativey
            bar_comb = [1,2,3;1,3,2;2,1,3;2,3,1;3,1,2;3,2,1];   % combinations of 3 positions
            bar_positions    = c(:,2);                          % conbination type: 1-6
            bar_positions    = bar_comb(bar_positions,:);       % trial by trial bar positions
            param.distractor_pos = bar_positions(:,1);              % cue-distractor position
            param.cue_target_pos = bar_positions(:,2);              % cue-target position
            subind = sub2ind(size(bar_positions), [1:nTot], param.condition');
            param.prob_target_pos= bar_positions(subind');  % prob-target position: depends on the conditions

            init_succeedTxt(dataPath);
            clear a b c idx
        elseif succeedNum % if sessionSucceedNum file exists for current condition
            readNset_params(dataPath, succeedNum);            
        end
        succeedNum = succeedNum+1;        
        param.saveTxtname = [dataPath filesep param.subjID '_sessionSucceed_' param.cond param.condNmPostfix '_SUCCEED_' num2str(succeedNum) '.txt' ];
    end 
 
    
	% set filenames 
    runNum = length(strfind(strcat(fileListSub.name),[param.cond, param.condNmPostfix, '_RUN_']))+1; % count num of runs
    
    param.saveFilename = [dataPath filesep param.subjID '_TPDpMap_' param.cond param.condNmPostfix '_RUN_' num2str(runNum) '.dv'];
    
    % remove unnecessary fields 
    param = rmfield(param, {'condNmPostfix', 'nx144'});
    
end % end of initialize_params()

function init_succeedTxt(dataPath)
    global param 
    
    SavedFields = {'condition';'isNonRet';'cue_target_tilt';'cue_distractor_tilt';...
                  'prob_target_tilt';'distractor_pos';'cue_target_pos';'prob_target_pos'};
    % open sessionSucceed.txt file: import txt filename from RUNfile (pp.savetxtname)
    txtid = fopen([dataPath filesep param.subjID '_init_' param.cond param.condNmPostfix '_SucceedTxt.txt'], 'w');
    
    fprintf(txtid, 'Date: %d-%d-%d %d:%d:%4.2f\r\n', clock);
    fprintf(txtid, 'SubjID: %s\r\n', param.subjID);
    fprintf(txtid, 'CondNm: %s\r\n', param.cond);
    fprintf(txtid, 'SavedFields: \n'); 
    fprintf(txtid, '%s;', SavedFields{:});
    fprintf(txtid, '\n');
    for ii = 1:length(SavedFields)
        fprintf(txtid, '%s: \n', SavedFields{ii});
        fprintf(txtid, '%d;', param.(SavedFields{ii}));
        fprintf(txtid, '\n');
    end
    
    fclose(txtid);
    
    for ii = 1:length(SavedFields)
        param.(SavedFields{ii}) = param.(SavedFields{ii})(1:param.numCal);
    end
    param.startTrialID = 0;
    
    ismade = 1; 
    if ismade; sprintf('%s is made successfully', [dataPath filesep param.subjID '_init_' param.cond param.condNmPostfix '_SucceedTxt.txt']);end
end % end of init_succeedTxt()

function readNset_params(dataPath, succeedNum)
    global param
    
    % Find the last tested ID: finTrialID
    sfname = [dataPath filesep param.subjID '_sessionSucceed_' param.cond param.condNmPostfix '_SUCCEED_' num2str(succeedNum) '.txt' ];
    txtid = fopen(sfname, 'r');
    text = fgetl(txtid);
    while ~any(strfind(lower(text),lower('finTrialID')))
        text = fgetl(txtid);
    end
    text = fgetl(txtid);
    startTrialID = strsplit(text,';');
    startTrialID = str2num(startTrialID{1});    
    
    % load pre-registered randomized parameters
    initfn = [dataPath filesep param.subjID '_init_' param.cond param.condNmPostfix '_SucceedTxt.txt'];
    txtid = fopen(initfn, 'r');
    text = fgetl(txtid);
    while ~any(strfind(lower(text),lower('SavedFields')))
        text = fgetl(txtid);
    end
    text = fgetl(txtid);
    SavedFields = strsplit(text,';');
    for ii = 1:length(SavedFields)-1
        text = fgetl(txtid);
        disp(text)
        text = fgetl(txtid);
        disp(text)
        splTxt = strsplit(text,';');
        % Stop the experiment if all the trials are all tested
        if ii ==1
            sprintf([num2str(startTrialID) ' trials among ' num2str(length(splTxt)-1) 'in condition' param.cond 'have been tested'])
            if startTrialID+1 >= length(splTxt)
                param.numCal = 0;
                sprintf(['All the trials in ' param.cond 'have been tested'])
                return;
            end
        end
        param.(SavedFields{ii}) = cellfun(@str2num,splTxt(startTrialID+1:min(startTrialID+param.numCal, end-1)))';
    end
    fclose(txtid);
    
    param.numCal = length(param.(SavedFields{1}));
    param.startTrialID = startTrialID;
        
end % end of readNset_params()
   