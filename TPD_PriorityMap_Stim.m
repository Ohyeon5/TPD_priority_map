function [myData] = TPD_PriorityMap_Stim(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Script is adapted from Oh-hyeon Choung's TPD_priority_map() ver 5.0%
%  										   BY Oh-hyeon Choung (2020-03-09)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% History 
% Ver 1.0 (2020-03-09)
% 	- First release 
% 	- Mainly cleaned TPD_priority_map() Ver 5.0 in Stim-Run form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tBlockStart = tic;
rng('shuffle')  %Shuffle seed of random number generator
global sci pp d tt

%% Input parameters 
pp = inputParser;

pp.addParamValue('subjID', 'OOO', @ischar);                 % subject IDs, default='OOO': all debuging files
pp.addParamValue('TestMode', 0, @isscalar);                 % 0: Experiment, 1: Debug
pp.addParamValue('saveFilename', 'test', @ischar);          % saving file's filename ex) 'OOO_GroupingVernier_cond1_RUN_1.dv'
pp.addParamValue('saveTxtname', 'test', @ischar);           % save sessionSucceed.txt filename and it's path
pp.addParamValue('cond', 'Demo', @ischar);                  % condition name
pp.addParamValue('task', 'priority_map', @ischar);          % task type: 'priority_map' in the moment, more will be included
pp.addParamValue('nDisks', 3, @isscalar);                   % the number of disks 

pp.addParamValue('isDebug', 0, @islogical);                 % Is debug mode?
pp.addParamValue('isSaveData', 1, @islogical);              % Is save data?
pp.addParamValue('isTimeoutFeedb', 0, @islogical);          % Is time out Feedback?
pp.addParamValue('isIncorrFeedb', 0, @islogical);           % Is in-corrrect response feedback?
pp.addParamValue('isMask', 0, @islogical);                  % Is mask?
pp.addParamValue('isISFix', 1, @islogical);                 % Is interstimuli fixation?
pp.addParamValue('startTrialID', 0, @isscalar);             % the last tested trial ID (which trial to start for this exp),saved in 'saveTxtname'

pp.addParamValue('OBS_DIST', 66, @isscalar);                % [cm] % observer distance
pp.addParamValue('numCal', 72, @isscalar);                  % [trials] number of trials per block

%Condition specific parameters 
pp.addParamValue('condition', [], @ismatrix);               % Experimental conditions: 1. Suppression(target position == d), 2. Enhancement(target position == p), 3. Baseline(target position == ~d & ~p) conditions
pp.addParamValue('cue_target_tilt', [], @ismatrix);         % Cue frame's target tilt orientation -1: left, 1: right
pp.addParamValue('cue_distractor_tilt', [], @ismatrix);     % Cue frame's distractor tilt orientation -1: left, 1: right
pp.addParamValue('prob_target_tilt', [], @ismatrix);        % Probe frame's target tilt orientation -1: left, 1: right
pp.addParamValue('isNonRet', [], @ismatrix);                % the non-ret (1) or ret (0) condition
pp.addParamValue('distractor_pos', [], @ismatrix);          % Cue frame's Distractor position 1/3
pp.addParamValue('cue_target_pos', [], @ismatrix);          % Cue frame's target position 1/2 (except distractor position. It's relative position compared to distractor position. So, the exact position values are specified in the trial sessions.
pp.addParamValue('prob_target_pos', [], @ismatrix);         % Prob frame's target position, it depends on condition and cue frame's distractor & target positions. 

%Repetition settings
pp.addParamValue('nStimPres', 2, @isscalar);                % number of stimulus presentations per trial %4
pp.addParamValue('nEmptyBefore', 1, @isscalar);             % number of disk presentations without dot before %1
pp.addParamValue('nEmptyBetween', [2], @ismatrix);            % number of disk presentations in-between two stimulus presentation %2 (should be an even number)
pp.addParamValue('nEmptyAfter', 0, @isscalar);              % number of disk presentations without dot after %1

%Temporal settings
pp.addParamValue('isiDursToBeUsedFrms', [12], @ismatrix);   % desired isi durations in number of screen refreshs 
pp.addParamValue('stimDursToBeUsedFrms', [12], @ismatrix);  % desired stimulus durations in number of screen refreshs 
pp.addParamValue('maskDursToBeUsedFrms', [2], @ismatrix);   % desired mask durations in number of screen refreshs 
pp.addParamValue('itiDurFrms', 60, @isscalar);              % iti (after stim break) duration in number of screen refreshs %60
pp.addParamValue('respTimeoutSecs', 3, @isscalar);          % bar width % 4
pp.addParamValue('freqHz', 60, @isscalar);                  % desired screen refresh rate %60

%Spatial settings
pp.addParamValue('fxRad', 1.5, @isscalar);                  % [arcmin] radius of fixationpoint
pp.addParamValue('diskRad', 90, @isscalar);                 % [arcmin] radius of each disk.
pp.addParamValue('stimShift', 0, @isscalar);                % [arcmin] frame shift (-) left, (+) right % for 1Disk conditions, stimShift should be specified
pp.addParamValue('stimElev', -180, @isscalar);              % [arcmin] stimulus elevation (fxpoint-center to stim-center, negative = up) 
pp.addParamValue('barWidth', 4.5, @isscalar);              	% [arcmin] width of the disk-enclosing frame %(approximately 3 pixel)

%Color settings
pp.addParamValue('backgrColor', [128, 128, 128], @ismatrix);%50cd/m2 % Clr.backgr = [0,0,0];           %0.4cd/m2 back ground should be black to make circle non-visible
pp.addParamValue('fxColor', [255, 0, 0], @ismatrix);        %20cd/m2
pp.addParamValue('diskColor', [0 0 0], @ismatrix);          %0.4cd/m2
pp.addParamValue('barColor', [255, 255, 255], @ismatrix);   %100cd/m2
pp.addParamValue('maskColor', [10,10,10], @ismatrix);       %???

% Eye tracker settings
pp.addParamValue('etDevice', 'none', @ischar);          % Which eyetracker to use: 'TET','SMI','GP3','none','mouse'
pp.addParamValue('etOnline', false, @isscalar);         % [bool] % whether online eyetracker fixation control (etCheck)

pp.parse(varargin{:});
pp = pp.Results;

%% Parameter shaping and ping possible errors
try
    if size(pp.condition,1)~=pp.numCal||size(pp.cue_target_tilt,1)~=pp.numCal||size(pp.cue_distractor_tilt,1)~=pp.numCal ...
       ||size(pp.prob_target_tilt,1)~=pp.numCal||size(pp.isNonRet,1)~=pp.numCal||size(pp.distractor_pos,1)~=pp.numCal...
       ||size(pp.cue_target_pos,1)~=pp.numCal||size(pp.prob_target_pos,1)~=pp.numCal
        error('Dimension mismatch in conditiion specific arguments.')
    end
    % some parameter related errors and settings 
    
    respKeys = {'left', 'right'}; % keys to listen for responses
    % Trial specific parameters 
    fdnms = {'condition';'cue_target_tilt';'cue_distractor_tilt';'prob_target_tilt';'isNonRet';...
             'distractor_pos';'cue_target_pos';'prob_target_pos';};
    for i = 1: numel(fdnms)
        tt.(fdnms{i}) = pp.(fdnms{i});
    end
    
    if length(pp.nEmptyBetween)==1
        pp.nEmptyBetween = repmat(pp.nEmptyBetween,pp.numCal,1);
    end
    
catch err
    commandwindow();
    rethrow(err);
end

%% 
try 
	%% Open screen and initialize screens
    if pp.TestMode
        lpsy.init('DontUseRT'); % relax reaction time criteria
    else
        lpsy.init();
    end
    
    sysInfo = lpsy.readSysInfo();
    
    % Screen info
    % Present stimuli on the main screen, where using multiple screens
    screens = Screen('Screens');
    screenNumber = 1;
    sci = Screen('Resolution',screenNumber,[],[],pp.freqHz); %get screen info
    sci.gamma = 1;             			% set to 1 to load linear gamma table using lpsy.setGammaTab ; %For use with lpsytb function: Scrn.gammaCalibCoeff = mean(str2num(sysInfo.colorcalibration.gamma)); % Desired gamma calibration coefficient 
	sci.gammaMode = 'grey';   			% if 'grey', lpsy.setGammaTab programs the gamma table with identical values for the R, G, and B channel
	sci.multisample = 0;       			% Desired anti-aliasing quality (brute force fullscreen method). Max value is hardware dependent, read 'help AntiAliasing'. MR suggests a max value of 8.

    if pp.TestMode
        wndScale = 0.4; % scale down to 40% of window
        sz = round(wndScale*[sci.width sci.height]);
        [sci.wnd, sci.Rect] = Screen('OpenWindow', 0, 128, [sci.width-sz(1) 0 sci.width-1 sz(2)], [],[],[],9); % last parameter: multisample, 9 -- antialiasing
        Screen('glScale', sci.wnd, wndScale, wndScale);
    else % open a regular fullscreen window
        [sci.wnd, sci.Rect] = Screen('OpenWindow',screenNumber, pp.backgrColor, [], [], [], [], sci.multisample);
        % Load monitor-specific, calibrated Colour LookUp Table (CLUT)        
        lpsy.setGammaTab(sci.gamma, sci.gammaMode);
    end
    
    % enable blending in order to make anti-aliasing work
    Screen('BlendFunction', sci.wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % drawing parameters
    sci.white = WhiteIndex(sci.wnd);    % pixel value for white
    sci.black = BlackIndex(sci.wnd);    % pixel value for black
    sci.red = [255 0 0];                % RGB value for red
    [sci.Xcenter, sci.Ycenter] = RectCenter(sci.Rect);  %Coordinates of the center of the window (RectCenter is a function from Psychtoolbox) 
    sci.vd_m = pp.OBS_DIST;             % subject's viewing distance as meter
    lpsy.setScreenDistCM(pp.OBS_DIST);  % set observation distanc
    disp(sci)
    
    % Stimuli spatial parameters set to pixel
    pp.fxRadPx    = ceil(lpsy.arcmin2pix(pp.fxRad));  % ceil to avoid pixelized spatial errors
    pp.diskRadPx   = ceil(lpsy.arcmin2pix(pp.diskRad)); % ceil to avoid pixelized spatial errors
    pp.stimShiftPx = ceil(lpsy.arcmin2pix(pp.stimShift));% ceil to avoid pixelized spatial errors
    pp.stimElevPx  = ceil(lpsy.arcmin2pix(pp.stimElev)); % ceil to avoid pixelized spatial errors
    pp.barWidthPx  = ceil(lpsy.arcmin2pix(pp.barWidth));% ceil to avoid pixelized spatial errors
    pp.intStimDistPx = ceil(2.5*pp.diskRadPx/2)*2;      % inter stimulus distance (center-to-center) =2*diskRadPx+21;
    if pp.stimShiftPx == 0 && pp.nDisks ==1
        pp.stimShiftPx = pp.intStimDistPx*0.5;          % automatically compensate and push the disk to the center
    end
    
    %Set position of fixationpoint
    fxCntr = [sci.Xcenter, sci.Ycenter] ;
    stimCntr = [fxCntr(1)+pp.stimShiftPx, fxCntr(2)+pp.stimElevPx];
    
    % timing parameters 
    tt.isiDurFrms = Shuffle(repmat(pp.isiDursToBeUsedFrms', pp.numCal/length(pp.isiDursToBeUsedFrms),1));    %randomize with equal probability (but no counterbalancing)                        
    tt.stimDurFrms = Shuffle(repmat(pp.stimDursToBeUsedFrms', pp.numCal/length(pp.stimDursToBeUsedFrms),1)); %randomize with equal probability (but no counterbalancing)                        
    tt.maskDurFrms = Shuffle(repmat(pp.maskDursToBeUsedFrms', pp.numCal/length(pp.maskDursToBeUsedFrms),1)); %randomize with equal probability (but no counterbalancing)                        
    tt.isiDurSecs = tt.isiDurFrms/pp.freqHz;
    pp.itiDurSecs = pp.itiDurFrms/pp.freqHz;
    tt.stimDurSecs = tt.stimDurFrms/pp.freqHz;
    tt.maskDurSecs = tt.maskDurFrms/pp.freqHz;
    tt.isLeftStart = ones(pp.numCal,1);  %Disks start always on the left (necessary to be able to instruct in the 2diskNonret conditions with only black disks)                            
    tt.nEmptyBetween = pp.nEmptyBetween;
    
    % Data saving parameters 
    d.isTimingErr = false(pp.numCal,1); %initialize d.isTimingErr. +1 for the flip clearing the screen after the last stimPres
    d.keyIdx = zeros(pp.numCal,1);
    d.isCorr = false(pp.numCal,1);
    d.rt = zeros(pp.numCal,1);
    d.tToFixSecs = zeros(pp.numCal,1); %for fixation control at trial start
    d.isEtOk = false(pp.numCal,1);    %for retrospective control of fixation during trial after the trial ended
    
    %Eyetracker settings
	etBaselineSecs = 0.2;       % Amount of data to be saved pre trialstart (can be used for post-hoc drift correction) 
	isOnlineFixCtrl = pp.etOnline;

	%Eyetracker settings for lpsy.etWaitForFixation at trialstart
	pp.maxNoiseDeg = [0.8];        % default: 0.8deg, independently for x and y direction [0.5]
	pp.maxDevDeg = [1];            % default: 1.5deg [0.25]
	pp.driftDelaySecs = [0.15];    % default: 0.15s (time for one saccade). Time to wait before drift correction is activated. 
	pp.tStartLookForFix = [-1];    % default: moment the function is called. GetSecs timestamp at which the function starts looking for fixations. Potentially useful to enforce an ITI period.
	pp.minFixDurSecs = 0.25;       % default: 0.25s.  Minimum duration of a central fixation before trial starts (don't choose this too long, or people will blink immediately after trialstart) %Function returns (and trials starts) if an interval of minFixDurSecs is free of blinks and at least >80% of samples are at the target position (within the other criteria, e.g., maxDevDeg and maxNoiseDeg)
	pp.maxWaitForFixSecs = [999];  % default: infinite. Function returns and script is aborted with an error if the function has not detected a valid fixation after maxWaitForFixSecs.

    %Get stim centers for use with MR's DrawRing.DrawRing function
    % Use of return variables:  
    % -diskCntrs(x/y, iDisk)  
    diskCntrs = getStimCntrs(stimCntr, pp.intStimDistPx);                          
        
    %dummy call to initialize the 'FillOval' function (first call is not actually drawn apparently)
    %This is a workaround for a bug in PTB    
    DrawRing.DrawRing(sci.wnd, fxCntr, pp.backgrColor, 2*pp.fxRadPx);
    
    % Connect to eyetracker 
    isEtConnected = lpsy.etConnect('Device', pp.etDevice, 'Recording', 'local', 'Eye', 'c'); 
    
    % Calibrate eyetracker, only for first trial of the experiment
    lpsy.etRunCalib('TargetCnt', 5, 'MaxWidth', 5, 'MaxHeight', 5, 'LoadOldCal', 'no');  
    
    %Hide the mouse cursor
    HideCursor;

%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- PTB PART --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%

    iTrial = 1;
    while iTrial <= length(tt.isLeftStart) 
        %% Setup the current trial        
        isEsc = 0; %reset isEsc 
        
        % Update each trial
        pp.nStimTot = pp.nStimPres + pp.nEmptyBefore + tt.nEmptyBetween(iTrial) + pp.nEmptyAfter; % Total number of stimulus Presentations
        isTimingErr = false(1,pp.nStimTot*2+1);
        
        % Get inner-stim positions for all disks and stimulus presentations of the trial (1=up, 2=right, 3=down, 4=left, 5=center, 6=none)
        innerStimPos = getInnerStimPos('bars', iTrial, stimCntr, pp.nStimPres);
        
        %draw fixationpoint 
        DrawRing.DrawRing(sci.wnd, fxCntr, pp.fxColor, 2*pp.fxRadPx);
        
        lpsy.releaseKeyWait();              %Wait until no keys or push buttons are pressed                
        lpsy.flip(sci.wnd);                                 
        WaitSecs(pp.itiDurSecs);   %enforce ITI (inter trial interval)
        
        %Wait until steady fixation on fixation point for at least
        %<minFixDurSecs> and perform drift correction
        tFix0 = GetSecs;        
        [isEsc, ~, err, ~, tStartOfFix] = lpsy.etWaitForFixation('MaxNoise', pp.maxNoiseDeg, 'MaxDev', pp.maxDevDeg, 'DriftDelay', pp.driftDelaySecs, 'StartTi', pp.tStartLookForFix, 'Duration', pp.minFixDurSecs, 'Timeout', pp.maxWaitForFixSecs);   %[isEsc, gazePos, err, targetIdx, startOfFix] = lpsy.etWaitForFixation('TargetPos', CntrPos.fovea, 'MaxNoise', maxNoiseDeg, 'MaxDev', maxDevDeg, 'DriftDelay', driftDelaySecs, 'StartTi', [], 'Timeout', [])        
        if isEsc             
            break   %break out of trial loop ending the experiment   
        elseif err
            error(err)
        else
            d.tToFixSecs(iTrial,1) = tStartOfFix-tFix0;
        end
        
        %%Loop through stimulus presentations
        for iStimPres = 1:pp.nStimTot
            %DRAW STIMFRAME
            %draw fixationpoint 
            DrawRing.DrawRing(sci.wnd, fxCntr, pp.fxColor, 2*pp.fxRadPx);   

                        
            %if stim will be presented left
            if (tt.isLeftStart(iTrial) && mod(iStimPres,2)) ...    
                    || (~tt.isLeftStart(iTrial) && ~mod(iStimPres,2))               
                
                % if Stim will be presented in left
                isLeftPresent = 1;
                
                %draw disks
                if pp.nDisks>=3      
                    DrawRing.DrawRing(sci.wnd, diskCntrs(:,1), pp.diskColor, pp.diskRadPx); %left disk
                end

                DrawRing.DrawRing(sci.wnd, diskCntrs(:,2), pp.diskColor, pp.diskRadPx);     %middle disk

                if pp.nDisks>=2
                    DrawRing.DrawRing(sci.wnd, diskCntrs(:,3), pp.diskColor, pp.diskRadPx); %right disk
                end
                
                %draw inner stimuli 
                if iStimPres == pp.nEmptyBefore+1 || iStimPres == pp.nEmptyBefore+tt.nEmptyBetween(iTrial)+2
                    % draw inner stim depends on the innerStimType
                    drawInnerStim('bars',sci.wnd,iStimPres,isLeftPresent,innerStimPos,pp.nDisks,tt.nEmptyBetween(iTrial))                    
                end

            %if stim will be presented right
            else       
                % if Stim will be presented in left
                isLeftPresent = 0;
                                
                %draw disks 
                DrawRing.DrawRing(sci.wnd, diskCntrs(:,2), pp.diskColor, pp.diskRadPx); %left disk
                               
                if pp.nDisks>=2
                    DrawRing.DrawRing(sci.wnd, diskCntrs(:,3), pp.diskColor, pp.diskRadPx);     %middle disk
                end
                
                if pp.nDisks>=3
                    DrawRing.DrawRing(sci.wnd, diskCntrs(:,4), pp.diskColor, pp.diskRadPx); %right disk
                end

                %draw dots
                if iStimPres == pp.nEmptyBefore+1 || iStimPres == pp.nEmptyBefore+tt.nEmptyBetween(iTrial)+2                       
                    % draw inner stim depends on the innerStimType
                    drawInnerStim('bars',sci.wnd,iStimPres,isLeftPresent,innerStimPos,pp.nDisks,tt.nEmptyBetween(iTrial))                    
                end
            end

            %FLIP STIM TO SCREEN
            if iStimPres == 1   %if this is the first stim pres of the trial
                %flip stim to screen as soon as lpsy.etWaitForFixation returns (and the stimulus drawing in the background is finished)                
                [tLstFlip,isEsc,isTimingErr(iStimPres)] = lpsy.flip(sci.wnd);

                %Start eyetracking recording as of etBaselineSecs before tLstFlip 
                %   -Local recording mode can go up to 1s back in time.         
                %   -The data collected during <etBaselineSecs> can be used to
                %   correct signal drift later, because the gaze is supposed to be
                %   always at the central fixation point (which is kind of certain
                %   since we make sure the trial doesn't start if it is not).
                lpsy.etRecordTrial(tLstFlip, etBaselineSecs, 'marker', 100+iTrial);
            else %if iStimPres>1 %If this is not the first stim pres of the trial                                
                if tt.isiDurFrms(iTrial)>0        %if ISI>0ms                    
                    %flip stim to screen after waiting isiDurSecs
                    [tLstFlip,isEsc,isTimingErr(iStimPres)] = lpsy.flip(sci.wnd, tLstFlip, tt.isiDurSecs(iTrial));  
                else%if p.isiDurFrms(iTrial)==0  %if ISI==0ms
                    %flip stim to screen after waiting stimDurSecs (i.e., without intermittent ISI)
                    [tLstFlip,isEsc,isTimingErr(iStimPres)] = lpsy.flip(sci.wnd, tLstFlip, tt.stimDurSecs(iTrial));
                end                    
            end
                
            if isEsc, break, end %if escape was pressed break out of stimPres loop                
            
            %DRAW ISI FRAME                      
            if tt.isiDurFrms(iTrial) > 0 && iStimPres < pp.nStimTot   %if ISI>0ms && this was not the last stim frame of the trial 

                %                     %For program testing: stops script until keypress
                %                     %before stimulus disappears for the ISI
                %                     lpsy.getKey(respKeys, 999);
                %                     lpsy.releaseKeyWait;
                %draw fixationpoint 
                if pp.isISFix
                    DrawRing.DrawRing(sci.wnd, fxCntr, pp.fxColor, 2*pp.fxRadPx); %Screen('FillRect', w, Clr.fx, fxRect);
                end
                
                if pp.isMask                    % mask presentation conditions
                    %flip stim to screen after waiting isiDurSecs
                    [tLstFlip,isEsc,isTimingErr(iStimPres)] = lpsy.flip(sci.wnd, tLstFlip, tt.maskDurSecs(iTrial));
                else                         % no mask condition
                    %FLIP ISI TO SCREEN after waiting for stimDurSecs                 
                    [tLstFlip,isEsc,isTimingErr(iStimPres)] = lpsy.flip(sci.wnd, tLstFlip, tt.stimDurSecs(iTrial));
                end 
                if isEsc, break, end %if escape was pressed break out of stimPres loop 
            end       
        end %stimPres loop
               
        %Clear screen after last stim is presented        
        [tLstFlip,isEsc,isTimingErr(iStimPres)] = lpsy.flip(sci.wnd, tLstFlip, tt.stimDurSecs(iTrial));                                
        d.isTimingErr(iTrial) = any(isTimingErr);
        if isEsc, break, end %if escape was pressed break out of trial loop 
        
        lpsy.etFinishTrial();   %Finish eyetracking recording for this trial
        
        %% Online fixation control                                
        if isOnlineFixCtrl && isEtConnected %~strcmpi(etDevice, 'none')            
            %check if the subject fixated well during the trial                
            dTrialSecs = pp.stimDurSecs(iTrial)*pp.nStimTot + pp.isiDurSecs(iTrial)*(nStimTot-1); % -1 because the trial starts with a stim onset and ends with a stim offset (ie., there is no ISI after the last disk)            
            fromTi = tLstFlip-dTrialSecs-pp.etBaselineSecs; %from last reference (timestamp of screen cleaning) back by the full trial duration and the duration of the baseline
            toTi = tLstFlip-pp.stimDurSecs(iTrial)*pp.nEmptyAfter-pp.isiDurSecs(iTrial)*pp.nEmptyAfter;      %don't analyze the recording during the dot-less disks at the end of the trial and the ISI presented before it (analysis ends with the offset of the last disk with dot)
            [d.isEtOk(iTrial), isEnufSmps4Anal, isTrialWiLongBadRun, isGoodFxTrial] = ctrlFx(sci.wnd, fromTi, toTi, sci, sysInfo);   
            if isEsc, break, end %if escape was pressed break out of trial loop
        else %if online fixation control is not used or no et is connected
            d.isEtOk(iTrial) = 1;   
        end %of isOnlineFixCtrl && isEtConnected
        
        %% Response  processing
        if d.isEtOk(iTrial) %if fixation was good
            
            %Get response: Wait max. <respTimeoutSecs> for one of the <respKeys> to be pressed
            [isEsc, keyIdx, t1] = lpsy.getKey(respKeys, pp.respTimeoutSecs);        
            if isEsc, break, end %if escape was pressed break out of trial loop 

            if keyIdx %if any button was pressed
                d.rt(iTrial) = t1-tLstFlip;
                d.keyIdx(iTrial) = keyIdx;
                if keyIdx ==2, d.keyIdx(iTrial) = -1; end
                                
                switch pp.task
                    case 'priority_map'
                        % priority map experiment is a discrimination task,
                        % pp.cue_target_tilt(iTrial) * pp.prob_target_tilt(iTrial) <0 :d.keyIdx(iTrial)==2 (-1)
                        % pp.cue_target_tilt(iTrial) * pp.prob_target_tilt(iTrial) >0 :d.keyIdx(iTrial)==1
                        if tt.cue_target_tilt(iTrial) * tt.prob_target_tilt(iTrial) * d.keyIdx(iTrial) > 0
                            d.isCorr(iTrial) = true;
                        else
                            d.isCorr(iTrial) = false;

                            if pp.isIncorrFeedb 
                                Beeper(600, 1, .1) %Beeper(frequency, [fVolume], [durationSec]);    
                                WaitSecs(1);
                            end                
                        end
                        
                    case 'LilacChaser'
                        % Lilac Chaser is a motion or green circle detection task. Either motion or green circle perceived?
                        if d.keyIdx(iTrial) == 1
                            d.isCorr(iTrial) = false;
                        elseif d.keyIdx(iTrial) == 2
                            d.isCorr(iTrial) = true;
                        end
                    otherwise    
                        error('Task %s unknown. Check spelling.', task)       
                end                                                
            else
                %nonresp handling
                if pp.isTimeoutFeedb
                    Beeper(300, 1, .8) %Beeper(frequency, [fVolume], [durationSec]);
                    WaitSecs(1);
                end

                % Repeat trial at random later moment and
                % continue with next trial without waiting for a response
                tt = repeat_nonresp(iTrial, tt);
            end %of if d.keyIdx(iTrial) (any button pressed)    
            
        else %if fixation was broken
            % Repeat trial at random later moment and
            % continue with next trial without waiting for a response
            tt = repeat_nonresp(iTrial, tt);                    
        end %of isEtOk(iTrial)
                             
        %Advance trial counter
        iTrial = iTrial+1;
        
    end % end of the trial block
    
    iTrial = iTrial-1; %undo the very last trial counter advancement
    
    %% Regular end of experiment (also executed when esc is pressed)    
    if pp.isSaveData
        myData = saveDvFile(iTrial);
        saveEpdFile();
        nTrials = iTrial-sum(d.keyIdx(1:iTrial)==0);
        if nTrials == pp.numCal 
            sessionSucceed(pp.startTrialID+nTrials); % The last tested trial = iTrial - missed trials (d.keyIdx==0)
        end
    end        
    
    lpsy.cleanup();    
    sprintf('\tElapsed time: %.0fm%.0fs',floor(toc(tBlockStart)/60), rem(toc(tBlockStart),60))
    
% if error, clean up
catch err
    lpsy.cleanup(err);
    
    %% Irregular end of experiment            
    if pp.isSaveData
        myData = saveDvFile(iTrial);    
        saveEpdFile();        
    end
        
    %print to command-prompt
    fprintf(1,'CATCH\n');
    
    %display error
    rethrow(err)  
end

end % end of TPD_PriorityMap_Stim()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		INNER FUNCTIONS 		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Saving data

function myData = saveDvFile(iTrial)
    global pp tt d sci
    ii = [1:iTrial];

    %create datatable to be written to dv-file (or be dumped to disk as .mat if dv file writing fails)
    myData.TRIAL_ID = ii';  
    myData.stimElev   = repmat(pp.stimElev,iTrial,1);
    myData.isiDurFrms = tt.isiDurFrms(ii);
    myData.isiDurSecs = tt.isiDurSecs(ii);        
    myData.maskDurFrms = tt.maskDurFrms(ii);
    myData.maskDurSecs = tt.maskDurSecs(ii);
    myData.stimDurFrms = tt.stimDurFrms(ii);
    myData.stimDurSecs = tt.stimDurSecs(ii); 
    myData.nEmptyBetween = tt.nEmptyBetween(ii);
    myData.isLeftStart = tt.isLeftStart(ii);
    myData.isNonRet    = tt.isNonRet(ii);
    myData.condition   = tt.condition(ii); 
    myData.distractor_pos = tt.distractor_pos(ii);
    myData.cue_target_pos = tt.cue_target_pos(ii);
    myData.prob_target_pos= tt.prob_target_pos(ii);
    myData.cue_target_tilt = tt.cue_target_tilt(ii);
    myData.cue_distractor_tilt = tt.cue_distractor_tilt(ii);
    myData.prob_target_tilt = tt.prob_target_tilt(ii);
    myData.keyIdx = d.keyIdx(ii);
    myData.isCorr = d.isCorr(ii);        
    myData.rt = d.rt(ii);        
    myData.tToFixSecs = d.tToFixSecs(ii);            
    myData.isEtOk = d.isEtOk(ii);
    myData.isAnyTimingErr = d.isTimingErr(ii);        

    %open file with 'w'rite permission, or create it if it doesn't exist
    fid = fopen(pp.saveFilename,'w');                 

    if fid>=0

        lpsy.writeDvHeader(fid); %Write header to file

        %% Collect setup parameters            
        %   lpsy.writeDvPmtr(description,value,unit,prec)
        %   lpsy.writeDvPmtr(description, string)

        %Program settings
        %!! bools converted to double                                    
        lpsy.writeDvPmtr('isDebug', double(pp.isDebug), 'log. true-false', 0); 
        lpsy.writeDvPmtr('isSaveData', double(pp.isSaveData), 'log. true-false', 0);       
        lpsy.writeDvPmtr('isMask', double(pp.isMask), 'log. true-false', 0);
        lpsy.writeDvPmtr('condNm', pp.cond);             
        
        lpsy.writeDvPmtr('task', pp.task); 
        lpsy.writeDvPmtr('nDisks', pp.nDisks, '', 0);
        
        lpsy.writeDvPmtr('isIncorrFeedb', double(pp.isIncorrFeedb), '(log. true-false)', 0);
        lpsy.writeDvPmtr('isTimeoutFeedb', double(pp.isTimeoutFeedb), '(log. true-false)', 0);            

        %Screen settings
        lpsy.writeDvPmtr('sci.gamma', sci.gamma, '(1 for linear)', 4);
        lpsy.writeDvPmtr('sci.gammaMode', sci.gammaMode); %string!
        lpsy.writeDvPmtr('sci.multisample', sci.multisample, '', 0);
        lpsy.writeDvPmtr('sci.wdthPx', sci.width, 'px', 0);
        lpsy.writeDvPmtr('sci.hghtPx', sci.height, 'px', 0); 
        lpsy.writeDvPmtr('sci.freqHz', sci.hz, 'hz', 0);
        lpsy.writeDvPmtr('pp.freqHz', pp.freqHz, 'hz', 0);
        lpsy.writeDvPmtr('Viewing distance', sci.vd_m, 'm', 2);

        %Repetition settings
        lpsy.writeDvPmtr('nTrials', pp.numCal, '', 0);
        lpsy.writeDvPmtr('nStimPres', pp.nStimPres, '', 0);
        lpsy.writeDvPmtr('nEmptyBefore', pp.nEmptyBefore, '', 0);
        lpsy.writeDvPmtr('nEmptyBetween', pp.nEmptyBetween(1), '', 0);
        lpsy.writeDvPmtr('nEmptyAfter', pp.nEmptyAfter, '', 0);
        lpsy.writeDvPmtr('nStimTot', pp.nStimTot, '', 0);

        %Temporal settings            
        lpsy.writeDvPmtr('isiDursToBeUsedFrms', [num2str(pp.isiDursToBeUsedFrms) '   frames']);
        lpsy.writeDvPmtr('isiDursToBeUsedSecs', [sprintf('%.3f   ', pp.isiDursToBeUsedFrms/pp.freqHz) 'secs']);            
        lpsy.writeDvPmtr('maskDursToBeUsedFrms', [num2str(pp.maskDursToBeUsedFrms) '   frames']);
        lpsy.writeDvPmtr('maskDursToBeUsedSecs', [sprintf('%.3f   ', pp.maskDursToBeUsedFrms/pp.freqHz) 'secs']);        
        lpsy.writeDvPmtr('stimDursToBeUsedFrms', [num2str(pp.stimDursToBeUsedFrms) '   frames']);
        lpsy.writeDvPmtr('stimDursToBeUsedSecs', [sprintf('%.3f   ', pp.stimDursToBeUsedFrms/pp.freqHz) 'secs']);
        lpsy.writeDvPmtr('itiDurFrms', pp.itiDurFrms, 'frames', 0);
        lpsy.writeDvPmtr('itiDurSecs', pp.itiDurFrms/pp.freqHz, 'secs', 3);
        lpsy.writeDvPmtr('respTimeoutSecs', pp.respTimeoutSecs, 'secs', 2);

        %Spatial settings
        lpsy.writeDvPmtr('fxRadPx', pp.fxRadPx, 'px', 0);
        lpsy.writeDvPmtr('diskRadPx', pp.diskRadPx, 'px', 0); 
        lpsy.writeDvPmtr('stimShiftPx', pp.stimShiftPx, 'px', 0);
        lpsy.writeDvPmtr('stimElevPx', pp.stimElevPx, 'px', 0);
        lpsy.writeDvPmtr('barWidthPx', pp.barWidthPx, 'px', 0);            

        %Color settings (matrices as strings)
        lpsy.writeDvPmtr('backgrColor', mat2str(pp.backgrColor));
        lpsy.writeDvPmtr('fxColor', mat2str(pp.fxColor));
        lpsy.writeDvPmtr('diskColor', mat2str(pp.diskColor));
        lpsy.writeDvPmtr('maskColor', mat2str(pp.maskColor));    

        %Eyetracker settings
        lpsy.writeDvPmtr('etDevice', pp.etDevice);
        lpsy.writeDvPmtr('isOnlineFixCtrl', double(pp.etOnline), 'logical', 0)

        %Eyetracker settings for lpsy.etWaitForFixation at trialstart
        lpsy.writeDvPmtr('maxNoiseDeg', pp.maxNoiseDeg, 'arcdeg', 4);
        lpsy.writeDvPmtr('maxDevDeg', pp.maxDevDeg, 'arcdeg', 4);
        lpsy.writeDvPmtr('driftDelaySecs', pp.driftDelaySecs, 'secs', 4);
        lpsy.writeDvPmtr('tStartLookForFix', pp.tStartLookForFix, 'arcdeg', 4);
        lpsy.writeDvPmtr('minFixDurSecs', pp.minFixDurSecs, 'secs', 4);
        lpsy.writeDvPmtr('maxWaitForFixSecs', pp.maxWaitForFixSecs, 'secs', 4);

        lpsy.writeDvPmtr(fid) %Write parameters to file                                    
        lpsy.writeDvTable(fid, 0, myData); %Write datatable to file

        %close the file
        ST = fclose(fid);

        %print confirmation of the saving to command prompt
        [dataDir, fn, ext] = fileparts(pp.saveFilename);
        fn = [fn, ext];            
        if ST == 0
            fprintf(1, '%s \n\t Dv-File saving succesful. \n\t File: %s \n\t Path: %s \n', ...
            sprintf('%d-%02.0f-%02.0f, %02.0fh%02.0fm%02.0fs', clock), ...
            fn, dataDir);
        else
            save(pp.saveFilename, myData);
            warning('Something went wrong while saving the dv file (ST~=0)!')
            warning('The myData structure has been dumped to the data directory.')
            warning('Please notify marc.lauffs@epfl.ch')
        end

    else
        %-create some redundancy in case file opening fails-
        %dlmwrite([ffn '.txt'], myData, '\t');
        save(pp.saveFilename, myData);
        fprintf(1,['DV-FILE SAVING FAILED!! \nThe myData '...
            'structure has been dumped to the data directory.' ...
            '\nPlease notify marc.lauffs@epfl.ch.']);
    end
end %of saveDvFile()

%--------------------------------------------------------------------

function epdSaveSuccess = saveEpdFile()
    global pp

    epdSaveSuccess = lpsy.etSaveRecording([pp.saveFilename(1:end-3)]); %.epd is appended automatically

    if strcmpi(pp.etDevice, 'none')
        fprintf(1, '\t No eyetracking was used. \n')
    elseif epdSaveSuccess % and an <etDevice> was specified
        fprintf(1, '\t Eyetracking data saving succesful. \n')
    else % an <etDevice> was specified and saving failed (~epdSaveSuccess)
        warning('Eyetracking data saving failed!')
        warning('Please notify marc.lauffs@epfl.ch')
    end
end


%% Initialization related functions
function diskCntrs = getStimCntrs(stimCntr, intStimDistPx)
    % Get stim centers for use with MR's DrawRing.DrawRing or Screen('DrawDots')
    % (Not needed if PTB's 'FillOval' function is used)                

    %get disk centers 
    diskCntrs = NaN(2,4);  % x/y (2) per disc (4)
    diskCntrs(:,1) = [stimCntr(1)-1.5*intStimDistPx; stimCntr(2)];
    diskCntrs(:,2) = [stimCntr(1)-0.5*intStimDistPx; stimCntr(2)];
    diskCntrs(:,3) = [stimCntr(1)+0.5*intStimDistPx; stimCntr(2)];
    diskCntrs(:,4) = [stimCntr(1)+1.5*intStimDistPx; stimCntr(2)];
          
end % end of getStimCntrs()

function [p] = repeat_nonresp(iTrial, p)
    %Shuffles the order of the remaining trials (incl. the current timeout trial)
    %  and pastes them under the current one, 
    %  thereby extending p.(pfields{i}) by 1.
    %  Leaves already run trials untouched.
    %
    % p is a struct array with one field per parameter. Each field
    % contains a matrix where each row contains the settings for one
    % trial.

    % Copy-pasted from TPD_TmpCrwd_v3

    pfields = fieldnames(p);
    %dfields = fieldnames(d);

    nTrialsNow = length(p.(pfields{1}));    
    
    newInd   = randi([min(iTrial+1,nTrialsNow), nTrialsNow]);
    newOrder = [iTrial+1:newInd,iTrial,newInd+1:nTrialsNow];
    
    for i = 1:numel(pfields)            
        p.(pfields{i})(iTrial+1 : nTrialsNow+1, :) = p.(pfields{i})(newOrder,:);  
    end

end %of repeat_nonresp() 
%% Draw stimuli related functions

function StimPos = getInnerStimPos(innerStimType, iTrial, stimCntr, nStimPres)
    global pp tt
    % Set Inner stimuli positions depends on its inner stimuli type (innerStimType)
    switch innerStimType
        case 'bars'
            StimPos = getBarPos(tt.distractor_pos(iTrial), tt.cue_target_pos(iTrial), tt.prob_target_pos(iTrial), tt.cue_target_tilt(iTrial), tt.cue_distractor_tilt(iTrial), tt.prob_target_tilt(iTrial),tt.isNonRet(iTrial), pp.diskRadPx*2, stimCntr', pp.diskRadPx*2.5, nStimPres);
        case 'barMask'
            StimPos = getBarMaskPos(pp.diskRadPx*2, stimCntr', pp.diskRadPx*2.5, pp.nStimPres);
    end 
end % end of getStimPos()

function StimPos = getBarPos(distractor_pos,cue_target_pos,prob_target_pos,cue_target_tilt, cue_distractor_tilt,prob_target_tilt,isNonRet,diskDiaPix,disk_center,disk_gap, nStimPres)
    % Set bar positions per disk, output fields are leftDisk, middleDisk, rightDisk
    % output StimPos: iStimPres x 2 x num_of_bars*2 (matrix configuration: [start end start end ...])

    if nStimPres ~= 2 % for getBarPos() only 2 frame case is implemented in this moment!
       error('Only 2 frame case is implemented in this moment! Please update the code to have more frames')
    end

    % bar unit 
    bar_unit     = diskDiaPix/15;
    inner_ring_r = 2*diskDiaPix/10;         % inner ring radius
    inner_ring_d = reshape(repmat([60,180,300,0,120,240],[2,1]),1,[]); % inner ring degrees, repeated to make start & end point of the bar % right-top:1,2 % left-top:3,4 % left-bottom:5,6 % right-bottom:7,8
    outer_ring_r = 4*diskDiaPix/10;         % outer ring radius
    outer_ring_d = reshape(repmat([0:30:359],[2,1]),1,[]);              % outer ring degrees, repeated to make start & end point of the bar

    % Bar positions based on [0;0] center
    inner_ring   = inner_ring_r*[sind(inner_ring_d) 0 0; cosd(inner_ring_d) 0 0]; % inner ring bars + center bar [x_start, x_end ...; y_start y_end ...]
    outer_ring   = outer_ring_r*[sind(outer_ring_d); cosd(outer_ring_d)];          % outer ring bars [x_start, x_end ...; y_start y_end ...]
    default_bar  = repmat([0  0; 1 -1]*bar_unit,[1,(length(inner_ring)+length(outer_ring))/2]);

    % setup overall bar positions for neutral, distractor, and target disk
    neutral_disk = cat(2,inner_ring,outer_ring)+default_bar; % no tilted bars [x_start, x_end ...; y_start y_end ...]
    cue_disk     = neutral_disk;
    prob_disk    = neutral_disk;

    % Cue and probe disk's bar positions 
    %   distractor_position, cue_target_position, probe_target_position
    %   cue_target_tilt, probe_tilt
    % Add distractor
    cue_disk (:,(distractor_pos-1)*2+1:distractor_pos*2)   = inner_ring(:,(distractor_pos-1)*2+1:distractor_pos*2) + bar_unit*[[-sind(45) +sind(45)]*cue_distractor_tilt; -cosd(45) +cosd(45)];
    % Add targets
    cue_disk (:,(cue_target_pos-1) *2+1:cue_target_pos*2)  = inner_ring(:,(cue_target_pos-1 )*2+1:cue_target_pos*2 ) + bar_unit*[[+sind(15) -sind(15)]*cue_target_tilt ; -cosd(15) +cosd(15)];
    prob_disk(:,(prob_target_pos-1)*2+1:prob_target_pos*2) = inner_ring(:,(prob_target_pos-1)*2+1:prob_target_pos*2) + bar_unit*[[+sind(15) -sind(15)]*prob_target_tilt; -cosd(15) +cosd(15)];

    % frame 1: Neutral, frame2: Neutral/Distractor, frame3: Target/Neutral, frame4: Neutral
    % input needed : disk_center xy coordinate (in between mid-left/mid-right), disk_gap 
    center_position = repmat(disk_center,[1,4])+[-1.5*disk_gap,-0.5*disk_gap,+0.5*disk_gap,+1.5*disk_gap;0,0,0,0]; % [x1,x2,x3,x4;y1,y2,y3,y4] left, mid-left, mid-right, right xy coordinate 

    % Set-up xy coordinates in StimPos.left/middle/right
    StimPos.leftDisk   = cat(3,center_position(:,2),center_position(:,1));
    StimPos.middleDisk = cat(3,center_position(:,3),center_position(:,2));
    StimPos.rightDisk  = cat(3,center_position(:,4),center_position(:,3));

    if isNonRet
        frame2MidDisk = cue_disk;
        frame2lefDisk = neutral_disk;
    else
        frame2MidDisk = neutral_disk;
        frame2lefDisk = cue_disk;
    end

    StimPos.leftDisk   = repmat(StimPos.leftDisk,   [1,length(neutral_disk),1]) + cat(3,frame2lefDisk,neutral_disk);
    StimPos.middleDisk = repmat(StimPos.middleDisk, [1,length(neutral_disk),1]) + cat(3,frame2MidDisk,prob_disk);
    StimPos.rightDisk  = repmat(StimPos.rightDisk,  [1,length(neutral_disk),1]) + cat(3,neutral_disk,neutral_disk);

end % end of getBarPos()

function StimPos = getBarMaskPos(diskDiaPix, disk_center, disk_gap, nStimPres)
    % Set bar positions per disk, output fields are leftDisk, middleDisk, rightDisk
    % output StimPos: iStimPres x 2 x num_of_bars*2 (matrix configuration: [start end start end ...])

    if nStimPres ~= 2 % for getBarPos() only 4 frame case is implemented in this moment!
       error('Only 2 frame case is implemented in this moment! Please update the code to have more frames')
    end

    % bar unit 
    bar_unit     = diskDiaPix/15;
    inner_ring_r = 2*diskDiaPix/10;         % inner ring radius
    inner_ring_d = reshape(repmat([25,155,205,335,90,270],[2,1]),1,[]); % inner ring degrees, repeated to make start & end point of the bar % right-top:1,2 % left-top:3,4 % left-bottom:5,6 % right-bottom:7,8
    outer_ring_r = 4*diskDiaPix/10;         % outer ring radius
    outer_ring_d = reshape(repmat([0:30:359],[2,1]),1,[]);              % outer ring degrees, repeated to make start & end point of the bar

    % Bar positions based on [0;0] center
    inner_ring   = inner_ring_r*[sind(inner_ring_d) 0 0; cosd(inner_ring_d) 0 0]; % inner ring bars + center bar [x_start, x_end ...; y_start y_end ...]
    outer_ring   = outer_ring_r*[sind(outer_ring_d); cosd(outer_ring_d)];          % outer ring bars [x_start, x_end ...; y_start y_end ...]

    % setup overall bar positions for mask disk
    mask_disk_pos= cat(2,inner_ring,outer_ring); % mask disk bar's positions [x_start, x_end ...; y_start y_end ...]

    % frame 1: Neutral, frame2: Neutral/Distractor, frame3: Target/Neutral, frame4: Neutral
    % input needed : disk_center xy coordinate (in between mid-left/mid-right), disk_gap 
    center_position = repmat(disk_center,[1,4])+[-1.5*disk_gap,-0.5*disk_gap,+0.5*disk_gap,+1.5*disk_gap;0,0,0,0]; % [x1,x2,x3,x4;y1,y2,y3,y4] left, mid-left, mid-right, right xy coordinate 

    % Set-up xy coordinates in StimPos.left/middle/right
    StimPos.leftDisk   = repmat(cat(3,center_position(:,1),center_position(:,2),center_position(:,1),center_position(:,2)), [1,length(mask_disk_pos),1]);
    StimPos.middleDisk = repmat(cat(3,center_position(:,2),center_position(:,3),center_position(:,2),center_position(:,3)), [1,length(mask_disk_pos),1]);
    StimPos.rightDisk  = repmat(cat(3,center_position(:,3),center_position(:,4),center_position(:,3),center_position(:,4)), [1,length(mask_disk_pos),1]);

    nBar2x = size(mask_disk_pos,2); % num of bars *2

    for ii = 1:nStimPres
        % get random bar orientations for each frames
        bar_ori = randi(360, 3, nBar2x/2);
        mask_disk_ori = zeros(2,nBar2x,3);
        mask_disk_ori(:,:,1) = reshape(reshape(cat(3, [-sind(bar_ori(1,:)); +sind(bar_ori(1,:));],[-cosd(bar_ori(1,:));cosd(bar_ori(1,:))]),2,[]),nBar2x,[])'; % clockwise orientation [[-sind;+sind];[-cosd;+cosd]] -> reshape(2,[]) [-sind, -cosd; sind, cosd] -> reshape(nBar2x,[])' [-sind, sind; -cosd, cosd]
        mask_disk_ori(:,:,2) = reshape(reshape(cat(3, [-sind(bar_ori(2,:)); +sind(bar_ori(2,:));],[-cosd(bar_ori(2,:));cosd(bar_ori(1,:))]),2,[]),nBar2x,[])'; 
        mask_disk_ori(:,:,3) = reshape(reshape(cat(3, [-sind(bar_ori(3,:)); +sind(bar_ori(3,:));],[-cosd(bar_ori(3,:));cosd(bar_ori(1,:))]),2,[]),nBar2x,[])'; 

        % xy coordinates for the mask bars
        StimPos.leftDisk(:,:,ii)   = StimPos.leftDisk(:,:,ii)   + mask_disk_pos + bar_unit*mask_disk_ori(:,:,1);
        StimPos.middleDisk(:,:,ii) = StimPos.middleDisk(:,:,ii) + mask_disk_pos + bar_unit*mask_disk_ori(:,:,2);
        StimPos.rightDisk(:,:,ii)  = StimPos.rightDisk(:,:,ii)  + mask_disk_pos + bar_unit*mask_disk_ori(:,:,3);
    end

end % end of getBarMaskPos()

function drawInnerStim(innerStimType,w,iStimPres,isLeftPresent,innerStimPos,nDisks,nEmptyBetween)
    global pp
    % Draw inner stimuli depends on innerStimType

    switch innerStimType

        case 'bars'
            disp(iStimPres)
            if iStimPres <= pp.nEmptyBefore + nEmptyBetween +1
                iStimPresFrame = iStimPres - pp.nEmptyBefore;
            elseif iStimPres > pp.nEmptyBefore + nEmptyBetween +1
                iStimPresFrame = iStimPres - pp.nEmptyBefore - nEmptyBetween;
            end 
            if isLeftPresent % Stimuli are presented left
                if pp.nDisks>=3
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.leftDisk(:,:,iStimPresFrame),pp.barWidthPx,pp.barColor,[0 0], 1); %  Screen(DrawLines, windowPtr, xy [,width=2] [,colors] [,center] [,smooth=2][,lenient]) % Bars in the Left disk % -nEmptyBefore because iStimPres is already (nEmptyBefore) when the loop is read for the first time; +iStimPres to advance the dot position after every stimulus presentation;
                end
                [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.middleDisk(:,:,iStimPresFrame),pp.barWidthPx,pp.barColor,[0 0], 1); % Bars in the middle disk
                if pp.nDisks>=2
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.rightDisk(:,:,iStimPresFrame),pp.barWidthPx,pp.barColor,[0 0], 1); % Bars in the right disk
                end
            else % Stimuli are presented right
                [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.leftDisk(:,:,iStimPresFrame),pp.barWidthPx,pp.barColor,[0 0], 1); % Bars in the left disk
                if pp.nDisks>=2
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.middleDisk(:,:,iStimPresFrame),pp.barWidthPx,pp.barColor,[0 0], 1); % Bars in the middle disk
                end
                if pp.nDisks>=3
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.rightDisk(:,:,iStimPresFrame),pp.barWidthPx,pp.barColor,[0 0], 1); % Bars in the right disk
                end
            end

         case 'barMask'
            if isLeftPresent % Stimuli are presented left
                if pp.nDisks>=3
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.leftDisk(:,:,iStimPres-nEmptyBefore),pp.barWidthPx,pp.maskColor,[0 0], 1); %  Screen(DrawLines, windowPtr, xy [,width=2] [,colors] [,center] [,smooth=2][,lenient]) % Bars in the Left disk % -nEmptyBefore because iStimPres is already (nEmptyBefore) when the loop is read for the first time; +iStimPres to advance the dot position after every stimulus presentation;
                end
                [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.middleDisk(:,:,iStimPres-nEmptyBefore),pp.barWidthPx,pp.maskColor,[0 0], 1); % Bars in the middle disk
                if pp.nDisks>=2
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.rightDisk(:,:,iStimPres-nEmptyBefore),pp.barWidthPx,pp.maskColor,[0 0], 1); % Bars in the right disk
                end
            else % Stimuli are presented right
                [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.leftDisk(:,:,iStimPres-nEmptyBefore),pp.barWidthPx,pp.maskColor,[0 0], 1); % Bars in the left disk
                if pp.nDisks>=2
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.middleDisk(:,:,iStimPres-nEmptyBefore),pp.barWidthPx,pp.maskColor,[0 0], 1); % Bars in the middle disk
                end
                if pp.nDisks>=3
                    [minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines',w,innerStimPos.rightDisk(:,:,iStimPres-nEmptyBefore),pp.barWidthPx,pp.maskColor,[0 0], 1); % Bars in the right disk
                end
            end  
    end
end

%% Eye tracking related
function [isEtOk, isEnufSmps4Anal, isTrialWiLongBadRun, isGoodFxTrial] = ctrlFx(w, fromTi, toTi, Scrn, sysInfo)

    %get data from buffer
    [etData, etTimestamps, etInfo] = lpsy.etGetSignal(fromTi, toTi);  %[signal, timestamps, info] = lpsy.etGetSignal(fromTi,toTi);   

    %% Settings                                
    sizeOnePxInM = sysInfo.stimulusmonitor.widthmm_fallback/Scrn.width/1000;  %sizeOnePxInM = 0.2768/1000; %pixel pitch on AsusVG248QE according to manufacturer
    sizeOnePxInDeg = 2*atand(sizeOnePxInM/ (2*Scrn.vd_m));                     %0.0240 (0.02402950232842317905289130831128) @0.66m (as well)
    scrnMaxXDeg = Scrn.width/2*sizeOnePxInDeg; %max distance of screen center to screen edge in x direction (arcdeg)
    scrnMaxYDeg = Scrn.height/2*sizeOnePxInDeg; %max distance of screen center to screen edge in y direction (arcdeg)

    etHz = etInfo.sampFreq;
    dSmpMs = 1000/etHz;                     %duration of 1smp in ms      
    vThresh4Artif = 750;                    %[arcdeg/sec]
    dBadRunToPadMs = 20;                    %Min. duration of a run of bad samples that is padded (in ms).
    dPadMs = 30;                            %Number of ms to pad before and after bad signal runs.
    pBadSmpThresh = 33;                     %Threshold for the proportion of samples with bad signal, before the trial is excluded [33]
    isExclTrialsWiRunsOfBadSignal = true;
    badRunThreshMs = 250;                   %exclude trial if it contains a run of >= <badRunThreshMs> ms with bad signal. Choose long enough so that a blink can happen in this time (so ca. 100-150ms PLUS PADDING, 250ms seems to work ok). Fixations are ca. 200-300ms.
    fixCtrlOkRadDeg = 1.5;                  %a valid (~isBadSignal) POR > abs(<fixCtrlOkRadDeg>) is counted as "gaze is not on fixationpoint" [1.5]
    fixCtrlTiThreshMs = 20;                 %fixation is considered broken when in >= <fixCtrlTiThreshMs> the POR is valid (~isBadSignal) and > abs(<fixCtrlOkRadDeg>) 

    %---               

    %convert settings from ms to smps
    dBadRunToPadSmps = dBadRunToPadMs/dSmpMs;   %Min. duration of a run of bad samples that is padded (in samples). [10]
    dPadSmps = dPadMs/dSmpMs;                   %Number of samples to pad before and after bad signal runs. [15]
    badRunThreshSmps = badRunThreshMs/dSmpMs;  %exclude trial if it contains a run of >= <badRunThreshSmps> samples with bad signal. Choose long enough so that a blink can happen in this time (so ca. 100-150ms PLUS PADDING, 250ms seems to work ok). Fixations are ca. 200-300ms. [75]
    fixCtrlTiThreshSmps = fixCtrlTiThreshMs/dSmpMs;  %fixation is considered broken when in >= <fixCtrlTiThreshSmps> the POR is valid (~isBadSignal) and > abs(<fixCtrlOkRadDeg>) [10]

    %get info from the etData        
    dEpochsSmps = size(etData,2);
    porX = etData(1,:);
    porY = etData(2,:);
    pupDia = etData(3,:);

    %calculate eye-velocity from non-smoothed position data
    v3d = sqrt((diff(porX).^2) + (diff(porY).^2)) *etHz; %speed in deg-per-second from (potentially) smoothed POR data                                    
    %prepend a single NaN to compensate for the sample lost in the diff() operation (the result of diff(x) is one shorter than x)
    v3d = [NaN, v3d];        
    %compare eye-velocity to threshold for biologically plausible velocities
    isPlausibleV = ~(abs(v3d)>vThresh4Artif);       %the complicated calc with ~ is because NaN > int is always false, but we want it true here            

    %check if por is on screen                        
    isPorOnScreen = (porX > -scrnMaxXDeg & porX < scrnMaxXDeg ...
        & porY > -scrnMaxYDeg & porY < scrnMaxYDeg);

    %per sample indication if data are unusable
    isBadSignal = (~isPorOnScreen | ~isPlausibleV | pupDia <= 0.001 ...
        | pupDia == etInfo.NOVAL ...
        | pupDia == etInfo.INVALID...
        | porX == etInfo.NOVAL ...
        | porX == etInfo.INVALID ...
        | porY == etInfo.NOVAL ...
        | porY == etInfo.INVALID);    

    %% Padding of bad samples
    %-Strategy: Identify runs of at least <dBadRunToPadSmps> and
    %pad around each sample that is member of the run 
    %(around each sample so we don't have to separate and identify the end of the runs) 
    badRunsStartIdxs = strfind(isBadSignal, true(1,dBadRunToPadSmps));                                                           

    %define start and end points for padding (in smps)
    startPad = badRunsStartIdxs - dPadSmps;
    stopPad  = badRunsStartIdxs+dBadRunToPadSmps + dPadSmps;
    %make sure we only work on samples within the range of the epoch (i.e., 1:dEpochsSmps). 
    %  This is necessary to not change the lenght of isBadSignal or produce errors because we try to use indices <1.
    startPad(startPad<1) = 1;                   %indexes smaller 1 start being padded at 1 and not before
    stopPad(stopPad>dEpochsSmps) = dEpochsSmps;   %indexes larger than dEpochsSmps (i.e., exceeding the no of samples in the trial) are padded up to dEpochsSmps                        
    %flag samples to be padded as bad
    for i = 1:length(startPad)
        isBadSignal(startPad(i):stopPad(i)) = true;                
    end

    %% Data quality check         
    pBadSmps = (1-(sum(~isBadSignal)/dEpochsSmps))*100;      %calculate percentage of smps with bad data quality
    if pBadSmps > pBadSmpThresh                              %if allowed threshold is exceeded
        isEnufSmps4Anal = false;                      
    else
        isEnufSmps4Anal = true;       
    end

    %% Find long periods of data loss (during which an eye movement might have occured)
    a = strfind(isBadSignal, true(1,badRunThreshSmps)); %find runs of at least <badRunThreshSmps> for which isBadSignal is true
    if ~isempty(a)
        isTrialWiLongBadRun = true;
    else
        isTrialWiLongBadRun = false;
    end
    clear a

    %% Fixation control
    %find runs of at least <fixCtrlTiThreshSmps> where the gaze data are valid but the POR is not within fixation+-fixCtrlOkRadDeg 
    isPorNotOnFx = ~isBadSignal & (abs(porX)>fixCtrlOkRadDeg | abs(porY)>fixCtrlOkRadDeg); %check for each sample if the POR is valid and outside <fixCtrlOkRadDeg> around 0                        
    a = strfind(isPorNotOnFx, true(1,fixCtrlTiThreshSmps)); %find runs
    if ~isempty(a)
        isGoodFxTrial = false;
    else
        isGoodFxTrial = true;
    end
    clear a

    %% Trial exclusion -> Set entire trial or only bad samples to NaN
    isEtOk = ~any([~isEnufSmps4Anal, ...
        isTrialWiLongBadRun*isExclTrialsWiRunsOfBadSignal, ...
        ~isGoodFxTrial]);


    %% Feedback if eyetracking data are not ok for this trial
    if ~isEtOk                    
        %determine reason and adapt message
        if ~isGoodFxTrial
            msg = 'Eyemovement detected! \n Keep your eyes on the red point! \n Trial will be repeated at a random later moment. \n Press a button to continue.';
            etFeedback.freq = 600;
            etFeedback.vol = 0.6;
            etFeedback.durSecs = 0.8;                    
            etFeedback.clr = [255 0 0];
        elseif ~isEnufSmps4Anal 
            msg = 'Oops, low eyetracking data quality. This is not your fault. \n Trial will be repeated at a random later moment. \n Press a button to continue.';
            etFeedback.freq = 1000;
            etFeedback.vol = 0.1;
            etFeedback.durSecs = 0.1;                    
            etFeedback.clr = [255 145 0];
        elseif isTrialWiLongBadRun
            msg = 'Oops, eyes closed or data quality issue. \n (Long run of bad data.) \n Trial will be repeated at a random later moment. \n Press a button to continue.';
            etFeedback.freq = 1000;
            etFeedback.vol = 0.1;
            etFeedback.durSecs = 0.1;
            etFeedback.clr = [255 145 0];
        end

        %deliver visual feedback                    
        [~, ~, bbox] = DrawFormattedText(w, msg,'center','center', [250 0 0], [], [],[],3); %dummy call to get bbox                    
        Screen('FillRect',w, [0 0 0], bbox + [-10, -30, +13, 0]);
        Screen('FrameRect',w, etFeedback.clr, bbox + [-10, -30, +13, 0], 5);
        [~, ~, bbox] = DrawFormattedText(w, msg,'center','center', etFeedback.clr, [], [],[],3); %actual call, to be seen by subject
        lpsy.flip(w);

        %deliver auditory feedback
        Beeper(etFeedback.freq, etFeedback.vol, etFeedback.durSecs) %Beeper(frequency, [fVolume], [durationSec]);            

        %Wait for a key press, but continue automatically after respTimeoutSecs
        lpsy.releaseKeyWait(); %wait until all buttons are released 
        [isEsc] = lpsy.getKey(respKeys, 5);                    
        lpsy.flip(w);          %flip to remove message, as visual confirmation to button press
        if isEsc, return, end %if escape was pressed return to invoking function                             
    end                                           
end %of ctrlFx()

%% Make sessionSucceed file: save TotTrialID and use as session indicator
function ismade = sessionSucceed(finTrialID)
    % save sessionSucceed.txt file only in the main session
    global pp
    
    ismade = 0;
    
    % if training session, 
    if any(strfind(lower(pp.cond),lower('train')))
        return;
    end
    
    % open sessionSucceed.txt file: import txt filename from RUNfile (pp.savetxtname)
    txtid = fopen(pp.saveTxtname, 'w');
    
    fprintf(txtid, 'Date: %d-%d-%d %d:%d:%4.2f\r\n', clock);
    fprintf(txtid, 'SubjID: %s\r\n', pp.subjID);
    fprintf(txtid, 'CondNm: %s\r\n', pp.cond);
    fprintf(txtid, 'finTrialID: \n');
    fprintf(txtid, '%d;', finTrialID);
    
    fclose(txtid);
    
    ismade = 1; 
    if ismade; sprintf('%s is made successfully', pp.saveTxtname);end
end % end of sessionSucceed()
    
