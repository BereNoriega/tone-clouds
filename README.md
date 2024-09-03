# tone-clouds
Python´s translation to the tone clouds Matlab code

Original Matlab code:

%% makes CORRECTED CONTINUOUS repeating tone cloud stimuli for tapping study
%
% corrected because grid_play should always be ones according to agus etal,
% since td x sd only specify the grid, and every cell should contain a tone
%
% UPDATE: 6 oct 2022... not possible to do tone clouds as above bc
% intervals need to be multiple of grid, and if grid is fine TC becomes
% very dense/impossible for gil, and if grid is wide, tempos get pushed to
% extremes that gil also can't do. CONCLUSION: stick with original stimuli
% with p tone per griddt
%
% NOTE: (Sept 2, 2024)
% - we MAY want to change the griddt system... I created this in order to
% accommodate being multiples of the tempos that the monkeys would be
% tapping, but here we will probably just ask humans to tap @ 500 ms (2hz).
% So instead of working in multiples of griddt, we can probably just define
% the grid directly based on the desired N tones per second.
%

% some params
clear all; clc;
setID = 'a';
outdir = ['testbere/set_',setID];
mkdir(outdir)
ntokenreps = 10; % how many times the target repeats in a trial
tonedur = .05;
griddt = 0.025; % time resolution of our grid
gridfs = 1/griddt;
minf = 200;
maxf = 12800;
noct = 6;
spectral_density = [4]; % tones per octave
temporal_density = [4]; % tones per second
temp_dens_ptonepergriddt = temporal_density./gridfs; % probability of a tone in each griddt
tempi = 0.5; % ms, the inter-repeat interval of tokens 
tokendur = tempi;

fs = 44100; % generate sounds with this fs
fs48 = 48828; % convert to this fs for TDT equipment
maxf48k = fs48/2;
trialseed = 888;
mintokenstart = 1; 
maxtokenstart = 1.5; % first token starts between 1-1.5s into (random) tone cloud
rng(trialseed);
ntrialspercond = 60;
stimdb = 85;
targetrms = ( 10.^(-94./20) ).*(10.^(stimdb./20));
targetrms_wav = 0.1; %0.4;

% compute hann envelope (for creating a smooth on and off ramp to sounds)
hanndur = 0.005; % s
thann = 1/fs:1/fs:hanndur;
ramp = cos(2*pi*0.5/hanndur.*thann);
ttone = 1/fs:1/fs:tonedur;
env = ones(size(ttone));
env(1:length(ramp)) = fliplr(ramp);
env(end-length(ramp)+1:end) = ramp;
env = (env+1)./2;



%% now loop through conditions and trials to generate
trialct = 0; sp_dens = []; tp_dens = []; trialnr = []; tempo = []; play_grid = {}; latency_grid = {}; 

% to save params for manipulations/variations
g_play = {}; g_latency = {}; g_freq = {}; t_play = {}; t_latency = {}; t_freq = {}; t_startid = {}; t_ioi = {};
g_freqlist = {}; g_seqdur = {}; g_setIDs = []; g_specdens = []; g_tempdens = [];
for s = 1:length(spectral_density)
    specdens = spectral_density(s);
    freq_list = 10.^linspace(log10(minf),log10(maxf),specdens*noct+1); % bins edges of freq
    for t = 1:length(temporal_density)
        tempdens = temp_dens_ptonepergriddt(t);
        for n = 1:ntrialspercond
            for N = 1:length(tempi)
                trialct = trialct+1;
                startt = round((maxtokenstart-mintokenstart)/griddt*rand) + mintokenstart/griddt; % time of first token start in units of n griddt
                
                % from inside if statement (if statement added 6nov2022)
                startt_ms = round(startt*griddt*1000); % should be integer...
                ioi = round(tempi(N)/griddt); % everything in units of griddt
                seqdur = startt+ntokenreps*ioi;
                tokenstartid = startt:ioi:seqdur-1;

                tseq = 1/fs:1/fs:seqdur*griddt; % time (s)
                grid_t = 0:seqdur-1; % time in ngriddt grid
                soundmat = zeros(length(freq_list)-1, round((seqdur+2)*griddt*fs)); % empty time vector for the whole sequence +2 griddt for tail end of thigns in last grid

                % first, rand background with correct temporal and spectral density

                % which cells in grid will contain a tone to play?
                tone_play = rand(length(freq_list)-1,seqdur); % generating freq x delta matrix filled with random numbers 0:1
                tone_play = tone_play < tempdens; % logical 1 if rand < p(tone per grid) based on temporal density
                %tone_play = ones(length(freq_list)-1,seqdur);

                % similarly, grid of latencies
                tone_latency = rand(length(freq_list)-1,seqdur).*griddt;
                tone_latency = tone_latency.*tone_play; % only care to keep latency where token_play=1

                % new variables for latency and play
                grid_play = tone_play;
                grid_latency = tone_latency;
                grid_freq = zeros(size(grid_play));
                for f = 1:size(soundmat,1) % ## create EXACT repeats! 24aug2022
                    fchanbins = [freq_list(f),freq_list(f+1)];
                    for d = 1:size(grid_play,2) % grid_t
                        if grid_play(f,d)==1
                            ftone = 10^(diff(log10(fchanbins))*rand+log10(fchanbins(1))); % logspaced rand
                            %ftone = diff(fchanbins)*rand+fchanbins(1);
                            grid_freq(f,d) = ftone;
                        end
                    end
                end 

                % ### save for later
                g_play{trialct} = grid_play;
                g_latency{trialct} = grid_latency;
                g_freq{trialct} = grid_freq;

                if ~strcmp(setID,'d')
                    % above will generate random tone cloud for duration of trial
                    % now select a section starting from startt to repeat
                    token_play = grid_play(:,startt:startt+tokendur(N)/griddt-1);
                    token_latency = grid_latency(:,startt:startt+tokendur(N)/griddt-1);
                    token_freq = grid_freq(:,startt:startt+tokendur(N)/griddt-1);

                    % ### save TOKEN for later
                    t_play{trialct} = token_play;
                    t_latency{trialct} = token_latency;
                    t_freq{trialct} = token_freq;

                    % ### save some other things for later
                    t_startid{trialct} = tokenstartid;
                    t_ioi{trialct} = ioi;
                    g_freqlist{trialct} = freq_list;
                    g_seqdur{trialct} = seqdur;
                    g_setIDs(trialct) = n;
                    g_specdens(trialct) = specdens;
                    g_tempdens(trialct) = temporal_density(t);

                    % check that there are 10 token starts
                    if ~isequal(length(tokenstartid),ntokenreps)
                        keyboard
                    end

                    % now replace sections of tone_play and tone_latency ## AND
                    % grid_freq
                    for ii = 1:ntokenreps
                        grid_play(:,tokenstartid(ii):tokenstartid(ii)+size(token_play,2)-1) = token_play;
                        grid_latency(:,tokenstartid(ii):tokenstartid(ii)+size(token_latency,2)-1) = token_latency;
                        grid_freq(:,tokenstartid(ii):tokenstartid(ii)+size(token_freq,2)-1) = token_freq;
                    end

                    %                 % visualise grid to check
                    %                 close all; figure; hold on;
                    %                 for gg = 0:seqdur % draw vertical lines
                    %                     plot([gg,gg]./gridfs,[freq_list(1),freq_list(end)],'linewidth',1.5,'Color',[.8,.1,.1])
                    %                 end
                    %                 for ff = 1:length(freq_list) % draw horizontal lines
                    %                     plot([0,seqdur]./gridfs,[freq_list(ff),freq_list(ff)],'linewidth',1.5,'Color',[.8,.1,.1])
                    %                 end
                    %                 for f = 1:length(freq_list)-1 % freq chans
                    %                     fchanbins = [freq_list(f),freq_list(f+1)];
                    %                     for d = 1:length(grid_t) % grid_t
                    %                         if grid_play(f,d)==1
                    %                             start_time = grid_t(d)*griddt + grid_latency(f,d);
                    %                             ftone = grid_freq(f,d);
                    %                             plot([start_time,start_time+.05],[ftone,ftone],'k','linewidth',2)
                    %                         end
                    %                     end
                    %                 end
                    %                 pause
                    %                 continue
                    %                 % looks good!

                else % control: just keep it random for whole trial
                    'control: not repeating token'
                end


                % GENERATE SIGNALS (always REP)
                for f = 1:size(soundmat,1) % freq chans
                    fchanbins = [freq_list(f),freq_list(f+1)];
                    for d = 1:length(grid_t) % grid_t
                        if grid_play(f,d)==1
                            start_time = grid_t(d)*griddt + grid_latency(f,d);
                            start_samp = start_time * fs;
                            if round(start_samp)==0
                                start_samp = 1;
                            end
                            ftone = grid_freq(f,d); %diff(fchanbins)*rand+fchanbins(1); % random freq inside this freq bin
                            amp = sin(2*pi*ftone.*ttone);
                            amp = amp.*env; % hanning env already calculated outside
                            if round(start_samp)+length(amp)-1>size(soundmat,2)
                                keyboard
                            end
                            soundmat(f,round(start_samp):round(start_samp)+length(amp)-1) = soundmat(f,round(start_samp):round(start_samp)+length(amp)-1) + amp;
                        end
                    end
                end
                stim = sum(soundmat,1); % summing across all freqs
                %stim = stim./rms(stim).*targetrms;

                % resample
                xorig = 1/fs:1/fs:length(stim)/fs;
                xnew = 1/fs48:1/fs48:length(stim)/fs;
                y48 = interp1(xorig,stim,xnew);
                tmp = find(isnan(y48));
                if length(tmp)>1
                    keyboard
                else
                    y48(tmp) = [];
                end
                y48 = y48./rms(y48).*targetrms_wav;

                % write to file
                % format: 08282022_550_t2s4_seed.wav
                sd = ['0',num2str(specdens)]; sd = sd(end-1:end);
                td = ['0',num2str(temporal_density(t))]; td = td(end-1:end);
                tr = ['0',num2str(n)]; tr = tr(end-1:end);
                tempo = num2str(tempi(N)*1000); % in ms
                exptseed = num2str(trialseed);
                fname = [outdir,'/',tr,'_',tempo,'_t',td,'s',sd,'_',exptseed,'_',num2str(startt_ms),'.wav'];
                %fname = [outdir,'/set',setID,'_seed',exptseed,'_',sd,'fpero_',td,'tpers_ioi',tempo,'ms_trial',tr,'.wav'];
                audiowrite(fname,y48,fs48);

                % mat file also (to be used by synapse at correct rms for 85 db spl)
                y48 = y48./rms(y48).*targetrms;
                fname = [outdir,'/',tr,'_',tempo,'_t',td,'s',sd,'_',exptseed,'_',num2str(startt_ms),'.mat'];
                save(fname,'y48');

                % also generate onsets file
                onsets_ms = round(grid_t.*griddt*1000);
                onset_is_played = zeros(size(onsets_ms));
                for gg = 1:length(tokenstartid) % where tokens start in units of n griddt
                    onset_is_played(tokenstartid(gg)) = 1;
                end
                %fname = [outdir,'/set',setID,'_seed',exptseed,'_',sd,'fpero_',td,'tpers_ioi',tempo,'ms_trial',tr];
                fname = [outdir,'/',tr,'_',tempo,'_t',td,'s',sd,'_',exptseed,'_',num2str(startt_ms)];
                fid = fopen(sprintf('%s.txt',fname),'w');
                for ii = 1:length(onsets_ms)
                    fprintf(fid,'%d,%d\n',onsets_ms(ii),onset_is_played(ii));
                end
                fclose(fid);

                % SAVE THINGS TO RECREATE LATER
                

            end
%             close all; figure;
%             subplot(1,2,1); imagesc([1:round(seqdur*griddt*fs)]./fs,1:size(soundmat,1),soundmat);
%             subplot(1,2,2); plot([1:round(seqdur*griddt*fs)]./fs,stim);
%             keyboard
            
        end
    end
end

save([outdir,'/t',td,'s',sd,'_',exptseed,'_',num2str(startt_ms),'.mat'],'g_play','g_latency','g_freq','t_freq','t_latency','t_play', ...
    'startt','fs48','trialseed','t_ioi','t_startid','g_freqlist','g_seqdur','g_setIDs','g_specdens','g_tempdens');
