% % EEG preprocessing with EEGLAB
% % author:finch
% % this pipline follows these steps:
%1. Importing raw EEG data and behavioral data, rejectd events according to your standard
%2. Resampling (down to 500 HZ)
%3. Filtering (high and low pass)
%4. Detection of bad channels
%5. Independent Components Analysis (ICA)
%6. Re-reference (bilateral earlobe or average reference)
%7. Baseline correction
%8. Detection of high amplitude noise
%9. Epoching of the continuous data

clear all;
close all;
clc;
eeglab;
%% define parameters
cd 'your dataset path'
sub_name='S1'
marks={'66','88','55','99'};
high_pass=0.01;
low_pass=30;
reref='earlobes'  % bilateral earlobes or average re-reference
first_epoch=[-0.8,1.5];
epoch=[-0.2,0.8];
baseline=[-0.2,0];
bad_range=[-70,70];
save_path='';

%% load data, merge behavioral data with EEG 
sub_name=['event_',sub_name];
EEG = pop_loadcnt() % dialogue open file
% load electrodes file
EEG = pop_chanedit(EEG, 'lookup','C:\Program Files\MATLAB\R2014a\toolbox\eeglab13_4_4b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );

p=1;
while p<length(EEG.event)+1
        if strcmpi(EEG.event(p).type,'boundary') 
           EEG = pop_editeventvals(EEG,'delete',p); 
           fprintf('boundary at %d \n has been deleted', p)           
           p=p-1;                                    
        end;
        p=p+1;
end;

event = load('D:\sub1.txt');

j=1;
while (j<length(EEG.event)+1)
     k = EEG.event(j).type;  
        if k~=num2str(event(j,6))
            fprintf('########### mismatched at ###########\n ID: %d \n', j) %         
            return;  
        end;
        j=j+1;
end;
fprintf('number of event: %d\n              J: %d\n', length(EEG.event), j-1); % if equal, all events are consistent

% create new event fields with initiative values of '0'
EEG = pop_editeventfield(EEG, 'rt', 0);  
EEG = pop_editeventfield(EEG, 'acc', 0); 
EEG = pop_editeventfield(EEG, 'item', 0); 
EEG = pop_editeventfield(EEG, 'sd', 0); 

for i=1:length(EEG.event)
    EEG.event(i).rt  = event(i,7);  % add rt value which is in the 7th column of the txt. file
    EEG.event(i).acc = event(i,8); 
    EEG.event(i).item = event(i,2);
    EEG.event(i).sd = event(i,9);
end;

% delete trials 
p = 1; 
while p<length(EEG.event)+1
     if EEG.event(p).rt < 400 ||  EEG.event(p).rt > 1500 || EEG.event(p).acc ~= 0
         EEG = pop_editeventvals(EEG,'delete',p);
         fprintf('trials deleted at %d \n', p) 
         p=p-1;                    % because events were deleted
      end
      p = p+1;
end
pop_saveset( EEG, 'filename',sub_name, 'filepath', save_path);

%% down sampling to 500hz
EEG = pop_resample( EEG, 500);
EEG = eeg_checkset( EEG );

%% high-pass and low-pass filter
sub_name=['filt_',sub_name];
EEG = pop_eegfiltnew(EEG, [], low_pass, [], 0, [], 1);
EEG = pop_eegfiltnew(EEG, [], high_pass, [], true, [], 1);
EEG = eeg_checkset( EEG );
pop_saveset( EEG, 'filename',sub_name, 'filepath', save_path);

%% plot data to check bad channel visually
pop_eegplot(EEG)

%% bad channel interpolation
sub_name=['badChannel_',sub_name];
badchannel={'F7'};
chanind=region_chan(EEG,badchannel);
EEG=eeg_interp(EEG,find(chanind==1));
pop_saveset( EEG, 'filename',sub_name, 'filepath', save_path);

%% ICA
sub_name=['ICA_',sub_name];
EEG = pop_select( EEG, 'nochannel',{'HEO' 'VEO'});
%compute rank for ICA
if isfield(EEG.etc,'clean_channel_mask')
    dataRank = min([rank(double(EEG.data(:,:,1)')) sum(EEG.etc.clean_channel_mask)]);
else
    dataRank = rank(double(EEG.data(:,:,1)'));
end
EEG = pop_epoch( EEG, marks2, first_epoch, 'newname', 'EGI new epochs', 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',dataRank);
EEG = eeg_checkset( EEG );
pop_ADJUST_interface(ALLEEG, EEG, CURRENTSET);
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
input('ICA selection finished?');
reject_component=find(EEG.reject.gcompreject==1);
show=sprintf('ica%d  ',reject_component);
display('ICA you selected£º');
display(show);
com='';
    % check the result,if not good enough, select ICA again
while isempty(com)
    [EEG,com] = pop_subcomp( EEG, reject_component, 1);
    if isempty(com)
        pop_ADJUST_interface(ALLEEG, EEG, CURRENTSET);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        input('ICA selection finished?');
        reject_component=find(EEG.reject.gcompreject==1);
        show=sprintf('ica%d  ',reject_component);         
        display('ICA you selected£º');
        display(show);
    end
end
pop_saveset( EEG, 'filename',sub_name, 'filepath', save_path);

%% re-reference
if reref == 'earlobes'
    % re-reference to average of both earlobes
    fprintf('earlobes')
    chan_idx=strcmpi('M2',{EEG.chanlocs.labels});
    EEG.data=bsxfun(@minus,EEG.data,0.5*EEG.data(chan_idx,:,:));
    EEG = pop_select( EEG, 'nochannel',{'M1' 'M2'});
elseif reref == 'average'
    % average reference
    fprintf('earlobesno')
    EEG = pop_select( EEG, 'nochannel',{'M1' 'M2'});
    EEG.nbchan = EEG.nbchan + 1;
    EEG.data(end+1,:,:) = zeros(1, EEG.pnts, EEG.trials);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select(EEG,'nochannel',{'initialReference'});
end

%% baseline correction and epoch
for i=1:size(marks,2)
    EEGn = pop_epoch( EEG, {mark{i}}, epoch, 'newname', ['_',num2str(marks{i})], 'epochinfo', 'yes');
    EEGn = pop_rmbase(EEGn, baseline);
    EEGn = pop_eegthresh(EEGn,1,1:62 ,bad_range(1),bad_range(2),epoch(1),epoch(2),0,0);
    EEGn = eeg_checkset( EEGn );
    reject_number=find(EEGn.reject.rejthresh==1);
    EEGn = eeg_rejsuperpose( EEGn, 1, 1, 1, 1, 1, 1, 1, 1);
    EEGn = pop_rejepoch( EEGn, reject_number ,0);
    pop_saveset( EEGn, 'filename',[sub_name,'_',num2str(file_name{i}),'.set'],'filepath',save_path);
    clear EEGn
end
