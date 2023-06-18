clear; clc
%% path Database
SUBJECTS_DIR = 'E:\Database_2\files';
SUBJECTS_DIR2 = 'F:\Physionet';
SUBJECTS_DIR3 = 'D:\Physionet';
% Organizar los datos.
COHORT   = 'S';
SUBJECTS = dir([SUBJECTS_DIR3 filesep COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,1:end)';
SS       = 1:numel(SUBJECTS);
% SS([88,92,100]) = [];
SUBJECTS_DIR_SAVE = 'D:\Physionet';

for s = 88:106
    load([SUBJECTS_DIR3 filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'])
    clear cv
    if size(y,2)>1
        y_ = [];
        for a = 1:size(y,2)
            y_ = [y_;y(:,a)];
        end
        clear y
        y = y_;
    end  
    cv = cvpartition(y,'KFold',10);
    save([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg' filesep 'cv.mat'],...
        'cv')
%     save([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'],...
%             'X','y','fs','channels','seg_start','seg_end')
end