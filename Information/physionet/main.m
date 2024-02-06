%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data physionet
% Information in Readme.txt
% link : https://physionet.org/content/eegmmidb/1.0.0/
% Data
% Each annotation:
%    T0 corresponds to rest
%    T1 corresponds to onset of motion (real or imagined) of
%       the left fist (in runs 3, 4, 7, 8, 11, and 12)
%        both fists (in runs 5, 6, 9, 10, 13, and 14)
%    T2 corresponds to onset of motion (real or imagined) of
%        the right fist (in runs 3, 4, 7, 8, 11, and 12)
%        both feet (in runs 5, 6, 9, 10, 13, and 14)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Subjects 88,92,100 with time different
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F. Y. Zapata C.
% Version 1.0
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc
%% path toolbox
% eeglab
% run('D:\Toolbox\eeglab14_1_2b\eeglab.m')
% biosig
run('D:\Dropbox\ERD\Toolbox\Biosig_ERD\biosig_installer.m')
close;

%% path Database
COHORT = 'S';
SUBJECTS_DIR = 'F:\Physionet_MI\files';
SUBJECTS = dir([SUBJECTS_DIR filesep COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,1:end-1)';
SUBJECTS([88,92,100,104]) = []; % Sujetos que no vamos a tener encuenta.
SS       = 1:numel(SUBJECTS);   % numeros de los sujetos seleccionados.

%% path SAVE
SUBJECTS_DIR_SAVE = 'D:\Physionet';

%% Parameters
t        = [0,8.1]; % time interval for analysis.
channels = 1:64;    % select channels.
T        = 2;       % task 1.rest, 2.left fist or rigth fist, 3.both fists or both feet
task     = 'MI';    % 'MI'. Motor Imagery, 'ME'. Motor Execution
task_r   = 1;       % if is task 1 Baseline, eyes open, task 2 Baseline, eyes closed.
t_r      = [0,2];   % time interval for analysis in resting.
tri      = 2;       % tri in 1 in trials, and tri in 2 continuo.

%% Load data
for s = SS
    % Tasks of subject
    SUBJECTS2 = dir([SUBJECTS_DIR filesep SUBJECTS{s} filesep COHORT '*']);
    SUBJECTS2 = struct2cell(SUBJECTS2);
    SUBJECTS2 = SUBJECTS2(1,1:2:end)';
    % clases de la base de datos.
    classes = [1,2];
    % se realiza la selección de que tipo de datos se van a seleccionar.
    if T == 1     % T0 - resting state
        SUBJECTS2 = SUBJECTS2(task_r);
    elseif T == 2 % T1 - (real or imagined) left or rigth fist
        if strcmp(task,'MI')
            SUBJECTS2 = SUBJECTS2([4,8,12]);
        elseif strcmp(task,'ME')
            SUBJECTS2 = SUBJECTS2([3,7,11]);
        end
    elseif T == 3 % T1 - (real or imagined) both fists or both feet
        if strcmp(task,'MI')
            SUBJECTS2 = SUBJECTS2([6,10,14]);
        elseif strcmp(task,'ME')
            SUBJECTS2 = SUBJECTS2([5,9,13]);
        end
    end
    data_ = cell(numel(SUBJECTS2),1);
    y_ = cell(numel(SUBJECTS2),1);
    % selección de la tarea realizada.
    if T ~= 1 
        for s1 = 1:numel(SUBJECTS2)
            %% read of data
            % lugar donde se ubica la base de datos.
            pathname = [SUBJECTS_DIR filesep SUBJECTS{s} filesep];
            % nombre del archivo con los datos.
            filename = SUBJECTS2{s1};
            % load data(EEG), h1 anotaciones del registro
            [h1, data]= edfread([pathname,filename],'targetSignals',channels);
            % load annotation, tiempos y marcaciones.
            [~, h]   = readEDF([pathname,filename]);
            % selección de las clases, duración del evento
            if s == 104 || s == 88
                Task_label = str2double(reshape(cell2mat(cellfun(@(x) x(2),h.annotation.event,'UniformOutput',false)),[numel(h.annotation.event),1]));
                Time_duration = h.annotation.duration;
                Task_sym   = h.annotation.data;
                strArray   = h.annotation.data;
            else
                % load marked of data
                [Task_label,Time_duration,Task_sym,strArray] = Eventread(pathname,filename);
            end
            % Select trial
            fs        = h.samplerate(1);                % frecuencia de muestreo
            triallen  = round((t(2) - t(1)) * fs) + 1;  % Trial length (in samples)
            indx      = (Task_label==1) + (Task_label==2); % indice del muestreo.
            indx(end) = 0; % elimino el último trial por motivos del tiempo del trial.
            % selección de los trials en la señal completa
            tmp       = trigg(data', h.annotation.starttime(ismember(Task_label, classes) & logical(indx))*fs, round(t(1)*fs)-2*fs, round(t(2)*fs)-2*fs);
            % organiza los datos en tre dimensiones (canales,tiempo,trial)
            tmp1      = reshape(tmp, size(tmp,1),triallen, floor(length(tmp)/triallen));  %
            tmp1(isnan(tmp1))=0; % coloca los datos en cero que sean NaN.
            % organiza los trials en celdas. y en cada celda se encuentra un (canales,tiempo)
            data_trial= cell(1,size(tmp1,3));
            for trial = 1:size(tmp1,3)
                data_trial{1,trial} = tmp1(:,1:triallen,trial);
            end
            data_{s1} = data_trial;
            y_{s1}    = Task_label(logical(indx));
        end
    else
        for s1 = 1:numel(SUBJECTS2)
            %% read of data
            % lugar donde se ubica la base de datos.
            pathname = [SUBJECTS_DIR filesep SUBJECTS{s} filesep];
            % nombre del archivo con los datos.
            filename = SUBJECTS2{s1};            
            % load data(EEG), h1 anotaciones del registro
            [h1, data]= edfread([pathname,filename],'targetSignals',channels);
            % load annotation, tiempos y marcaciones.
            [~, h]    = readEDF([pathname,filename]);
            % Select trial      
            fs        = h.samplerate(1); % frecuencia de muestreo
            %h.annotation.starttime = 0:2*fs:size(data,2)-2*fs;
            %Task_label= ones(floor(size(data,2)/(2*fs)),1);
            triallen     = 2*fs;%round((t(2) - t(1)) * fs) + 1;  % Trial length (in samples)
            %indx      = ones(floor(size(data,2)/(2*fs)),1);    % indice del muestreo.
            %indx(end) = 0;              % elimino el último trial por motivos del tiempo del trial.
            % selección de los trials en la señal completa
            %tmp       = trigg(data', h.annotation.starttime(ismember(Task_label, 1) & logical(indx)), round(t_r(1)*fs), round(t_r(2)*fs));
            % organiza los datos en tre dimensiones (canales,tiempo,trial)
            %tmp1      = reshape(tmp, size(tmp,1),triallen, floor(length(tmp)/triallen));  %
            %tmp1(isnan(tmp1))=0; % coloca los datos en cero que sean NaN.
            % organiza los trials en celdas. y en cada celda se encuentra un (canales,tiempo)
            if tri == 1
                data_trial = cell(1,1); pos_t = 1;
                for trial = 1:2*fs:size(data_,1)
                    if trial+(2*fs)-1 <= length(data_)
                        data_trial{1,pos_t} = data_(trial:trial+(2*fs)-1,:)';% tmp1(:,1:triallen,trial);
                        pos_t = pos_t+1;
                    end
                end
            else
                data_{s1} = data;
            end
        end
    end
    % Datos organizados en X - Data ,y clases, fs frecuencias de muestreo
    % Seg_start punto de inicio de imaginación motora,
    % Seg_end punto de fin de imaginación motora.
    % channels contiene el numbre de los canales.
    
    if  T ~= 1
        y = cell2mat(y_);
        X = [];
        for s1 = 1:numel(SUBJECTS2)
            X  = [X,data_{s1}];
        end
        channels = h1.label;
        seg_start= 2*fs;
        seg_end  = 2+Time_duration(2)*fs;
        %Almacena los datos organizados. mkdir es para saber si esta el lugar
        %de guardado de los datos.
        if mkdir([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg'])
            save([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'],...
                'X','y','fs','channels','seg_start','seg_end')
        else
            mkdir([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg'],'newdir')
            save([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'],...
                'X','y','fs','channels','seg_start','seg_end')
        end
    else
        X = data_{s1};
        channels = h1.label;
        if mkdir([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg'])
            save([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg' filesep 'rest.mat'],...
                'X','fs','channels')
        else
            mkdir([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg'],'newdir')
            save([SUBJECTS_DIR_SAVE filesep SUBJECTS{s} filesep 'eeg' filesep 'rest.mat'],...
                'X','fs','channels')
        end
    end
end

% almacenar X signal {trial}(ch,time), y labels de clases, fs frecuencia de
% la señal,
% codes of load edf
% [hdr, record] = pop_readedf([pathname,filename]);
% [h,data] = edfread([SUBJECTS_DIR filesep SUBJECTS{s} filesep SUBJECTS2{s1}]);
% [Task_label,Time_duration,Task_sym,strArray] = Eventread([SUBJECTS_DIR filesep SUBJECTS{s} filesep], SUBJECTS2{s1});

