% clear 
% clc;

filename = 'F:\Toolbox\EEGLAB\eeglab14_1_2b\eeglab14_1_2b\sample_locs\Standard-10-20-Cap81.ced';
EEG=readlocs(filename);
EEG(2).labels = 'FPz';
EEG(8).labels = 'AFz';
EEG(17).labels = 'F';
EEG(17).labels = 'Fz';
EEG(50).labels = 'CPz';
EEG(61).labels = 'Pz';

load labels.mat
% chan = EEG;
for ch = 1:numel(labels)
%     posi = {};
    for ch2 = 1:81
       if find(strcmp(EEG(ch2).labels,M1.lab(ch)),1)
           posi = ch2;
       end
    end
%     pos = cell2mat(posi);
%     chan.lab{ch,1}  = EEG(posi).labels;
%     chan.phase(ch,1)  = EEG(posi).theta;
%     chan.phase(ch,2)  = EEG(posi).radius;
%     chan.xyz(ch,1)  = EEG(posi).X;
%     chan.xyz(ch,2)  = EEG(posi).Y;
%     chan.xyz(ch,3)  = EEG(posi).Z;
%     chan.sph(ch,1)  = EEG(posi).sph_theta;
%     chan.sph(ch,2)  = EEG(posi).sph_phi;
%     chan.sph(ch,3)  = EEG(posi).sph_radius;
%     chan.besa(ch,1) = EEG(posi).sph_theta_besa;
%     chan.besa(ch,2) = EEG(posi).sph_phi_besa;
%     chan.type{ch,1} = 'EEG';%EEG(posi).type;
    chan(ch).labels  = EEG(posi).labels;
    chan(ch).theta  = EEG(posi).theta;
    chan(ch).radius  = EEG(posi).radius;
    chan(ch).X  = EEG(posi).X;
    chan(ch).Y  = EEG(posi).Y;
    chan(ch).Z  = EEG(posi).Z;
    chan(ch).sph_theta  = EEG(posi).sph_theta;
    chan(ch).sph_phi  = EEG(posi).sph_phi;
    chan(ch).sph_radius  = EEG(posi).sph_radius;
    chan(ch).sph_theta_besa = EEG(posi).sph_theta_besa;
    chan(ch).sph_phi_besa = EEG(posi).sph_phi_besa;
    chan(ch).type = 'EEG';%EEG(posi).type;
end
clear M1
M1 = chan;
save('pos_struct_BCI_CIV_2a.mat','M1')

for  ch = 1:numel(labels)
    elec_pos(ch,1) = chan(ch).X;
    elec_pos(ch,2) = chan(ch).Y;
    elec_pos(ch,3) = chan(ch).Z; 
    Channels{ch} =  chan(ch).labels;
end

save('electrodesBCICIV2a.mat','elec_pos','Channels')