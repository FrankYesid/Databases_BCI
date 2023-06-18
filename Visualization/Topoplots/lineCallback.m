function lineCallback(snr,evt,database)
if database == 1
    load('BCICIV_1\electrodesBCICIV1.mat')
elseif database == 2
    load('BCICIV_2a\electrodesBCICIV2a.mat')
    electrodes = Channels;
    pos = elec_pos;
elseif database == 3
    load('Gigasc\posiciones.mat')
    electrodes = M1.lab;
    pos = M1.xy;
elseif database == 4
    load('BCICIIIIVa\electrodesBCICIIIIVa.mat')
end
if snr.Color == [0 0 0]
    snr.Color = 'red';
    b = strfind(electrodes,snr.String,'ForceCellOutput',true);
    for i = 1:length(electrodes)
        if b{i} == 1
            a = i;
        end
    end
    snr.String = num2str(a);
else
    snr.Color = 'black';
    snr.String = electrodes{str2num(snr.String)};
end