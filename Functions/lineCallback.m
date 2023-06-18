function lineCallback(snr,~)
load('BCICIV_1\electrodesBCICIV1.mat')
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