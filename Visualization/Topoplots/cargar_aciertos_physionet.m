clear; clc
SS = 1:109;SS([88,92,100,104]) = [];
a =1;
for s = SS
    for time = [5,1,15,2]
        if time== 5
            if s<10
                load(['C:\Users\frany\Desktop\escritorio\Variando_freqs_time\Data_05seg_6cspS00',num2str(s),'.mat'])
            elseif s<100
                load(['C:\Users\frany\Desktop\escritorio\Variando_freqs_time\Data_05seg_6cspS0',num2str(s),'.mat'])
            else
                load(['C:\Users\frany\Desktop\escritorio\Variando_freqs_time\Data_05seg_6cspS',num2str(s),'.mat'])
            end
        else
            if s<10
                load(['C:\Users\frany\Desktop\escritorio\Variando_freqs_time\Data_',num2str(time),'time_2frec_sS00',num2str(s),'.mat'])
                
            elseif s<100
                load(['C:\Users\frany\Desktop\escritorio\Variando_freqs_time\Data_',num2str(time),'time_2frec_sS0',num2str(s),'.mat'])
            else
                load(['C:\Users\frany\Desktop\escritorio\Variando_freqs_time\Data_',num2str(time),'time_2frec_sS',num2str(s),'.mat'])
            end
        end
        [acc_,pos] = max(Data.Result.mean_test_acc);
        std = Data.Result.std_test_acc(pos);
        acc = [acc_,std];
        if time == 5
            Acc{1,1}(a,:) = acc;
        elseif time == 1
            Acc{1,2}(a,:) = acc;
        elseif time == 15
            Acc{1,3}(a,:) = acc;
        elseif time ==2
            Acc{1,4}(a,:) = acc;
        end
    end
    a = a+1;
end
for time = 1:4
    Dat{time} = Acc{1,time}(:,1).*100;
end