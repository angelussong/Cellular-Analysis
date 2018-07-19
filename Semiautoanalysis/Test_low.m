clear
clc
Neuronlist={'Apr18IR4a'};


neuron_count=1;

file_low=sprintf([Neuronlist{neuron_count},'-LowRn.mat']);

% Section 3: analyze low_RN protocol. 
cd(Neuronlist{neuron_count});
load(file_low);
% the purpose for analyzing High/Low RN is to compute ISI, ins_freq, delay
% and AP num. In this case, we also ignore the doublets since during
% doublets, the system is not stably oscillating. 
mat_low=matfile(file_low);
Name_low=whos(mat_low);
clear mat_low
N_low=length(Name_low)/2;
Trace_temp=eval(Name_low(2*1-1).name);
l_low=length(Trace_temp(:,1));
T_low=table([1:N_low]',zeros(N_low,l_low),zeros(N_low,l_low),'VariableNames',{'Trace','t','v'});
Cur=zeros(N_low,1);
str_temp=Name_low(1).name;
pos_underscore=strfind(str_temp,'_');
hiname_template=str_temp(1:pos_underscore(3));
f_low = fopen(['../',Neuronlist{neuron_count},'_low.txt'], 'a+');
fprintf(f_low,'Current\tNum_AP\tdelay\tMean_ISI\tstd_ISI\tCV_ISI\tFreq_ISI\tfreq_count\tmAHP\tsAHP\n');   
num_aplevel=0;
min_v=100;
for i_hi=1:N_low
TraceLow=eval([hiname_template,num2str(i_hi),'_1']);
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
Cur_temp=roundn(mean(Trace_cur(5000:10000,2))*10^12,1);
v_ap=TraceLow(:,2).*1000;
if((max(v_ap)>0)&&(Cur_temp>0))
    num_aplevel=num_aplevel+1;
end
min_v=min(min_v,min(v_ap));
end

pos_ap=0;
Cur_ap=zeros(num_aplevel,1);
Freq_ISI=zeros(num_aplevel,1);
Freq_count=zeros(num_aplevel,1);
for i_hi=1:N_low
Tracelow=eval([hiname_template,num2str(i_hi),'_1']);
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
%How to get current: isolate the current trace and take the average between
%5000 and 10000ms, round up to eliminate the single digit
v_ap=Tracelow(:,2).*1000;
t_ap=Tracelow(:,1).*1000;
Cur(i_hi)=roundn(mean(Trace_cur(5000:10000,2))*10^12,1);
if((max(v_ap)>0)&&(Cur(i_hi)>0))
    pos_ap=pos_ap+1;
    l=length(v_ap);
    dv_ap=get_dVdt(t_ap,v_ap);
    d2v_ap=get_dVdt(t_ap(2:end-1),dv_ap(2:end));
    t_bound=t_ap(end);
    dv_thre=3;
    
    [delay APnum ISI freq m_ahp s_ahp]=hilo_ana(Neuronlist{neuron_count},2,Cur(i_hi),t_ap,v_ap,dv_ap,d2v_ap,dv_thre);
    Cur_ap(pos_ap)=Cur(i_hi);
    Freq_ISI(pos_ap)=mean(1000./ISI);
    Freq_count(pos_ap)=APnum/2;
    
    fprintf(f_low,'%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Cur,APnum,delay,mean(ISI),std(ISI),std(ISI)/mean(ISI),Freq_ISI(pos_ap),Freq_count(pos_ap),mean(m_ahp),mean(s_ahp));
    subplot(num_aplevel+1,2,2*pos_ap-1)
    plot(t_ap,v_ap)
    xlim([0 3000]);
    ylim([min(v_ap)-1 max(v_ap)+1]);
    title([num2str(Cur(i_hi)),'PA']);
    subplot(num_aplevel+1,2,2*pos_ap)
    plot(ISI,'o')
else    
    subplot(num_aplevel+1,2,2*num_aplevel+1)
    plot(t_ap,v_ap);
    xlim([0 3000]);
    ylim([min_v-2 max(v_ap)+1]);
    hold on;    
end
end

subplot(num_aplevel+1,2,2*(num_aplevel+1))
    plot(Cur_ap,Freq_ISI,'-ro')
    hold on;
    plot(Cur_ap,Freq_count,'-bo')
    legend('F-I from ISI','F-I from AP_num');
%fclose(f_low);
cd ..



