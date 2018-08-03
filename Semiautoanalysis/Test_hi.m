clear
clc
Neuronlist={'Jun14IR3f'};
for neuron_count=1:length(Neuronlist)

file_High=sprintf([Neuronlist{neuron_count},'-HighRn.mat']);

% Section 3: analyze High_RN protocol. 
cd(Neuronlist{neuron_count});
mkdir(['../',Neuronlist{neuron_count},'_figs']);
hi_mahp=0;
hi_sahp=0;
%save AP count for -100:20:120 pA
High_apcount=zeros(12,2);
High_apcount(:,1)=-100:20:120;
%%%%
if exist(file_High, 'file')
load(file_High);
% the purpose for analyzing High/High RN is to compute ISI, ins_freq, delay
% and AP num. In this case, we also ignore the doublets since during
% doublets, the system is not stably oscillating. 
mat_High=matfile(file_High);
Name_High=whos(mat_High);
if(isempty(Name_High)~=1)
clear mat_High
N_High=length(Name_High)/2;
Trace_temp=eval(Name_High(2*1-1).name);
l_High=length(Trace_temp(:,1));
T_High=table([1:N_High]',zeros(N_High,l_High),zeros(N_High,l_High),'VariableNames',{'Trace','t','v'});
Cur=zeros(N_High,1);
str_temp=Name_High(1).name;
pos_underscore=strfind(str_temp,'_');
hiname_template=str_temp(1:pos_underscore(3));
f_High = fopen(['../',Neuronlist{neuron_count},'_High.txt'], 'a+');
fprintf(f_High,'Current\tNum_AP\tdelay\tMean_ISI\tstd_ISI\tCV_ISI\tFreq_ISI\tfreq_count\tfAHP\tmAHP\tsAHP\tspktime\n');   
num_aplevel=0;
min_v=100;
for i_hi=1:N_High
TraceHigh=eval([hiname_template,num2str(i_hi),'_1']);
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
Cur_temp=ceil(mean(Trace_cur(5000:10000,2))*10^12/20)*20;
v_ap=TraceHigh(:,2).*1000;
if((max(v_ap)>0)&&(Cur_temp>0))
    num_aplevel=num_aplevel+1;
end
min_v=min(min_v,min(v_ap));
end
pos_ap=0;
Cur_ap=zeros(num_aplevel,1);
Freq_ISI=zeros(num_aplevel,1);
Freq_count=zeros(num_aplevel,1);
for i_hi=1:N_High
TraceHigh=eval([hiname_template,num2str(i_hi),'_1']);
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
%How to get current: isolate the current trace and take the average between
%5000 and 10000ms, round up to eliminate the single digit
v_ap=TraceHigh(:,2).*1000;
t_ap=TraceHigh(:,1).*1000;
Cur(i_hi)=ceil(mean(Trace_cur(5000:10000,2))*10^12/20)*20;
if((max(v_ap)>0)&&(Cur(i_hi)>0))
    pos_ap=pos_ap+1;
    l=length(v_ap);
    dv_ap=get_dVdt(t_ap,v_ap);
    d2v_ap=get_dVdt(t_ap(2:end-1),dv_ap(2:end));
    t_bound=t_ap(end);
    dv_thre=3;
    [delay APnum ISI freq f_ahp m_ahp s_ahp]=hilo_ana(Neuronlist{neuron_count},1,Cur(i_hi),t_ap,v_ap,dv_ap,d2v_ap,dv_thre);
    hi_mahp=hi_mahp+mean(m_ahp);
    hi_sahp=hi_sahp+mean(s_ahp);
    Cur_ap(pos_ap)=Cur(i_hi);
    Freq_ISI(pos_ap)=mean(1000./ISI);
    Freq_count(pos_ap)=APnum/2;
    %save apcount for HighRN to High_apcount for later file output!
    highid_temp=find(High_apcount(:,1)==Cur(i_hi));
    High_apcount(highid_temp,2)=APnum;
    % save complete, if no AP then APnum=0
    fprintf(f_High,'%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t',Cur(i_hi),APnum,delay,mean(ISI),std(ISI),std(ISI)/mean(ISI),Freq_ISI(pos_ap),Freq_count(pos_ap),f_ahp,mean(m_ahp),mean(s_ahp));
    % attach the spiking times to each line!
    spkt=zeros(APnum,1);
    spkt(1)=delay;
    for apnumtemp=2:APnum
        spkt(apnumtemp)=spkt(apnumtemp-1)+ISI(apnumtemp-1);
    end
    for apnumtemp=1:APnum
        fprintf(f_High,'%f\t',spkt(apnumtemp));
    end
    fprintf(f_High,'\n');
    %
    subplot(num_aplevel+1,2,2*pos_ap-1)
    plot(t_ap,v_ap)
    xlim([0 3000]);
    ylim([min(v_ap)-1 max(v_ap)+1]);
    xlabel('Time (ms)','FontSize', 5);
    ylabel('Vm (mV)','FontSize', 5);
    title([num2str(Cur(i_hi)),'PA'],'FontSize', 5);
    subplot(num_aplevel+1,2,2*pos_ap)
    plot(ISI,'o')
    ylim([0 max(ISI)+20])
    title([num2str(Cur(i_hi)),'PA'],'FontSize', 5);
    xlabel('Spks Count','FontSize', 5);
    ylabel('ISI (ms)','FontSize', 5);
else    
    subplot(num_aplevel+1,2,2*num_aplevel+1)
    plot(t_ap,v_ap);
    xlim([0 3000]);
    ylim([min_v-2 max(v_ap)+1]);
    xlabel('Time (ms)','FontSize', 5);
    ylabel('Vm (mV)','FontSize', 5);
    hold on;   
    
end
end

subplot(num_aplevel+1,2,2*(num_aplevel+1))
    plot(Cur_ap,Freq_ISI,'-ro')
    hold on;
    plot(Cur_ap,Freq_count,'-bo')
    title('F-I Red:ISI; Blue:spkcount','FontSize', 5);
    xlabel('Current (pA)','FontSize', 5);
    ylabel('Freq (Hz)','FontSize', 5);
fclose(f_High);

cd(['../',Neuronlist{neuron_count},'_figs']);
saveas(gcf,[Neuronlist{neuron_count},'_High_RN.png']);
saveas(gcf,[Neuronlist{neuron_count},'_High_RN.fig']);
close;
% % Section 3.5 Compute sag
% i_hi=1;
% TraceHigh=eval([hiname_template,num2str(i_hi),'_1']);
% CurHigh=eval([hiname_template,num2str(i_hi),'_2']);
% Cur_sag=roundn(mean(CurHigh(5000:10000,2))*10^12,1);
% v_ap=TraceHigh(:,2).*1000;
% t_ap=TraceHigh(:,1).*1000;
% [sag_v id_sag]=min(v_ap);
% rnfit_start=21439;%floor(0.45*l_rn);
% rnfit_end=23949;%floor(0.55*l_rn);
% stable_v=mean(v_ap(rnfit_start:rnfit_end));
% sag=sag_v-stable_v;
% figure()
% plot(t_ap,v_ap);
% hold on;
% plot(t_ap(id_sag),v_ap(id_sag),'ro');
% text(t_ap(id_sag),v_ap(id_sag)+20,'This is the sag point');
% title(['Sag is computed at ',num2str(Cur_sag),'pA'])
% xlabel('Time (ms)');
% ylabel('Vm (mV)');
% saveas(gcf,[Neuronlist{neuron_count},'sag.png']);
% close;
hi_mahp=hi_mahp/num_aplevel;
hi_sahp=hi_sahp/num_aplevel;
else
    figure()
    title('HighRn empty')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'HighRn empty.png');
    close;
    hi_mahp=NaN;
    hi_sahp=NaN;
end
else
    figure()
    title('No HighRn protocol found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No HighRn.png');
    close;
    hi_mahp=NaN;
    hi_sahp=NaN;
end
cd ..
end


