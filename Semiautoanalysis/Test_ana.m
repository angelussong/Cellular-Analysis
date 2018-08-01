clear
clc
Neuronlist={'Apr18IR4d'};

for neuron_count=1:length(Neuronlist)
% set default values to all the outputs to 0 so in the file if anything is
% 0, that protocol is bad!!
thre=NaN;
pe=NaN;
ampl=NaN;
riset=NaN;
decayt=NaN;
halft=NaN;
tau=NaN;
b1=[NaN NaN];
RN_allfit=NaN;
RN_fifit=NaN;
sag=NaN;
rheo=NaN;
f_ahp=NaN;

file_AP = sprintf([Neuronlist{neuron_count},'-AP.mat']); 
file_RN = sprintf([Neuronlist{neuron_count},'-Rn.mat']);
file_High=sprintf([Neuronlist{neuron_count},'-HighRn.mat']);
file_low=sprintf([Neuronlist{neuron_count},'-LowRn.mat']);
file_rheo = sprintf([Neuronlist{neuron_count},'-ramp.mat']); 
mkdir([Neuronlist{neuron_count},'_figs']);
cd(Neuronlist{neuron_count})


%%%% section 1 RN analysis
%Resting: linear regression-y-intersection.
if exist(file_RN, 'file')
load(file_RN);
mat_tmp=matfile(file_RN);
% This is to read in the traces for tau_fitting as base
% N_trace: # of traces for RN sweep
% N_tau: # of trace for tau_fitting!!! 
Name_temp=whos(mat_tmp);
if(isempty(Name_temp)~=1)
N_trace=length(Name_temp)/2;
N_tau=4;
exp_tv=eval(Name_temp(1+2*(N_tau-1)).name).*1000;
l_rn=length(exp_tv);
%setup a table to store data, t v the original, vf the jump fixed 
%vfa the fixed and re-aligned data
T_RN=table([1:N_trace]',zeros(N_trace,l_rn),zeros(N_trace,l_rn),zeros(N_trace,l_rn),zeros(N_trace,l_rn),'VariableNames',{'Trace','t','v','vf','vfa'});
rnfit_start=2125;%floor(0.45*l_rn);
rnfit_end=2475;%floor(0.55*l_rn);
fit_i=-40:10:40;
fit_v=zeros(N_trace,1);
fff=figure(3);
movegui(fff,[25,25])
for i=1:N_trace
    RN_tmp=eval(Name_temp(1+2*(i-1)).name);
    T_RN.t(i,:)=RN_tmp(:,1).*1000;
    T_RN.v(i,:)=RN_tmp(:,2).*1000;
    %legendInfo{i} = ['Trace# ',num2str(i),' represents ',num2str(-40+10*(i-1)),'pA']; 
end
%legend(legendInfo);
for i=1:N_trace
    %This way, everything is with the units ms and mv.     
    T_RN.vf(i,:)=fix_ephys_jumps_simple(T_RN.t(i,:),T_RN.v(i,:),[15 215]);
end
base_vm=mean(T_RN.vf(5,1:200));
min_v=0;
max_v=-100;
for i=1:N_trace
    T_RN.vfa(i,:)=T_RN.vf(i,:)-T_RN.vf(i,1)+base_vm;
    text(30+(i-1)*15,mean(T_RN.vfa(i,1000:2000)),num2str(i))
    min_v=min(min_v,min(T_RN.vfa(i,:)));
    max_v=max(max_v,max(T_RN.vfa(i,:)));
    fit_v(i)=mean(T_RN.vfa(i,rnfit_start:rnfit_end));
end
rest_vm=mean(T_RN.vfa(5,1:200));
min_v=min_v-1;
max_v=max_v+1;
vvv=[190 min_v;215 min_v;215 max_v;190 max_v];
fff=[1 2 3 4];

min_RN=min([min(T_RN.v(1,:))-1 min(T_RN.vf(1,:))-1 min(T_RN.vfa(1,:))-1]);
max_RN=max([max(T_RN.v(N_trace,:))+1 max(T_RN.vf(N_trace,:))+1 max(T_RN.vfa(N_trace,:))+1]);
subplot(2,2,1)
for i=1:N_trace
    plot(T_RN.t(i,:),T_RN.v(i,:))
    hold on;
    %legendInfo{i} = ['Trace# ',num2str(i),' represents ',num2str(-40+10*(i-1)),'pA']; 
end
%legend(legendInfo);
title([Neuronlist{neuron_count},' Original']);
xlim([0 300])
ylim([min_RN max_RN])
xlabel('Time (ms)');
ylabel('Vm (mV)');
subplot(2,2,2)
for i=1:N_trace
    plot(T_RN.t(i,:),T_RN.vf(i,:))
    hold on;
end
xlim([0 300])
ylim([min_RN max_RN])
xlabel('Time (ms)');
ylabel('Vm (mV)');
title('After jumpfixed');
subplot(2,2,3)
min_v=0;
max_v=-100;
for i=1:N_trace
    plot(T_RN.t(i,:),T_RN.vfa(i,:))
    text(30+(i-1)*15,mean(T_RN.vfa(i,1000:2000)),num2str(i))
    hold on;    
end
patch('Faces',fff,'Vertices',vvv,'FaceColor','red','FaceAlpha',.2);
text(220,min_v+2,'Fitting region')
xlim([0 300])
ylim([min_RN max_RN])
xlabel('Time (ms)');
ylabel('Vm (mV)');
title('After jumpfixed and re-aligned');

subplot(2,2,4)
%first we use all nine traces for fitting%
fit_x=[ones(length(fit_i),1) fit_i'];
b=fit_x\fit_v;
RN_allfit=b(2)*1000;
scatter(fit_i,fit_v);
hold on;
for i=1:length(fit_i)
    text(fit_i(i)+2,fit_v(i),num2str(i))
end
vCalc=fit_x*b;
plot(fit_i,vCalc,'b--')
xlabel('Time (ms)');
ylabel('Averaged Vm (mV)');
%%%%% Method 2, select which traces we want to include for RN
%%%%% fitting.
%%%%% Trace 1 -40pA and Trace 9 +40pA
list = {'Trace 1','Trace 2','Trace 3',...                   
'Trace 4','Trace 5','Trace 6',...
'Trace 7','Trace 8','Trace 9'};
[indx,tf] = listdlg('PromptString',[Neuronlist{neuron_count},' Select traces for RN fitting'],...
                           'ListString',list);
%This way indx contains all the traces we want for RN fitting. 
fit_i_opt=fit_i(indx);
fit_v_opt=fit_v(indx);
%End of Method 2
%%%%%%%%%%%

cd(['../',Neuronlist{neuron_count},'_figs']);
saveas(gcf,[Neuronlist{neuron_count},'_Fixandfit.png']);
close;

figure(11)
fit_x=[ones(length(fit_i),1) fit_i'];
b=fit_x\fit_v;
RN_allfit=b(2)*1000;
scatter(fit_i,fit_v);
hold on;
for i=1:length(fit_i)
    text(fit_i(i)+2,fit_v(i),num2str(i))
end
vCalc=fit_x*b;
plot(fit_i,vCalc,'b--')
fit_x1=[ones(length(fit_i_opt),1) fit_i_opt'];
b1=fit_x1\fit_v_opt;
RN_fifit=b1(2)*1000;
v2Calc_1=fit_x1*b1;
plot(fit_i_opt,v2Calc_1,'r-')
legend('mean Vm','Fit for all','Final fit','Location','northwest');
text(-30,mean(fit_v),num2str(fit_i_opt))
text(-30,mean(fit_v)+1,'We selected these traces for RN fitting (red circles)')
plot(fit_i_opt,fit_v_opt,'ro','MarkerFaceColor', 'r')
ttl1 = sprintf('RN fit all traces: %g M-Ohm;Best fitting RN: %g M-Ohm\n',RN_allfit,RN_fifit);
ttl2 = sprintf('Resting membrane potential: %g mV',rest_vm);
title([ttl1 ttl2]);
xlabel('Time (ms)');
ylabel('Vm (mV)');
saveas(gcf,[Neuronlist{neuron_count},'_RN_fit.png']);
close;


figure(12)
for i=1:N_trace
    plot(T_RN.t(i,:),T_RN.v(i,:),'r')
    hold on;
    plot(T_RN.t(i,:),T_RN.vf(i,:),'b')
    plot(T_RN.t(i,:),T_RN.vfa(i,:),'black') 
end
legend('Original','After Fixed','Fixed and Re-aligned');
xlabel('Time (ms)');
ylabel('Vm (mV)');
saveas(gcf,[Neuronlist{neuron_count},'_RN_all.png']);
saveas(gcf,[Neuronlist{neuron_count},'_RN_all.fig']);
close;
cd ..
%%%% section 1 end

%%%% section 2 exponential fitting of tau
% This is to read in the fourth trace of RN for tau_fitting as base
% N_trace: # of traces for RN sweep
% N_tau: # of trace for tau_fitting!!! 
% Baseline: -100 to 0
% Slope: arbitrary
% Tau: 
N_tau=4;
% t=exp_tv(250:2750,1);
% v=exp_tv(250:2750,2);
tau_ti=eval(Name_temp(2*N_tau).name).*1000;
t_tau=tau_ti(:,1);
i_tau=tau_ti(:,2);
tau_cur=round(10^9*(mean(i_tau(250:2550))-mean(i_tau(1:20))));
t=T_RN.t(N_tau,250:2550);
v=T_RN.vf(N_tau,250:2550);
%use exp_tv for exponential fit!!!
cd([Neuronlist{neuron_count},'_figs']);
[fitresult, gof] = createFit(Neuronlist{neuron_count},t, v,tau_cur);
tau=fitresult.tau;

%%%% section 2 end
else
    figure()
    title('RN empty')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'RN empty.png');
    close;
end
else
    figure()
    title('No RN protocol found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No RN.png');
    close;
end
cd ..

% Section 3: analyze high_RN protocol. 
cd(Neuronlist{neuron_count});
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
Cur_temp=roundn(mean(Trace_cur(5000:10000,2))*10^12,1);
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
Cur(i_hi)=roundn(mean(Trace_cur(5000:10000,2))*10^12,1);
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

% Section 4: analyze low_RN protocol. 
cd(Neuronlist{neuron_count});
%save AP count for -100:20:120 pA
low_mahp=0;
low_sahp=0;
Low_apcount=zeros(12,2);
Low_apcount(:,1)=-220:50:330;
%%%%
if(exist(file_low,'file'))
load(file_low);
% the purpose for analyzing High/Low RN is to compute ISI, ins_freq, delay
% and AP num. In this case, we also ignore the doublets since during
% doublets, the system is not stably oscillating. 
mat_low=matfile(file_low);
Name_low=whos(mat_low);
if (isempty(Name_low)~=1)
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
Tracelow=eval([hiname_template,num2str(i_hi),'_1']);
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
Cur_temp=roundn(mean(Trace_cur(5000:10000,2))*10^12,1);
v_ap=Tracelow(:,2).*1000;
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
    [delay APnum ISI freq f_ahp m_ahp s_ahp]=hilo_ana(Neuronlist{neuron_count},2,Cur(i_hi),t_ap,v_ap,dv_ap,d2v_ap,dv_thre);
    low_mahp=low_mahp+mean(m_ahp);
    low_sahp=low_sahp+mean(s_ahp);
    Cur_ap(pos_ap)=Cur(i_hi);
    Freq_ISI(pos_ap)=mean(1000./ISI);
    Freq_count(pos_ap)=APnum/2;
    %save apcount for LowRN to Low_apcount for later file output!
    lowid_temp=find(Low_apcount(:,1)>=Cur(i_hi),1);
    Low_apcount(lowid_temp,2)=APnum;
    % save complete, if no AP then APnum=0
    fprintf(f_low,'%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t',Cur(i_hi),APnum,delay,mean(ISI),std(ISI),std(ISI)/mean(ISI),Freq_ISI(pos_ap),Freq_count(pos_ap),f_ahp,mean(m_ahp),mean(s_ahp));
    % attach the spiking times to each line!
    spkt=zeros(APnum,1);
    spkt(1)=delay;
    for apnumtemp=2:APnum
        spkt(apnumtemp)=spkt(apnumtemp-1)+ISI(apnumtemp-1);
    end
    for apnumtemp=1:APnum
        fprintf(f_low,'%f\t',spkt(apnumtemp));
    end
    fprintf(f_low,'\n');
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
fclose(f_low);
cd(['../',Neuronlist{neuron_count},'_figs']);
saveas(gcf,[Neuronlist{neuron_count},'_Low_RN.png']);
saveas(gcf,[Neuronlist{neuron_count},'_Low_RN.fig']);
close;

% Section 4.5 Compute sag
i_low=1;
TraceLow=eval([hiname_template,num2str(i_low),'_1']);
CurLow=eval([hiname_template,num2str(i_low),'_2']);
Cur_sag=roundn(mean(CurLow(5000:10000,2))*10^12,1);
v_ap=TraceLow(:,2).*1000;
t_ap=TraceLow(:,1).*1000;
[sag_v id_sag]=min(v_ap);
rnfit_start=21439;%floor(0.45*l_rn);
rnfit_end=23949;%floor(0.55*l_rn);
stable_v=mean(v_ap(rnfit_start:rnfit_end));
sag=sag_v-stable_v;
figure()
plot(t_ap,v_ap);
hold on;
plot(t_ap(id_sag),v_ap(id_sag),'ro');
text(t_ap(id_sag),v_ap(id_sag)+20,'This is the sag point');
title(['Sag is computed at ',num2str(Cur_sag),'pA'])
xlabel('Time (ms)');
ylabel('Vm (mV)');
saveas(gcf,[Neuronlist{neuron_count},'sag.png']);
close;

low_mahp=low_mahp/num_aplevel;
low_sahp=low_mahp/num_aplevel;
else
    figure()
    title('LowRn empty')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'LowRn empty.png');
    close;
    low_mahp=NaN;
    low_sahp=NaN;
end
else
    figure()
    title('No LowRn protocol found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No LowRn.png');
    close;
    low_mahp=NaN;
    low_sahp=NaN;
end
cd ..


%%%section 5:AP analysis single trace
cd(Neuronlist{neuron_count});
if (exist(file_AP,'file'))
load(file_AP);
mat_ap=matfile(file_AP);
Name_ap=whos(mat_ap);
if (isempty(Name_ap)~=1)
clear mat_ap
N_ap=length(Name_ap)/2;
Trace_temp=eval(Name_ap(2*1-1).name);
l_ap=length(Trace_temp(:,1));
T_ap=table([1:N_ap]',zeros(N_ap,l_ap),zeros(N_ap,l_ap),'VariableNames',{'Trace','t','v'});
Cur=zeros(N_ap,1);
str_temp=Name_ap(1).name;
pos_underscore=strfind(str_temp,'_');
hiname_template=str_temp(1:pos_underscore(3)); 
num_aplevel=0;
min_v=100;
i_hi=1;
%seperate the trace we want for AP analysis
marker_ap=0;
while (i_hi<N_ap)
pos_ap=0;
Cur_ap=zeros(num_aplevel,1);
Traceap=eval([hiname_template,num2str(i_hi),'_1']);
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
l_cur=length(Trace_cur);
Cur_temp=roundn(mean(Trace_cur(round(0.3*l_cur):round(0.5*l_cur),2))*10^12,1);

v_ap=Traceap(:,2).*1000;
t_ap=Traceap(:,1).*1000;
if((max(v_ap)>0)&&(Cur_temp>0))
    [peak_temp id_peaktemp]=findpeaks(v_ap);
    ap_peaktemp=find(peak_temp>0);
    num_ptemp=length(find(peak_temp>0));
    t_aptemp=t_ap(id_peaktemp(ap_peaktemp));
    if ((num_ptemp>=3)&&(max(diff(t_aptemp))>20))
    marker_ap=1;
    break
    end
end
i_hi=i_hi+1;
end    

if(marker_ap~=0)
Trace_cur=eval([hiname_template,num2str(i_hi),'_2']);
Cur=roundn(mean(Trace_cur(1500:2500,2))*10^12,1);
l=length(v_ap);
dv_ap=get_dVdt(t_ap,v_ap);
d2v_ap=get_dVdt(t_ap(2:end-1),dv_ap(2:end));
t_bound=t_ap(end);
dv_thre=3;

cd(['../',Neuronlist{neuron_count},'_figs']);
subplot(1,2,1)
plot(t_ap,v_ap);
xlabel('Time (ms)');
ylabel('Vm (mV)');
title(['Full trace at ',num2str(Cur),'pA']);

[thre pe ampl riset decayt halft]=AP_ana(Neuronlist{neuron_count},t_ap,v_ap,dv_ap,d2v_ap,dv_thre);
else
    figure()
    title('No eligible AP trace found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No AP trace.png');
    close;
end

else
    figure()
    title('AP empty')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'AP empty.png');
    close;
end
else
    figure()
    title('No AP protocol found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No AP.png');
    close;
end   
cd ..
%%%%section 5 end

%%%section 6 Rheobase
cd(Neuronlist{neuron_count});
if(exist(file_rheo,'file'))
load(file_rheo);
mat_rheo=matfile(file_rheo);
Name_rheo=whos(mat_rheo);
if(isempty(Name_rheo)~=1)
N_rheo=length(Name_rheo)/2;
str_temp=Name_rheo(1).name;
pos_underscore=strfind(str_temp,'_');
rheoname_template=str_temp(1:pos_underscore(3));
Trace_temp=eval([rheoname_template,'1_1']);
ramp_temp=eval([rheoname_template,'1_2']);
t_rheo=Trace_temp(:,1)*1000;
v_rheo=Trace_temp(:,2)*1000;
i_rheo=ramp_temp(:,2)*10^12;
if (max(v_rheo)<=0)
    Trace_temp=eval([rheoname_template,'2_1']);
    ramp_temp=eval([rheoname_template,'2_2']);
    t_rheo=Trace_temp(:,1)*1000;
    v_rheo=Trace_temp(:,2)*1000;
    i_rheo=ramp_temp(:,2)*10^12;
end
if (max(v_rheo)>0)
    spk_st=find(v_rheo>=0,1);
    v_temp=v_rheo(spk_st:end);
    spk_end=find(v_temp<=0,1)+spk_st+1;
    [v_1spk id_1spk]=max(v_rheo(spk_st:spk_end));
    %locate the index for the first spike in Vm
    id_rheo=spk_st+id_1spk;
    %compute the ramping slope
    %Linear ramp should continue until the first spike, take id_rheo-20
    %One time step would be 0.33ms.
    %The ramping start at 500ms, ~1500 time steps, use 1550 to make sure
    %the ramp is already being applied! This could change among other
    %protocols! 
    slope_i=(i_rheo(id_rheo-20)-i_rheo(1550))/(t_rheo(id_rheo-20)-t_rheo(1550));
    base_i=mean(i_rheo(find((t_rheo<=500)&(t_rheo)>=200)));
    ramp_startid=find(t_rheo>500,1);
    rheo=base_i+slope_i*(t_rheo(id_rheo)-t_rheo(ramp_startid));
    figure()
    subplot(2,1,1)
    plot(t_rheo,v_rheo,'b');
    hold on;
    plot(t_rheo(id_rheo),v_rheo(id_rheo),'ro')
    xlabel('Time (ms)')
    ylabel('Vm (mV)')
    subplot(2,1,2)
    plot(t_rheo,i_rheo,'b');
    hold on;
    plot(t_rheo(id_rheo),rheo,'ro')
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    title(['Rheobase measured at slope',num2str(slope_i),'pA']);
else
    figure()
    subplot(2,1,1)
    plot(t_rheo,v_rheo,'b');
    subplot(2,1,2)
    plot(t_rheo,i_rheo,'b');
    title('Rheobase protocol does not work');
    rheo=0;
end
    
cd(['../',Neuronlist{neuron_count},'_figs']);
saveas(gcf,[Neuronlist{neuron_count},'_Rheo.png']);
close;
else
    figure()
    title('Rheo empty')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'Rheo empty.png');
    close;
end
else
    figure()
    title('No Rheo protocol found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No Rheo.png');
    close;
end 
cd ..


fout = fopen('output.txt', 'a+');
str_date=date;

fprintf(fout,'%s\t%s\t%f\t%f\t%f\t%f\t%f\t%s\t',str_date,Neuronlist{neuron_count},tau,RN_allfit,RN_fifit,b1(1),sag,rheo,f_ahp);
fprintf(fout,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t\t',0.5*(hi_mahp+low_mahp),0.5*(hi_sahp+low_sahp),thre,ampl,riset,decayt,halft);
for n_out=1:length(High_apcount)
    fprintf(fout,'%d\t',High_apcount(n_out,2));
end
fprintf(fout,'\t');
for n_out=1:length(Low_apcount)
    fprintf(fout,'%d\t',Low_apcount(n_out,2));
end
fprintf(fout,'\n');
% Another handy functino to find local peak 
% [pks2,locs2] = findpeaks(dv_ap,'MinPeakDistance',6);
% figure(2);plot(t_ap(2:end-1),dv_ap(2:end),t_ap(locs2),pks2,'or')

end


