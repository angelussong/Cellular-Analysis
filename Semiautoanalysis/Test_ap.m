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


fout = fopen('output_ap.txt', 'a+');
str_date=date;

fprintf(fout,'%s\t%s\t%f\t%f\t%f\t%f\t%f\n',str_date,Neuronlist{neuron_count},thre,ampl,riset,decayt,halft);

end


