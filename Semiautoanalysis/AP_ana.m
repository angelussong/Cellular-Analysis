function [thre pe ampl riset decayt halft]=AP_ana(Neuron_name,t_ap,v_ap,dv_ap,d2v_ap,dv_thre)

[peak_temp id_peaktemp]=findpeaks(v_ap);
ap_peaktemp=find(peak_temp>0);
num_ptemp=length(find(peak_temp>0));
t_aptemp=t_ap(id_peaktemp(ap_peaktemp));
v_aptemp=v_ap(id_peaktemp(ap_peaktemp));

id_peakap=id_peaktemp(ap_peaktemp);
AP_start=zeros(num_ptemp,1);

for i=1:num_ptemp
    [max_dv id_maxdv]=max(dv_ap(id_peakap(i)-10:id_peakap(i)));
    AP_start(i)=id_peakap(i)-10+id_maxdv-5;
end
% plot(t_ap,v_ap);
% hold on;
% plot(t_ap(AP_start),v_ap(AP_start),'o')

l_tmp=length(AP_start);
ind_tmp=[];
for i=1:l_tmp-1
    if ((t_ap(AP_start(i+1))-t_ap(AP_start(i))<20)&&(t_ap(AP_start(i+1))-t_ap(AP_start(i))>0.5))
        ind_tmp=[i i+1 ind_tmp];
    end
    if (t_ap(AP_start(i+1))-t_ap(AP_start(i))<0.5)
        ind_tmp=[i ind_tmp];
    end
end
AP_start(ind_tmp)=[];

l_tmp=length(AP_start);
ind_tmp3=[];
for i=1:l_tmp
    if((t_ap(AP_start(i))>=2000)||(t_ap(AP_start(i))<=20))
        ind_tmp3=[i ind_tmp3];
    end
end
AP_start(ind_tmp3)=[];

n_AP=length(AP_start);
l_AP=length(1:n_AP);
v_threh=zeros(l_AP,1);
ap_peak=zeros(l_AP,1);
ap_rise=zeros(l_AP,1);
ap_decay=zeros(l_AP,1);
ap_amp=zeros(l_AP,1);
ap_half=zeros(l_AP,1);

f_ap = fopen(['../',Neuron_name,'_APoutput.txt'], 'a+');
fprintf(f_ap,'AP#\tThreshold\tPeak\tAmp\tRiseT\tDecayT\tHalfDur\n');    

for i=1:n_AP
    
    if(i<=n_AP-1)
        if((v_ap(AP_start(i)-2)-v_ap(AP_start(i)-3)>=4)&&(v_ap(AP_start(i)-2)>-40))
            v_thre=v_ap(AP_start(i)-3);
        else
            v_thre=v_ap(AP_start(i)-2);
        end
    vap_temp=v_ap(AP_start(i):AP_start(i+1));
    else
        v_thre=v_ap(AP_start(i)-1);
        vap_temp=v_ap(AP_start(i):end);
    end
    id_downthre=find(vap_temp<v_thre,1);
    v_downthre=v_ap(AP_start(i)+id_downthre);
    [peak max_ind]=max(v_ap(AP_start(i):AP_start(i)+id_downthre));
    amp=peak-v_thre;
    subplot(2,2,2)
    plot(v_ap(AP_start(i):AP_start(i)+id_downthre))
    hold on;    
    t_rise=t_ap(AP_start(i)+max_ind)-t_ap(AP_start(i));
    t_decay=t_ap(AP_start(i)+id_downthre)-t_ap(AP_start(i)+max_ind);
    id_halfamp=find(v_ap(AP_start(i):AP_start(i)+id_downthre)>=0.5*amp+v_thre,1);
    id_halfend=find(v_ap(AP_start(i):AP_start(i)+id_downthre)>=0.5*amp+v_thre,1,'last');
    t_half=t_ap(AP_start(i)+id_halfend+1)-t_ap(AP_start(i)+id_halfamp);
    v_threh(i)=v_thre;
    ap_peak(i)=peak;
    ap_amp(i)=amp;
    ap_rise(i)=t_rise;
    ap_decay(i)=t_decay;
    ap_half(i)=t_half;
    fprintf(f_ap,'%d\t%f\t%f\t%f\t%f\t%f\t%f\n',i,v_thre, peak, amp, t_rise, t_decay, t_half);    
end
    fprintf(f_ap,'Average\t%f\t%f\t%f\t%f\t%f\t%f\n',mean(v_threh),mean(ap_peak),mean(ap_amp),mean(ap_rise),mean(ap_decay),mean(ap_half)); 
    fprintf(f_ap,'STD\t%f\t%f\t%f\t%f\t%f\t%f\n',std(v_threh),std(ap_peak),std(ap_amp),std(ap_rise),std(ap_decay),std(ap_half)); 
    fclose(f_ap);
%just return the first index, which corresponds to the second AP (first one
%being a doublet or unstable). 
if(length(v_threh)>=2)
    thre=v_threh(2);
    pe=ap_peak(2);
    ampl=ap_amp(2);
    riset=ap_rise(2);
    decayt=ap_decay(2);
    halft=ap_half(2);
    i=2;
else
    thre=v_threh(1);
    pe=ap_peak(1);
    ampl=ap_amp(1);
    riset=ap_rise(1);
    decayt=ap_decay(1);
    halft=ap_half(1);
    i=1;
end
ttl0=sprintf('Measurement for 2nd AP ');
ttl1=sprintf('threshold%g\t peak%g\t',thre,pe);
ttl2=sprintf('amp%g\t rise%g\t',ampl,riset);
ttl3=sprintf('decay%g\t half%g',decayt,halft);
text(12,0,{ttl0;ttl1;ttl2;ttl3});
title('Superimposed APs');
xlabel('Time (ms)');
ylabel('Vm (mV)');

v_thre=v_ap(AP_start(i));
if (length(AP_start)>2)
vap_temp=v_ap(AP_start(i):AP_start(i+1));
end
if (length(AP_start)==2)
vap_temp=v_ap(AP_start(i):end);
end
id_downthre=find(vap_temp<v_thre,1);
v_downthre=v_ap(AP_start(i)+id_downthre);
[peak max_ind]=max(v_ap(AP_start(i):AP_start(i)+id_downthre));
amp=peak-v_thre;
t_rise=t_ap(AP_start(i)+max_ind)-t_ap(AP_start(i));
    t_decay=t_ap(AP_start(i)+id_downthre)-t_ap(AP_start(i)+max_ind);
    id_halfamp=find(v_ap(AP_start(i):AP_start(i)+id_downthre)>=0.5*amp+v_thre,1);
    id_halfend=find(v_ap(AP_start(i):AP_start(i)+id_downthre)>=0.5*amp+v_thre,1,'last');
    t_half=t_ap(AP_start(i)+id_halfend)-t_ap(AP_start(i)+id_halfamp);


subplot(2,2,4)    
plot(t_ap(AP_start(i):AP_start(i)+id_downthre),v_ap(AP_start(i):AP_start(i)+id_downthre))
hold on;
plot(t_ap(AP_start(i)),v_thre,'bo')
text(t_ap(AP_start(i)),v_thre,'\leftarrow Threshold')
plot(t_ap(AP_start(i)+max_ind),peak,'ro')
text(t_ap(AP_start(i)+max_ind)-0.1,peak,'\leftarrow Peak')
plot([t_ap(AP_start(i)+max_ind-1) t_ap(AP_start(i)+max_ind-1)],[v_thre peak]);
text(t_ap(AP_start(i)+max_ind),0.5*(v_thre+peak)+10,'\leftarrow Amplitude')
plot(t_ap(AP_start(i)+id_halfamp)-0.1,0.5*amp+v_thre,'go')
text(t_ap(AP_start(i))-0.6,-5,'Half rise point')
plot(t_ap(AP_start(i)+id_halfend),0.5*amp+v_thre,'go')
text(t_ap(AP_start(i)+id_halfend)+0.2,-5,'Half decay point')
plot([t_ap(AP_start(i)+id_halfamp) t_ap(AP_start(i)+id_halfend)],[0.5*amp+v_thre 0.5*amp+v_thre]);
text(t_ap(AP_start(i)+id_halfamp)+0.3,0.5*amp+v_thre-20,'Half activation time')
title('example AP (2nd)')
xlabel('Time (ms)');
ylabel('Vm (mV)');
saveas(gcf,[Neuron_name,'_AP_analysis.png'])

close



