function [thre pe ampl riset decayt halft]=AP_ana(Neuron_name,t_ap,v_ap,dv_ap,d2v_ap,dv_thre)
id_dvl = dv_ap>=dv_thre;
id_dvl(1) = 0;
id_dv = find(id_dvl);
yest = dv_ap(id_dv-1)<dv_thre; 
AP_starttemp=id_dv(yest);
lAP_temp=length(AP_starttemp);
AP_start=zeros(lAP_temp,1);
for j=1:lAP_temp
    [maxd2v id_d2v]=max(d2v_ap(AP_starttemp(j):min((AP_starttemp(j)+15),length(d2v_ap))));
    if(maxd2v>=75)
    AP_start(j)=AP_starttemp(j)+id_d2v-1;
    end
end
AP_start=nonzeros(AP_start);


%%%%%%%%Criteria for doublet
%%%%%%%%between two threshold points, the time does not exceed 20ms. 
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
ind_tmp2=[];
for i=1:l_tmp-1
temp_peak=findpeaks(v_ap(AP_start(i):AP_start(i+1)));
if (length(find(temp_peak>0))>1)
    ind_tmp2=[i ind_tmp2];
end
end
% t_ap(AP_start)
% ind_tmp
% stop;
AP_start(ind_tmp2)=[];

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

thre=v_threh(2);
pe=ap_peak(2);
ampl=ap_amp(2);
riset=ap_rise(2);
decayt=ap_decay(2);
halft=ap_half(2);
ttl0=sprintf('Measurement for 2nd AP ');
ttl1=sprintf('threshold%g\t peak%g\t',thre,pe);
ttl2=sprintf('amp%g\t rise%g\t',ampl,riset);
ttl3=sprintf('decay%g\t half%g',decayt,halft);
text(12,0,{ttl0;ttl1;ttl2;ttl3});
title('Superimposed APs');
xlabel('Time (ms)');
ylabel('Vm (mV)');

i=2;
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



