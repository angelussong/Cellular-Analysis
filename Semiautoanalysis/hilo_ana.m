function [delay APnum ISI freq f_ahp m_ahp s_ahp]=hilo_ana(Neuronname,Hi_lo,Cur,t_ap,v_ap,dv_ap,d2v_ap,dv_thre)
f_ahp='no';
if (Hi_lo==2)
    hilo_header='LowRN';
else
    hilo_header='HighRN';
end

id_dvl = dv_ap>=dv_thre;
id_dvl(1) = 0;
id_dv = find(id_dvl);
yest = dv_ap(id_dv-1)<dv_thre; 
AP_starttemp=id_dv(yest);
lAP_temp=length(AP_starttemp);
AP_start=zeros(lAP_temp,1);
% for j=1:lAP_temp
%     [maxd2v id_d2v]=max(d2v_ap(AP_starttemp(j):min((AP_starttemp(j)+15),length(d2v_ap))));
%     if(maxd2v>=15)
%     AP_start(j)=AP_starttemp(j)+id_d2v-1;
%     end
% end
% AP_start=nonzeros(AP_start);
AP_start=AP_starttemp;
% t_ap(AP_start)

f=figure(15);
plot(t_ap,v_ap);
% t_ap(AP_start)
% xlim([0 500])
% figure(2);plot(t_ap(1:end-1),dv_ap)
% xlim([0 500])
% stop;
hold on;
% 
% [pks2,locs2] = findpeaks(dv_ap,'MinPeakDistance',6);
% 
% AP_start=locs2(find(pks2>50));

%%%%%%%%Criteria for doublet
%%%%%%%%between two threshold points, the time does not exceed cri_d ms. 
[v_pk id_pk]=findpeaks(v_ap);
v_pos=find(v_pk>0);
sus_dif=diff(t_ap(id_pk(v_pos)));
su_dif=sort(sus_dif);
id_cri=find(su_dif>=prctile(sus_dif,10),1);
cri_d=min(20,su_dif(id_cri));
l_tmp=length(AP_start);
ind_tmp=[];
for i=2:l_tmp
    % this part is to identify the doublets, not the precise location, but
    % enough to eliminate that from AP_ counts. Put this back in if
    % doublets are no longer interesting... 
%     if ((round(t_ap(AP_start(i))-t_ap(AP_start(i-1)))<cri_d)&&(t_ap(AP_start(i))-t_ap(AP_start(i-1))>0.5))
%         ind_tmp=[i i-1 ind_tmp];
%     end
    if (t_ap(AP_start(i))-t_ap(AP_start(i-1))<0.5)
        ind_tmp=[i ind_tmp];
    end
end

AP_start(ind_tmp)=[];
% figure(3)
% plot(t_ap,v_ap);
% hold on;
% plot(t_ap(AP_start(:)),v_ap(AP_start(:)),'o')

% l_tmp=length(AP_start);
% ind_tmp=[];
% for i=1:l_tmp-1
% temp_peak=findpeaks(v_ap(AP_start(i):AP_start(i+1)));
% if (length(find(temp_peak>0))>1)
%     ind_tmp=[i ind_tmp];
% end
% end
% % t_ap(AP_start)
% % ind_tmp
% % stop;
% AP_start(ind_tmp)=[];
% figure(2)
% plot(t_ap,v_ap);
% hold on;
% plot(t_ap(AP_start(:)),v_ap(AP_start(:)),'o')

l_tmp=length(AP_start);
ind_tmp=[];
for i=1:l_tmp-1
if (max(v_ap(AP_start(i):AP_start(i+1)))<-5)
    ind_tmp=[i ind_tmp];
end
end
AP_start(ind_tmp)=[];
% figure(2)
% plot(t_ap,v_ap);
% hold on;
% plot(t_ap(AP_start(:)),v_ap(AP_start(:)),'o')
% stop;


l_tmp=length(AP_start);
ind_tmp=[];
max_temp=max(v_ap(AP_start(l_tmp):end));
if (max_temp<=-5)
    ind_tmp=[l_tmp ind_tmp];
end

% t_ap(AP_start)
% ind_tmp
% stop;
AP_start(ind_tmp)=[];
% figure(1)
% plot(t_ap,v_ap);
% hold on;
% plot(t_ap(AP_start(:)),v_ap(AP_start(:)),'o')
% stop


n_AP=length(AP_start);
vv=0;
if(n_AP>1)
    APnum=n_AP;
    delay=t_ap(AP_start(1));
    freq=zeros(n_AP,1);
    ISI=zeros(n_AP-1,1);
    v_ahp=zeros(n_AP,1);
    v_thre=zeros(n_AP,1);
    m_ahp=zeros(n_AP,1);
    
    for i=1:n_AP
        
        if(i<=n_AP-1)
            v_thre(i)=v_ap(AP_start(i)-1);
            vap_temp=v_ap(AP_start(i):AP_start(i+1));
        else
            v_thre(i)=v_ap(AP_start(i)-1);
            vap_temp=v_ap(AP_start(i):end);
        end
    end
    for i=1:n_AP
        if(i<=n_AP-1)
            vap_temp=v_ap(AP_start(i):AP_start(i+1));
        else
            vap_temp=v_ap(AP_start(i):end);
        end
        id_downthre=find(vap_temp<v_thre(i),1);
        if(length(id_downthre)==0)
            id_downthre=find(vap_temp<v_thre(i+1),1);
        end
        v_downthre=v_ap(AP_start(i)+id_downthre);
        [peak max_ind]=max(v_ap(AP_start(i):AP_start(i)+id_downthre));
        j=0;
        while(v_ap(AP_start(i)+id_downthre+j+1)-v_ap(AP_start(i)+id_downthre+j)<0)
            j=j+1;
        end
        plot(t_ap(AP_start(i)+id_downthre+j+1),v_ap(AP_start(i)+id_downthre+j+1),'ro');vv=vv+1;
        f_ahp='yes';
        hold on;
        if(i<=n_AP-1)
            v_ahp(i)=min(v_ap(AP_start(i):AP_start(i+1)));
            m_ahp(i)=v_ahp(i)-v_thre(i);                   
        else
            if (min(v_ap(AP_start(i):25190))~=[])
            v_ahp(i)=min(v_ap(AP_start(i):25190));
            m_ahp(i)=v_ahp(i)-v_thre(i);
            else
                v_ahp(i)=0;
                m_ahp(i)=0;
            end
        end
        if (i<=n_AP-1)
            ISI(i)=t_ap(AP_start(i+1))-t_ap(AP_start(i)+id_downthre);
        end 
    end
        freq=1000/ISI(1);
        [sahp_v sahp_id]=min(v_ap);
        plot(t_ap(sahp_id),sahp_v,'bo');
        xlabel('Time (ms)');
        ylabel('Vm (mV)');
        s_ahp=min(sahp_v)-v_ap(1);
        m_ahp=nonzeros(m_ahp);
elseif(n_AP==1)
    APnum=1;
    delay=0;
    ISI=0;
    freq=0;
    m_ahp=0;
    s_ahp=0;
else
    APnum=0;
    delay=0;
    ISI=0;
    freq=0;
    m_ahp=0;
    s_ahp=0;  
end
%xlim([0 500])
title('Vm Trace: Red: fAHP, Blue: sAHP');
saveas(gcf,['../',Neuronname,'_figs/',Neuronname,'_',hilo_header,'_',num2str(Cur),'pA_fAHP.fig']);
close;
cd(['../',Neuronname,'/'])









