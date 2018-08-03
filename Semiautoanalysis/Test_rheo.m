clear
clc
Neuronlist={'Jun19IR4c'};

neuron_count=1;
    
file_rheo = sprintf([Neuronlist{neuron_count},'-ramp.mat']); 
cd(Neuronlist{neuron_count});
f_rheo = fopen('../Rheobase.txt', 'a+');
  
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
num_rheo=1;
marker_rheo=0;
while (max(v_rheo)<=0)&&(exist([rheoname_template,num2str(num_rheo+1),'_1']))
    Trace_temp=eval([rheoname_template,num2str(num_rheo+1),'_1']);
    ramp_temp=eval([rheoname_template,num2str(num_rheo+1),'_2']);
    t_rheo=Trace_temp(:,1)*1000;
    v_rheo=Trace_temp(:,2)*1000;
    i_rheo=ramp_temp(:,2)*10^12;
    if (max(v_rheo)>0)
        marker_rheo=1;
    end
    num_rheo=num_rheo+1;
end
if ((max(v_rheo)>0)||(marker_rheo==1))
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
    fprintf(f_rheo,'%s\t%d\n',Neuronlist{neuron_count},rheo);  
else
    figure()
    subplot(2,1,1)
    plot(t_rheo,v_rheo,'b');
    subplot(2,1,2)
    plot(t_rheo,i_rheo,'b');
    title('Rheobase protocol does not work');
    fprintf(f_rheo,'%s\tNan\n',Neuronlist{neuron_count});  
end
mkdir(['../',Neuronlist{neuron_count},'_figs']);    
cd(['../',Neuronlist{neuron_count},'_figs']);
saveas(gcf,[Neuronlist{neuron_count},'_Rheo.png']);
close;
else
    figure()
    title('Rheo empty')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'Rheo empty.png');
    close;
    fprintf(f_rheo,'%s\tNan\n',Neuronlist{neuron_count});  
end
else
    figure()
    title('No Rheo protocol found')
    cd(['../',Neuronlist{neuron_count},'_figs']);
    saveas(gcf,'No Rheo.png');
    close;
    fprintf(f_rheo,'%s\tNan\n',Neuronlist{neuron_count});  
end 
cd ..


