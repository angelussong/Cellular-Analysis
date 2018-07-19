%   script test_RN:  for computing the Vr and RN from FitMaster data, which
%   has been saved as .mat files.
%
%   code written for the Luebke Lab by Dr Hanbing Song, June 2018.  If you
%   have questions, contact Hanbing at  hsong1@fandm.edu, or Dr Christina
%   Weaver at christina.weaver@fandm.edu .
%


clear
clc
%
%  INSTRUCTIONS FOR USER:  Type in the names of neurons to be analyzed in the cell array
%  like this:
%       Neuronlist={'Mar10IR1c','Mar10IR1d','Apr3IR2e'};
%
Neuronlist={'Mar10IR1c'};

for neuron_count=1:length(Neuronlist)
    % set default values to all the outputs to 0 so in the file if anything is
    % 0, that protocol is bad!!
    RN_allfit=NaN;
    RN_fifit=NaN;
    b1=[NaN NaN];
    % the code below only runs the RN protocol.  Don't look for the other
    % files.
%     file_AP = sprintf([Neuronlist{neuron_count},'-AP.mat']);
    file_RN = sprintf([Neuronlist{neuron_count},'-Rn.mat']);
%     file_High=sprintf([Neuronlist{neuron_count},'-HighRn.mat']);
%     file_low=sprintf([Neuronlist{neuron_count},'-LowRn.mat']);
%     file_rheo = sprintf([Neuronlist{neuron_count},'-ramp.mat']);
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
            xlabel('Injected Current (pA)');
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
            xlabel('Injected Current (pA)');
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
    
    % Hanbing:  if the output_RN.txt file doesn't exist yet, please create
    % a new one that includes the column headers.
    %
    %   also compute the Vr from the RN data -- the value of the regression
    %   when the injected current is zero. 
    %   Vr is stored at b1(1)
    fout = fopen('output_RN.txt', 'a+');
    str_date=date;
    fprintf(fout,'%s\t%s\t%f\t%f\t%f\n',str_date,Neuronlist{neuron_count},RN_allfit,RN_fifit,b1(1));
    
end
