clear
clc


    tmp=importdata('Rheo.dat');
    t=tmp(:,1);
    v=tmp(:,2);
    l=length(t);
    vprime=zeros(l-1,1);
    for i=1:l-1
        vprime(i)=(v(i+1)-v(i))/(t(i+1)-t(i));
    end
%     plot(vprime)
    for i=2:l-1
        if (vprime(i)>=0.2)

            break;
        end
    end
    rheo=((t(i)-215)*0.0005-0.0005)*1000;
% WTD2
input_v=[-90.719559 
-88.246588 
-85.139182 
-80.997041 
-74.719869    ];
%computed RN = 135.9823, Rheo=109.15
%empirical RN=134.9, Rheo=83.2
% mean WTD2 RN=172.4, Rheo=109.8

%WTD1
input_v=[-90.073608 
-88.608264 
-86.96457 
-85.074987 
-82.889562];
%computed RN=89.8006, Rheo=156.3
%mean WTD1 RN=144.284, Rheo=156.944

%HETD1
input_v=[-93.278306 
-91.923424 
-90.341276 
-88.446223 
-86.094111];
% computed RN=89.8024, Rheo=129.55
% mean HETD1 RN = 186.8876, Rheo=127.4494

%HETD2
input_v=[-90.822567 
-88.484266 
-85.68567 
-82.250093 
-77.919052 ];
% computed RN=161.2939, Rheo=89.5
%mean RN=226.5159, Rheo=109.4349

% summary: RN from WTD2 to HETD2, 
% mean empirical increase 31%, model increase 20%

% mean rheo: WTD2: really good to mean: 109.15 to 109.8
% HETD2: computed 89.5 to emp mean: 109.4

% RN from WTD1 to HETD1, 
% mean empirical increase 31%, model almost the same

% mean rheo: WTD1: really good to mean: 156.3 to 156.944
% HETD1: also good, 129.55 to emp 127.45


RN=-mean(diff(input_v)./(-0.02));
RN
rheo
