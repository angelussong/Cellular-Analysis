clear
clc
Ntype='WTD2';
Nname='msnWTD2_May3IR2a';

A=[30 80 130];

for i=1:3
ephysname=importdata(['ephys/',Nname,'/',Nname,'_2000ms_',num2str(A(i)),'.dat']);
t_ephys=ephysname(2:end,1);
v_ephys=ephysname(2:end,2);

fname = sprintf('%dpA.dat',A(i));
dstrg = 'native';
fin = fopen(fname,'r',dstrg);
[npts] = fread(fin,1,'double');
[t, csz] = fread(fin,npts,'double');
[idat, csz] = fread(fin,npts,'double');
subplot(2,2,i)
plot(t_ephys,v_ephys,'b','linewidth',1);
hold on;
plot(t,idat,'r','linewidth',2)
%xlim([0 2215])
hold on;
if (i==2)
    subplot(2,2,4)
    plot(t_ephys,v_ephys,'b','linewidth',1);
hold on;
plot(t,idat,'r','linewidth',2)
xlim([200 650])
hold on;
end
end
title(Nname);
xlabel('time (ms)');
ylabel('V (mV)');

FIGNAME=[Nname,'.fig'];
saveas(gcf,FIGNAME)

