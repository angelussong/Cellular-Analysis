function [vFix] = fix_ephys_jumps_simple(tvec, vvec, manBds)
%   fix_ephys_jumps_simple
%
%   BE SURE THE DATA SENT TO THIS FUNCTION IS REPORTED IN MS AND MV.  OFTEN
%   THIS MEANS YOU MUST SCALE THE DATA BY 1000 IN BOTH time AND voltage.
%
%   Jennie noted that, at the very start and end of a current step, the
%   series resistance causes an artifact of an unphysiological voltage
%   jump.  She corrects it by just subtracting out the start
%   and end voltage during a short window.  I have thought about more
%   elaborate ways to correct for it (Cengiz Gunay's 2011 CNS poster may
%   have done so too), but she says not to worry about it.
%
%   input   tvec        vector of times, MEASURED IN MS
%           vvec        vector of voltage data, MEASURED IN MV
%           manBds      boundaries of the current step, MEASURED IN MS,
%                       e.g. [15 215]
%
%   output  vFix        the adjusted voltage trace
%
%   modified from fix_ephys_jumps, to make this function easier to use.
%
%
%   Christina Weaver, christina.weaver@fandm.edu
%   June 2012

% autoBds = [];

dv = get_dVdt(tvec,vvec);
dt = tvec(2)-tvec(1);
% 'tst' stores index of the last time point before the time window specified in 'manBds'
% cw 6/7/18
tst = max(find(tvec <= manBds(1) ));

% t0:  index of the first time point that is 0.25 ms after time 'tst'
% cw 6/7/18
t0 = min(find(tvec-tvec(tst) > 0.25 ));
% compute mean voltage for 10 ms before the manBds time window - cw 6/7/18
tprev = find(tvec>tvec(tst)-10 & tvec < tvec(tst));
vmean = mean(vvec(tprev));
fprintf('Jump starts at %g ms and ends at %g ms\nPrevious Vmean = %g\n',...
    tvec(tst),tvec(t0),vmean);
voff = vmean - vvec(t0);
fprintf('V at end = %g\nVoffset = %g\n',vvec(t0),voff);

% 'tend' stores index of the last time point within the 'manBds' time window 
% cw 6/7/18
tend = max(find(tvec <= manBds(2) ));

t1 = min(find(tvec-tvec(tend) > 0.25 ));
tprev2 = find(tvec>tvec(tend)-10 & tvec < tvec(tend));
fprintf('End of pulse detected at %g ms; fix until %g ms\npulse length %g ms',...
    tvec(tend),tvec(t1),tvec(t1)-tvec(t0));

fprintf('V(%g) = %g\tV(%g) = %g\n',tvec(tend),vvec(tend),tvec(t1),vvec(t1));
voff = vmean - vvec(t0);
fprintf('V at end = %g\nVoffset = %g\n',vvec(t0),voff);

vFix = vvec;
vFix(tst:t0) = vmean;
vFix(t0+1:t1) = vFix(t0+1:t1) + voff;
vmean2 = mean(vFix(tprev2));

vFix(tend:t1) = vmean2;
dvFix = get_dVdt(tvec,vFix);

return;