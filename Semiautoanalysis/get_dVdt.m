function [dVdt] = get_dVdt(tvec, vvec)
%
%   get_dVdt    use the central difference formula to estimate the
%               derivative of the second vector with respect to the first. 
%
%   INPUT       TVEC, vector of times
%               VVEC, vector of membrane potential, for example.
%   OUTPUT      dVdt, the derivative of VVEC with respect to TVEC.
%
%   Christina Weaver, christina.weaver@mssm.edu, July 2005
%

% fprintf('Computing dV/dt: ');
for i = 2:length(tvec)-1 
    if( floor(i/50000)==i/50000) fprintf('%d...',i);  end;
    dVdt(i) = (vvec(i+1)-vvec(i-1))/(tvec(i+1)-tvec(i-1));
%     dVall(i,1) = (vvec(i+1)-vvec(i-1));
%     dVall(i,2) = (tvec(i+1)-tvec(i-1));
%     dVall(i,3) = (vvec(i+1)-vvec(i-1))/(tvec(i+1)-tvec(i-1));
%     if( mod(i,5000) == 0 ) fprintf('%d...',i); end;
end;

% fprintf('done\n');