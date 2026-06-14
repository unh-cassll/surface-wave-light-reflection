function [deltak]=Elfouhaily_spread(k,u10,omega,us)
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% Delft University of Technology published work.
% This function has been developped in the frame of Franco Fois's PhD.
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% Function Name:                  Elfouhaily_spread.m

% Description:                    This function computes the angular
%                                 spreading functionof the Elfouhaily sea 
%                                 surface spectrum.
%
% Original Author:                Franco Fois.

% Revision History:

%   Rev:      Author:       Date:         Description of change:
%   r0        FF            04/09/13      Preliminary Version.

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%______________________________________________________
% INPUT PARAMETERS DESCRIPTION
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% 1) k           : vector of wavenumbers [rad/m];
% 2) u10         : wind speed at 10 m height [m/s];
% 3) omega       : wave age [-];

%______________________________________________________
% OUTPUT PARAMETERS DESCRIPTION
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% 1) deltak     : angular spreading function;


% CD10=1e-3*(0.8+0.065*u10);
% us=sqrt(CD10)*u10;

cp=u10./omega;
g=9.81;
km=363;
c=sqrt(g./k.*(1+(k./km).^2));
cm=sqrt(2*g/km);
a0=log(2)/4;
ap=4;
am=0.13*us/cm;

deltak=tanh(a0+ap.*(c./cp).^2.5+am.*(cm./c).^2.5);

% Elfouhaily spreading function is: {PH 06Jan2015}
%  D(k,th)=[1+deltak cos(2*th)]; eq. (49), JGRv102p1578Elfouhaily
%  th=[-pi,pi]

return

