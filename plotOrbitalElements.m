function [coe] = plotOrbitalElements(time,posvel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.4415;

nofStates = length(time);
coe = zeros(nofStates,6);
for i = 1:nofStates
    coe(i,:) = pv2po(posvel(1:3,i),posvel(4:6,i),mu);
end

figure;
subplot(2,3,1); plot(time,coe(:,1)); hold on; xlabel('t'); ylabel('a [km]');
subplot(2,3,2); plot(time,coe(:,2)); hold on; xlabel('t'); ylabel('e [-]');
subplot(2,3,3); plot(time,rad2deg(coe(:,3))); hold on; xlabel('t'); ylabel('i [deg]');
subplot(2,3,4); plot(time,rad2deg(coe(:,4))); hold on; xlabel('t'); ylabel('\Omega [deg]');
subplot(2,3,5); plot(time,rad2deg(coe(:,5))); hold on; xlabel('t'); ylabel('\omega [deg]');
subplot(2,3,6); plot(time,rad2deg(coe(:,6))); hold on; xlabel('t'); ylabel('\nu [deg]');
subplot(2,3,6); plot(time,rad2deg(coe(:,6))+rad2deg(coe(:,5))); hold on; xlabel('t'); ylabel('\nu [deg]');

end

