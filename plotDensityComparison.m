function [orbitAvgRMSerr,dailyAvgRMSerr] = plotDensityComparison(time,rho_real,rho_rom,rho_jb2,rho_msise,rho_rom_std,latitudes,satName,varargin)

rho_rom_absErr = abs(rho_rom-rho_real)./rho_real*100;
rho_jb2_absErr = abs(rho_jb2-rho_real)./rho_real*100;
rho_msise_absErr = abs(rho_msise-rho_real)./rho_real*100;

% plotTime = (time-time(1))/86400;
plotTime = time;

densityPlot = figure;
plot(plotTime,rho_real); hold on;
plot(plotTime,rho_rom); hold on;
plot(plotTime,rho_jb2); hold on;
plot(plotTime,rho_msise); hold on;
xlabel('Day of year'); ylabel('\rho [kg/m^3]'); xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
legend(satName,'ROM','JB2008','NRLMSISE-00');

if ~isempty(rho_rom_std)
    covPlot = figure;
    plot(plotTime,rho_real); hold on;
    plot(plotTime,rho_rom); hold on;
    plot(plotTime,rho_rom - 3*rho_rom_std,'k--'); hold on;
    plot(plotTime,rho_rom + 3*rho_rom_std,'k--'); hold on;
    xlabel('Day of year'); ylabel('\rho [kg/m^3]'); xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
    legend(satName,'ROM','ROM \pm 3\sigma');
end

errPlot = figure;
plot(plotTime,rho_rom_absErr); hold on;
plot(plotTime,rho_jb2_absErr); hold on;
plot(plotTime,rho_msise_absErr); hold on;
xlabel('Day of year'); ylabel('Error \rho [%]');  xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
legend('ROM','JB2008','NRLMSISE-00');

% avgErrPlot = plotOrbitAverage(latitudes,plotTime,rho_rom_absErr,rho_jb2_absErr);
% xlabel('Day of year'); ylabel('Orbit-averaged error \rho [%]'); xlim([0 240]); xticks([0:24:240]);
% legend('ROM','JB2008');

% errPlot = figure;
% plot(movmean(plotTime,100),movmean(rho_rom_absErr,100)); hold on;
% plot(movmean(plotTime,100),movmean(rho_jb2_absErr,100)); hold on;
% xlabel('Day of year'); ylabel('Moving average error \rho [kg/m^3]'); xlim([0 240]); xticks([0:24:240]);
% legend('ROM','JB2008');

[avgErrPlot,orbitAvgRMSerr] = plotOrbitAveragedDensityError(latitudes,plotTime,rho_real,rho_rom,rho_jb2,rho_msise);
xlabel('Day of year'); ylabel('Orbit-averaged \rho error [%]');  xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
legend('ROM','JB2008','NRLMSISE-00');

avgDensPlot = plotOrbitAverage(latitudes,plotTime,rho_real,rho_rom,rho_jb2,rho_msise);
xlabel('Day of year'); ylabel('Orbit-averaged \rho [kg/m^2]');  xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
legend(satName,'ROM','JB2008','NRLMSISE-00');

[dailyAvgErrPlot,dailyAvgRMSerr] = plotDailyAveragedDensityError(plotTime,rho_real,rho_rom,rho_jb2,rho_msise);
xlabel('Day of year'); ylabel('Daily-averaged density error [%]');
legend('ROM','JB2008','NRLMSISE-00');  xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);

if nargin > 7
    SAVE_PLOTS = varargin{1};
    if SAVE_PLOTS
        filenameBase = varargin{2};
        savefig(densityPlot,[filenameBase satName '_Density.fig']);
        savefig(errPlot,[filenameBase satName '_DensityErr.fig']);
        savefig(avgErrPlot,[filenameBase satName '_DensityAvgErr.fig']);
        savefig(dailyAvgErrPlot,[filenameBase satName '_DensityDailyAvgErr.fig']);
        savefig(avgDensPlot,[filenameBase satName '_AvgDensity.fig']);
        if ~isempty(rho_rom_std)
            savefig(covPlot,[filenameBase satName '_DensityCov.fig']);
        end
    end
end

end

