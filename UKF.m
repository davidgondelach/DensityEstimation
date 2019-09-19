function [X_est,Pv] = UKF(X_est,Meas,time,stateFnc,measurementFcn,P,RM,Q)
%UKF Unscented Kalman filter
%   X_est: initial state guess
%   Meas: measurements
%   time: measurement times
%   stateFnc: function to propagate state
%   measurementFnc: function to convert from state to measurement space
%   P: state covariance matrix
%   RM: observation noise
%   Q: process noise

% Unscented Filter Parameter
% Compute the Sigma Points
[Wm,Wc,L,lam] = Unscented_Transform(X_est);
SR_Wc = sqrt(Wc); SR_Wm = sqrt(Wm);

S=chol(P)';
SR_R = sqrt(RM); % measurement noise
SR_Q = sqrt(Q); % process noise
eta = sqrt(L+lam);

try

m = size(Meas,2);
for i = 1:m-1
    
    fprintf('%.0f of %.0f \n',i,m-1)
    sigv=real([eta*S -eta*S]);
    xx=[X_est(:,i) sigv+kron(X_est(:,i),ones(1,2*L))];
    
    % Time Update
    [Xp] = stateFnc(xx,time(i),time(i+1));
    
    X_est(:,i+1) = Wm(1) * Xp(:,1) + Wm(2) * sum(Xp(:,2:end),2);
    
    % Get Propagated Square Root
    [~,S_minus]=qr([(SR_Wc(2)*(Xp(:,2:end)-kron(X_est(:,i+1),ones(1,2*L)))) SR_Q]',0);
    %     if S_minus(1,1) < 0; S_minus = -1.*S_minus; end
    S_minus=cholupdate(S_minus,Wc(1)*(Xp(:,1)-X_est(:,i+1)))';
    
    % Measurement function
    [Ym] = measurementFcn(Xp);
    ym = Wm(1) * Ym(:,1) + Wm(2) * sum(Ym(:,2:end),2);
    
    DY = Ym(:,1)-ym;
    DY(6:6:end) = wrapToPi(DY(6:6:end));
    
    DY2 = Ym(:,2:end)-kron(ym,ones(1,2*L));
    DY2(6:6:end) = wrapToPi(DY2(6:6:end));
    
    
    % Measurement Update
    [~,S_y]=qr([SR_Wc(2)*(DY2) SR_R]',0);
    %     if S_y(1,1) < 0; S_y = -1.*S_y;end
%     S_y=cholupdate(S_y,Wc(1)*(Ym(:,1)-ym))';
    S_y=cholupdate(S_y,Wc(1)*DY)';
    
    % Calculate Pxy
    %     Pyy0=Wc(1)*(Ym(:,1)-ym)*(Ym(:,1)-ym)';
%     Pxy0=Wc(1)*(Xp(:,1)-X_est(:,i+1))*(Ym(:,1)-ym)';
    Pxy0=Wc(1)*(Xp(:,1)-X_est(:,i+1))*DY';
    Pyymat=DY2;
    Pmat=Xp(:,2:end)-kron(X_est(:,i+1),ones(1,2*L));
    Pxy=Pxy0+Wc(2)*(Pmat*Pyymat');
    
    yres = Meas(:,i+1)-ym;
    yres(6:6:end) = wrapToPi(yres(6:6:end));
    
%     figure(123);
% subplot(2,3,1); plot(time(i+1),Meas(1:6:end,i+1),'b.'); hold on; plot(time(i+1),ym(1:6:end),'r.'); xlabel('t'); ylabel('p [km]');
% subplot(2,3,2); plot(time(i+1),Meas(2:6:end,i+1),'b.'); hold on; plot(time(i+1),ym(2:6:end),'r.'); xlabel('t'); ylabel('f');
% subplot(2,3,3); plot(time(i+1),(Meas(3:6:end,i+1)),'b.'); hold on; plot(time(i+1),ym(3:6:end),'r.'); xlabel('t'); ylabel('g');
% subplot(2,3,4); plot(time(i+1),(Meas(4:6:end,i+1)),'b.'); hold on; plot(time(i+1),ym(4:6:end),'r.'); xlabel('t'); ylabel('h');
% subplot(2,3,5); plot(time(i+1),(Meas(5:6:end,i+1)),'b.'); hold on; plot(time(i+1),ym(5:6:end),'r.'); xlabel('t'); ylabel('k');
% subplot(2,3,6); plot(time(i+1),rad2deg(Meas(6:6:end,i+1)),'b.'); hold on; plot(time(i+1),rad2deg(ym(6:6:end)),'r.'); xlabel('t'); ylabel('L [deg]');
% subplot(2,3,6); hold on; plot(time(i+1),X_est(7:7:end-10,i+1),'r.'); xlabel('t'); ylabel('BC');
    
    % Gain and Update
    KG = real(Pxy/S_y')/S_y;
%     X_est(:,i+1) = X_est(:,i+1) + KG * (Meas(:,i+1)-ym);
    X_est(:,i+1) = X_est(:,i+1) + KG * yres;
    U = KG * S_y;
    S = S_minus;
    for j = 1:length(ym)
        S = cholupdate(S',U(:,j),'-')';
    end
    
    HH(:,:,i+1)=(pinv(S*S')*KG*RM)';
    Pv(:,i+1) = diag(S*S')';
    
end

catch errMsg
    nowTimeStr = datestr(now,'yymmddHHMMSS');
    save(['workspace_UKF_' nowTimeStr]);
    rethrow(errMsg);
end

end

