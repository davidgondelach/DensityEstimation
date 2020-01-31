function [mee] = fullmee2mee(Xp,nop,svs)
%FULLMEE2MEE - Return only objects states without BCs or ROM state

mee = zeros(nop*6,size(Xp,2));
for k = 1:nop
    mee(6*(k-1)+1:6*k,:) = Xp(svs*(k-1)+1:svs*(k-1)+6,:);
end

end

