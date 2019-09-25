function [mee] = fullmee2mee(Xp,nop,svs)

mee = zeros(nop*6,size(Xp,2));
for k = 1:nop
    mee(6*(k-1)+1:6*k,:) = Xp(svs*(k-1)+1:svs*(k-1)+6,:);
end

end

