%produces an animated surface plot of Temperatue or cure

for count=1:te/recordstep


figure(1)

surf(Tm(:,:,count))
%plot(0:deltaz:hz,Tm(:,1,count))
%plot(Tm(:,1,count))
title(num2str(count))
colorbar
%colorbar([0 250])
%axis([-inf inf -inf inf 25 50])
%figure(2)
%surf(alpham(:,:,count))
%pause(0.2)
%colormap(gray)

end

%figure(2)
%contourf(alpham(:,:,count),3000)
%hold on
%contourf(alpham(:,:,count),5000)