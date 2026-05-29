
Re_lam = 2000;  Re_tur = 3000;
Re = linspace(100,1e6,1e4);
%r = 0;
r = 0;
D_inner = 0.009;

f_1 = 8*((8./Re).^12 + 1./((2.457*log((7./Re).^0.9 + 0.27*(r/D_inner))).^16 + (37530./Re).^16).^(3/2)).^(1/12);

f_2 = zeros(size(Re));
for i = 1:size(Re,2)
if Re(i) <= Re_lam
    f_2(i) = 64 / Re(i);
elseif Re(i) > Re_tur
    f_2(i) = 0.25*(log10(150.39/Re(i)^0.98865-152.66/Re(i)))^(-2);
else
    f_2(i) = (1.1525*Re(i) + 895)*1e-5;
end

end
figure(10)
hold on 
plot(Re,f_1)
plot(Re,f_2)
legend(["1";"2"])

hAxes = findobj(gcf,"Type","axes");         % 先获取图的对象
fontsize1 = 18;                             % 参数化·

for i = 1:1                                 % 进入循环，将所有子图的共性内容给处理掉，非共性的可以在外边处理
    hAxes(i).XScale = "log";
    %hAxes(i).YScale = "log";
    hAxes(i).FontName = "Times New Roman";
    hAxes(i).FontSize = fontsize1;
    hAxes(i).Box = "on";
    hAxes(i).BoxStyle = "full";
    hAxes(i).TickLength = [0.01 0.025];
end
