% function PrecisionVSDistance
STDy = zeros(1,10);
dis = zeros(1,10);
hw = 31;
for i = 1:10
    pd = fitdist([Ceny{i,:}]','Normal');
    STDy(i) = 1/pd.sigma^2;
    Coor = [-16 15.5;-16 15.5];
    Coor(:,2) = Coor(:,2)-3*i+2;
    dis(i) = Distance_xy(Coor(1),Coor(2),Coor(3),Coor(4),2*hw+1);
end
%%
plot(dis,STDy,'o-','LineWidth',1.5)
xlabel('Distance ','fontSize',10)
ylabel('Precision','fontSize',10)
% end