function plot_llt(collocation_points,vertices,b,c,offset,particoes)
clf
hold on
xlabel('Metros')
ylabel('Metros')
xlim([0 1.5])
ylim([-0.75 0.75])
aux = 0;
offset = [0 offset];
aux_off = offset(1);

plot([0 0],[0 -c(1)],'k');
for i=1:particoes
    % BORDO DE ATAQUE
    plot([aux b(i)+aux],[-offset(i) -offset(i+1)],'k');
    plot([-aux -(b(i)+aux)],[-offset(i) offset(i+1)],'k');
    
    % CORDAS
    plot([aux+b(i) b(i)+aux],[-(offset(i+1)) -(c(i+1)+offset(i+1))],'k');
    plot([-(aux+b(i)) -b(i)-aux],[-(offset(i+1)) -(c(i+1)+offset(i+1))],'k');
    
    % BORDO DE FUGA
    plot([aux b(i)+aux],[-(c(i)+offset(i)),-(offset(i+1)+c(i+1))],'k');
    plot([-aux -(b(i)+aux)],[-(c(i)+offset(i)),-(offset(i+1)+c(i+1))],'k');
    
    aux = aux + b(i);
    aux_off = aux_off + offset(i);
end
%     plot(collocation_points(:,2),-collocation_points(:,1),'k');
scatter(collocation_points(:,2),-collocation_points(:,1),'r');
scatter(vertices(:,2),-vertices(:,1),'g');
yo=1;
end