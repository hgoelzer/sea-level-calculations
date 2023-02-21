% plot column evoltion in time

%SURFACE
%THICK
%BASE
%BED
%GROUND_MASK

nb = nt+1;

set(groot,'defaultBarEdgeColor','w')
figure
hold on; box on;
for n = 1:nb
    bar(n, -3,'k')
    if BED(1,n)<0
        bar(n, BED(1,n),'FaceColor','b')
        bar(n, BASE(1,n),'FaceColor',[0.8,0.8,0.8])
        bar(n, SURFACE(1,n),'FaceColor',[0.8,0.8,0.8])
    else
        bar(n, BASE(1,n),'FaceColor',[0.8,0.8,0.8])
        bar(n, SURFACE(1,n),'FaceColor',[0.8,0.8,0.8])
        bar(n, BED(1,n),'FaceColor','k')
    end
end

plot([0,nb+1],[0,0],'--k','LineWidth', 1)
axis([0 nb+1 -3 3])

% save figure
%print -dpng -r300 p1.png


