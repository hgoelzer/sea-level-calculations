% plot column evoltion in time
function plot_columns(BED,BASE,SURFACE,index)

%SURFACE
%THICK
%BASE
%BED
%GROUND_MASK

nb = length(BED);

set(groot,'defaultBarEdgeColor','w')
figure
hold on; box on;
for n = 1:nb
    bar(n, -3,'k')
    if BED(index,n)<0
        bar(n, BED(index,n),'FaceColor','b')
        bar(n, BASE(index,n),'FaceColor',[0.8,0.8,0.8])
        bar(n, SURFACE(index,n),'FaceColor',[0.8,0.8,0.8])
    else
        bar(n, BASE(index,n),'FaceColor',[0.8,0.8,0.8])
        bar(n, SURFACE(index,n),'FaceColor',[0.8,0.8,0.8])
        bar(n, BED(index,n),'FaceColor','k')
    end
end

plot([0,nb+1],[0,0],'--k','LineWidth', 1)
axis([0 nb+1 -3 3])

% save figure
%print -dpng -r300 p1.png


