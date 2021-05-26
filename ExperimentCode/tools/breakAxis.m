function [ax, ax2] = breakAxis(ax, gap)

sep = .02;

ax2 = copyobj(ax, ax.Parent);
ylim('auto');
yl = ylim;
xl = xlim;
gap = sort(gap);

lorange = [yl(1) gap(1)];
hirange = [gap(2) yl(2)];

loamt = diff(lorange);
hiamt = diff(hirange);

lofrac = loamt / (loamt + hiamt);
hifrac = hiamt / (loamt + hiamt);

x = ax.Position(1); y = ax.Position(2);
w = ax.Position(3); h = ax.Position(4);

ax.Position = [x y+h*lofrac+sep/2 w h*hifrac-sep/2];
ylim(ax, hirange);
xlim(ax, xl);
set(ax, 'XTick', [], 'XColor', 'none');

ax2.Position = [x y w h*lofrac-sep/2];
ylim(ax2, lorange);
xlim(ax2, xl);

end