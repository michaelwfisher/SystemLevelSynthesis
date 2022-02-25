
function fixPlot(number)

  figure(number);
  
  set(findall(gcf,'type','line'),'linewidth',2);
  set(gca,'Fontsize',16);
  ax = get(gca);
  set(ax.XLabel,'Fontsize',16);
  set(ax.YLabel,'Fontsize',16);
  set(ax.ZLabel,'Fontsize',16);

end