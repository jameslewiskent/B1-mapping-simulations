function plot_tx(Maps,settings)
    imagesc(imtile(abs(Maps),'Gridsize',[1 settings.Modes])) % Plot B1Tx FA maps for each mode
    set(gca,'Ydir','normal')
    axis image
    set(gca,'YTick',[]); set(gca,'XTick',[]);
    xlabel('Transmit Modes')
    xticks(((1:settings.Modes)).*size(Maps,2) - size(Maps,2)/2);
    xticklabels(1:settings.Modes);
end

