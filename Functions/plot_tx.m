function plot_tx(results,settings)
    imagesc(imtile(abs(settings.Tx_FA_map),'Gridsize',[1 settings.Modes])) % Plot B1Tx FA maps for each mode
    set(gca,'Ydir','normal')
    axis image
    set(gca,'YTick',[]); set(gca,'XTick',[]);
    xlabel('Transmit Modes')
    xticks(((1:settings.Modes)).*size(settings.Tx_FA_map,2) - size(settings.Tx_FA_map,2)/2);
    xticklabels(1:settings.Modes);
end

