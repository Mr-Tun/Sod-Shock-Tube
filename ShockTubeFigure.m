%% Visualization of results in different time 
figure(1)                          %Pressure Distribution on the Space at Different Time
subplot(2,2,1)
epoch = 1:100:timestep;                   %select the epoch time
len_epoch = length(epoch);
legend_str = cell(1, len_epoch);
for  idx = 1:length(epoch)
    legend_str{idx} =[num2str(time(epoch(idx))*1000),'ms'];
    plot(space,phis(:,epoch(idx))/1000000,'LineWidth', 1)
    xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
    ylabel('\fontname{Times New Roman}\itp/\rmMPa','FontSize', 16);
    hold on
end
legend(legend_str,'location','best');
title('\fontname{Times New Roman}Pressure Distribution in the Space at Different Time');

subplot(2,2,2)                           % Density Distribution on the Space at Different Time
for  idx = 1:length(epoch)
    legend_str{idx} =[num2str(time(epoch(idx))*1000),'ms'];
    plot(space,rhohis(:,epoch(idx)),'LineWidth', 1)
    xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
    ylabel('\fontname{Times New Roman}\itρ/\rmkg/m³','FontSize', 16);
    hold on
end
legend(legend_str,'location','best');
title('\fontname{Times New Roman} Density Distribution in the Space at Different Time');

% subplot(2,2,3)                           % Temperature Distribution on the Space at Different Time
% for  idx = 1:length(epoch)
%     legend_str{idx} =[num2str(time(epoch(idx))*1000),'ms'];
%     plot(space,This(:,epoch(idx)),'LineWidth', 1)
%     xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
%     ylabel('\fontname{Times New Roman}\itT/\rmK','FontSize', 16);
%     hold on
% end
% legend(legend_str,'location','best');
% title('\fontname{Times New Roman} Temperature Distribution on the Space at Different Time');

subplot(2,2,4)                           % Velocity Distribution on the Space at Different Time
for  idx = 1:length(epoch)
    legend_str{idx} =[num2str(time(epoch(idx))*1000),'ms'];
    plot(space,uhis(:,epoch(idx)),'LineWidth', 1)
    xlabel('\fontname{Times New Roman}\itx/\rmm','FontSize', 16);
    ylabel('\fontname{Times New Roman}\itu/\rmm/s','FontSize', 16);
    hold on
end
legend(legend_str,'location','best');
title('\fontname{Times New Roman} Velocity Distribution on the Space at Different Time');