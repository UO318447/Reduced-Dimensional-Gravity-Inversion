function dibuja_anomali(gobs, gpre, x, y)

% Plot observed and predicted gravity anomalies using a common color scale.

figure('Position', [100, 100, 600, 800])

% Common color limits
vmin = min([gobs(:); gpre(:)]);
vmax = max([gobs(:); gpre(:)]);

% Observed gravity anomaly
g_medida = reshape(gobs, size(x));

subplot(2, 1, 1)
contourf(x, y, g_medida)
colorbar
clim([vmin vmax])
axis('equal')
ylabel('Y (m)', 'FontWeight', 'bold')
title('g^{obs} (\muGal)')
set(gca, 'FontSize', 12)

% Predicted gravity anomaly
grave = reshape(gpre, size(x));

subplot(2, 1, 2)
contourf(x, y, grave)
colorbar
clim([vmin vmax])
axis('equal')
xlabel('X (m)', 'FontWeight', 'bold')
ylabel('Y (m)', 'FontWeight', 'bold')
title('g^{pre} (\muGal)')
set(gca, 'FontSize', 12)

% Apply consistent font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);

end