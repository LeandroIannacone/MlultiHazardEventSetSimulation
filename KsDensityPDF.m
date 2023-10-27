
% Mainshock-Mainshock
gridx1 = min(z.mm(:,1)):.05:max(z.mm(:,1));
gridx2 = min(z.mm(:,2)):.05:max(z.mm(:,2));
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.mm(:,1),z.mm(:,2)];
figure
ksdensity(x,xi);
xlim([4.45 6])
ylim([4.45 6])
view(0,90)
xlabel('1st hazard severity, $m_m$','Interpreter','latex','FontSize',18)
ylabel('2nd hazard severity, $m_m$','Interpreter','latex','FontSize',18)
colorbar
colormap(flipud(copper))
shading interp

% Mainshock-Aftershock
gridx1 = min(z.ma(:,1)):.05:max(z.ma(:,1));
gridx2 = min(z.ma(:,2)):.05:max(z.ma(:,2));
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.ma(:,1),z.ma(:,2)];
figure
ksdensity(x,xi);
xlim([4.45 6])
ylim([4.45 6])
view(0,90)
xlabel('1st hazard severity, $m_m$','Interpreter','latex','FontSize',18)
ylabel('2nd hazard severity, $m_a$','Interpreter','latex','FontSize',18)
colorbar
colormap(flipud(copper))
shading interp

% Aftershock - Mainshock
gridx1 = min(z.am(:,1)):.05:max(z.am(:,1));
gridx2 = min(z.am(:,2)):.05:max(z.am(:,2));
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.am(:,1),z.am(:,2)];
figure
ksdensity(x,xi);
title("Aftershock - Mainshock")

% Mainshock - Rain
gridx1 = linspace(min(z.mr(:,1)),max(z.mr(:,1)),100);
gridx2 = linspace(min(z.mr(:,2)),max(z.mr(:,2)),100);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.mr(:,1),z.mr(:,2)];
figure
ksdensity(x,xi);
xlim([4.45 6])
ylim([0 100])
view(0,90)
xlabel('1st hazard severity, $m_m$','Interpreter','latex','FontSize',18)
ylabel('2nd hazard severity, $m_r$ [mm/hr]','Interpreter','latex','FontSize',18)
colorbar
colormap(flipud(copper))
shading interp

% Aftershock - Rain
gridx1 = min(z.ar(:,1)):.05:max(z.ar(:,1));
gridx2 = min(z.ar(:,2)):.05:max(z.ar(:,2));
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.ar(:,1),z.ar(:,2)];
figure
ksdensity(x,xi);
title("Aftershock - Rain")

% Rain - Rain
gridx1 = linspace(min(z.rr(:,1)),max(z.rr(:,1)),100);
gridx2 = linspace(min(z.rr(:,2)),max(z.rr(:,2)),100);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.rr(:,1),z.rr(:,2)];
figure
ksdensity(x,xi);
title("Rain - Rain")

% Rain - Mainshock
gridx1 = min(z.rm(:,1)):.05:max(z.rm(:,1));
gridx2 = min(z.rm(:,2)):.05:max(z.rm(:,2));
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.rm(:,1),z.rm(:,2)];
figure
ksdensity(x,xi);
title("Rain - Mainshock")

% Rain - Aftershock
gridx1 = min(z.ra(:,1)):.05:max(z.ra(:,1));
gridx2 = min(z.ra(:,2)):.05:max(z.ra(:,2));
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
x = [z.ra(:,1),z.ra(:,2)];
figure
ksdensity(x,xi);
title("Rain - Aftershock")