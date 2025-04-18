% --- Article: Time-reversible synchronization of analog and digital chaotic systems
% --- Authors: Artur Karimov, Vyacheslav Rybin, Ivan Babkin, Timur Karimov, Veronika Ponomareva and Denis Butusov
% --- Time-reversible synchronization of Sprott Case B system

clear all;
% Integration step
h = 1.0e-2;
% Computational iterations
CT_iter = ceil(250/h);
% Transient iterations
TT_iter = ceil(250/h);
% Initial conditions
y0 = [0.1;0.1;0.1];
% Number of synchronization iterations
nglobal = 200;
% Number of window points
ns = ceil(1/h);
% Forward synchronization coefficient
K0 = [3 0 0];
% Backward synchronization coefficient
Kb = [3 0 0];
% Decimation of master array (for calculation acceleration)
decimator = 20;

% All points
ym = zeros(3,CT_iter+TT_iter + ns);

% Calculation of master signal
X = y0;
for i = 1:CT_iter+TT_iter + ns
    h1 = h*0.5; h2 = h*0.5;
    a1 = 0.5;
    a2 = 1.6;
    a3 = -1.5;
    a4 = 1;
    a5 = - 1;
    X(1) = X(1) + h1 *  (a1*X(2)*X(3));
    X(2) = X(2) + h1 *  (a2*X(1) + a3*X(2));
    X(3) = X(3) + h1 *  (a4 + a5*X(1)*X(2));
    X(3) = X(3) + h2 *  (a4 + a5*X(1)*X(2));
    X(2) = (X(2) + h2 * (a2*X(1)))/(1 - a3*h2);
    X(1) = (X(1) + h2 * (a1*X(2)*X(3)));
    ym(:, i) = X;

end

% Removing transient points
ym = ym(:,TT_iter:end);



N = length(ym(1,:)) - ns;
ym_starts = zeros(3, N);
err = zeros(1, N);
niters =  zeros(1, N);

hw = waitbar(0,'Please wait...');

for nstart = 1:decimator:N
    waitbar(nstart/N,hw,'Processing...');

    nts = nstart;
    ym_starts(:, nstart) = ym(:, nts);
    ymaster = ym(:, nts:nts + ns - 1);
    y0 = ym(:,1)+0.01;
    %y0 = [0;0;0];
    [err(nstart), niters(nstart)] = caseB_synchro(y0, h, ymaster, ns, nglobal,K0,Kb);
end
close(hw);
err(isnan(err)) = 1000;
err10 = log10(err);
err10 = err10(1:decimator:end);
ym_starts = ym_starts(:,1:decimator:end);

figure
cm1 = [linspace(0.4796,0.8,100)' linspace(0.01583,0.8,100)' linspace(0.01055,0.8,100)'];
cm = [colormap(turbo(450)); cm1];
LineWdth = 1.5;
alpha = 1;
set(gcf, 'Position', [350 50 700 500]);
set(gcf,'renderer','Painters');
ax(1) = axes;
ax1pos = get(gca,'position');
ax1pos(3) = ax1pos(3)*0.7; %length x
ax1pos(4) = 1*ax1pos(4); %length y
set(gca,'position', ax1pos);
ym_starts = [nan nan nan; ym_starts'; nan nan nan]';
err10 = [nan err10() nan];
p1 = patch(ym_starts(1,1:end),ym_starts(2,1:end),ym_starts(3,1:end),err10(1:end),'EdgeColor','interp','FaceColor','interp','LineWidth',LineWdth);
p1.EdgeAlpha = alpha;
view(-37.5,30); grid on;
set(gca,'FontSize',18);
set(gca, 'XDir', 'normal'); hold on;
set(gca, 'YDir', 'reverse'); hold on;
set(gca,'TickLabelInterpreter','latex','FontSize',20);
zlabel('$z$','interpreter','latex','FontSize',24)
ylabel('$y$','interpreter','latex','FontSize',24)
xlabel('$x$','interpreter','latex','FontSize',24)
colormap(cm);
c = colorbar('eastoutside');
clim([-12 2]);
c.Label.String = '$log_{10}(||Error||)$';
c.Label.Interpreter = 'latex';
set(c,'TickLabelInterpreter','latex','FontSize',20);
c.Location = 'eastoutside';
c.Position = c.Position .* [1.27 1 1.1 1];
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
c.Ruler.TickLabelFormat = '$%g$';


% Calcualtion of synchronization
function [reachederror, niter] = caseB_synchro(y0, h, ym, ns,nglobal,K0,Kb)
%slave CD simulation

Ys = y0(1:3);
N = [0, 0, 0];
ys = zeros(3,ns);
ys(:,1) = Ys;

for j = 1:nglobal
    %forward
    K = K0;
    for i = 2:ns
        %slave
        X = Ys;
        ex = (ym(1, i-1) - X(1));
        ey = (ym(2, i-1) - X(2));
        ez = (ym(3, i-1) - X(3));

        N(1) = K(1)*ex;
        N(2) = K(2)*ey;
        N(3) = K(3)*ez;

        if j == nglobal && i == ns
            reachederror = sqrt(ex*ex + ey*ey + ez*ez);
        end
        h1 = h*0.5; h2 = h*0.5;


        a1 = 0.5;
        a2 = 1.6;
        a3 = -1.5;
        a4 = 1;
        a5 = - 1;

        X(1) = X(1) + h1 *  (N(1) + a1*X(2)*X(3));
        X(2) = X(2) + h1 *  (N(2) + a2*X(1) + a3*X(2));
        X(3) = X(3) + h1 *  (N(3) + a4 + a5*X(1)*X(2));
        X(3) = X(3) + h2 *  (N(3) + a4 + a5*X(1)*X(2));
        X(2) = (X(2) + h2 * (N(2) + a2*X(1)))/(1 - a3*h2);
        X(1) = (X(1) + h2 * (N(1) + a1*X(2)*X(3)));

        Ys = X;
        ys(:, i) = Ys;

        if isnan(X)
            reachederror = nan;
            break;
        end
    end

    if j < nglobal
        %backward
        K = -Kb;
        for i = ns-1:-1:1

            X = Ys;
            N(1) = K(1)*(ym(1, i+1) - X(1));
            N(2) = K(2)*(ym(2, i+1) - X(2));
            N(3) = K(3)*(ym(3, i+1) - X(3));
            h1 = -h*0.5; h2 = -h*0.5;

            a1 = 0.5;
            a2 = 1.6;
            a3 = -1.5;
            a4 = 1;
            a5 = - 1;

            X(1) = X(1) + h1 *  (N(1) + a1*X(2)*X(3));
            X(2) = X(2) + h1 *  (N(2) + a2*X(1) + a3*X(2));
            X(3) = X(3) + h1 *  (N(3) + a4 + a5*X(1)*X(2));
            X(3) = X(3) + h2 *  (N(3) + a4 + a5*X(1)*X(2));
            X(2) = (X(2) + h2 * (N(2) + a2*X(1)))/(1 - a3*h2);
            X(1) = (X(1) + h2 * (N(1) + a1*X(2)*X(3)));

            Ys = X;
            if isnan(X)
                reachederror = nan;
                break;
            end
        end
    end
end
niter = j;
end





