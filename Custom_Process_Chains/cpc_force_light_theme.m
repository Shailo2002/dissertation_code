function cpc_force_light_theme(fig)
% CPC_FORCE_LIGHT_THEME  Force a clean white-background, black-text
% appearance on every graphics object inside `fig` so that `print`
% produces thesis-ready output regardless of MATLAB's current theme
% (newer releases default to a dark theme that bleeds into saved JPEGs).
%
% Call this immediately before `print(fig, ...)`.

if nargin < 1; fig = gcf; end

set(fig, 'Color', 'w', 'InvertHardcopy', 'off');

% --- all axes (subplots, overlays) ---
all_ax = findall(fig, 'Type', 'axes');
for i = 1:numel(all_ax)
    ax = all_ax(i);
    set(ax, 'Color', 'w', ...
            'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
            'GridColor', [0.30 0.30 0.30], ...
            'MinorGridColor', [0.50 0.50 0.50]);
    try; ax.Title.Color  = 'k'; catch; end
    try; ax.XLabel.Color = 'k'; catch; end
    try; ax.YLabel.Color = 'k'; catch; end
    try; ax.ZLabel.Color = 'k'; catch; end
end

% --- every text object (annotations, sgtitle) ---
all_txt = findall(fig, 'Type', 'text');
for i = 1:numel(all_txt)
    try; set(all_txt(i), 'Color', 'k'); catch; end
end

% --- legends ---
all_leg = findall(fig, 'Type', 'legend');
for i = 1:numel(all_leg)
    try
        set(all_leg(i), 'TextColor', 'k', 'Color', 'w', ...
                        'EdgeColor', [0.6 0.6 0.6]);
    catch
    end
end

% --- colorbars (Type can be 'ColorBar' depending on release) ---
all_cb = [findall(fig, 'Type', 'colorbar'); findall(fig, 'Type', 'ColorBar')];
for i = 1:numel(all_cb)
    cb = all_cb(i);
    try; set(cb, 'Color', 'k'); catch; end
    try; cb.Label.Color = 'k'; catch; end
end
end
