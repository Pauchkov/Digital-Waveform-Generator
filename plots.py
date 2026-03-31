import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Wedge


RADIUS = 2.0
WIDTH = 0.6
MARGIN = 0.5
WOM1_CENTER = (0, 2.5)
WOM2_CENTER = (0, -2.5)
CHANNEL_ANGLES = [
    (90, 135),
    (135, 180),
    (180, 225),
    (225, 270),
    (270, 315),
    (315, 360),
    (0, 45),
    (45, 90),
]


def plot_wom_channels(wom1_values, wom2_values, cmap_name='viridis', label = 'WOM Channel Values'):
    if len(wom1_values) != 8 or len(wom2_values) != 8:
        raise ValueError('Expected 8 values for WOM1 and 8 values for WOM2')

    all_values = list(wom1_values) + list(wom2_values)
    vmin = min(all_values)
    vmax = max(all_values)
    if vmin == vmax:
        vmin -= 0.5
        vmax += 0.5

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    fig, (ax, cax) = plt.subplots(
        ncols=2,
        figsize=(5.5, 8),
        gridspec_kw={'width_ratios': [1, 0.06]},
        constrained_layout=True,
    )

    WOM1_channels = [
        (a1, a2, value) for (a1, a2), value in zip(CHANNEL_ANGLES, wom1_values)
    ]

    WOM2_channels = [
        (a1, a2, value) for (a1, a2), value in zip(CHANNEL_ANGLES, wom2_values)
    ]

    for a1, a2, value in WOM1_channels:
        ax.add_patch(Wedge(WOM1_CENTER, RADIUS, a1, a2, width=WIDTH,
                           facecolor=cmap(norm(value)), edgecolor='black'))

    for a1, a2, value in WOM2_channels:
        ax.add_patch(Wedge(WOM2_CENTER, RADIUS, a1, a2, width=WIDTH,
                           facecolor=cmap(norm(value)), edgecolor='black'))

    x_min = min(WOM1_CENTER[0], WOM2_CENTER[0]) - RADIUS - MARGIN
    x_max = max(WOM1_CENTER[0], WOM2_CENTER[0]) + RADIUS + MARGIN
    y_min = min(WOM1_CENTER[1], WOM2_CENTER[1]) - RADIUS - MARGIN
    y_max = max(WOM1_CENTER[1], WOM2_CENTER[1]) + RADIUS + MARGIN

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal', adjustable='box')
    ax.set_anchor('C')
    ax.axis('off')

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array(all_values)
    fig.colorbar(sm, cax=cax, label='Channel value')

    plt.suptitle(label)

    plt.show()

plot_wom_channels([1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], 'plasma', 'Example WOM Channel Values')
