import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Wedge
from pathlib import Path
import ROOT


ROOT.TH1.AddDirectory(False)

RADIUS = 2.0
WIDTH = 0.6
MARGIN = 0.5
WOM1_CENTER = (0, 2.5)
WOM2_CENTER = (0, -2.5)
N_SAMPLES = 1024
TIME_WINDOW_NS = 320.0
TIME0_BINS = 128
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


def save_hist_to_root(histogram, output_path, hist_name=None, mode="RECREATE"):
    output_path = Path(output_path)
    output_file = ROOT.TFile.Open(str(output_path), mode)
    if not output_file or output_file.IsZombie():
        raise OSError(f"Unable to open output ROOT file: {output_path}")

    target_name = hist_name or histogram.GetName()
    histogram_copy = histogram.Clone(target_name)
    histogram_copy.SetDirectory(output_file)

    output_file.cd()
    written = histogram_copy.Write(target_name, ROOT.TObject.kOverwrite)
    output_file.Close()

    if written <= 0:
        raise OSError(f"Unable to write histogram '{target_name}' to {output_path}")

    return output_path


def plot_wom_channels(wom1_values, wom2_values, cmap_name='plasma', label = 'WOM Channel Values'):
    if len(wom1_values) != 8 or len(wom2_values) != 8:
        raise ValueError('Expected 8 values for WOM1 and 8 values for WOM2')

    all_values = list(wom1_values) + list(wom2_values)
    vmin = min(all_values)
    vmax = max(all_values)

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



# Load files
f = ROOT.TFile.Open("mutestfile.root")
main_hist3D = f.Get("evt_quadrant_time")
main_hist3D.SetDirectory(0)
f.Close()


norm1_pe = []
f = ROOT.TFile.Open("norm_1pe.root")
for ch in range(64):
    graph = f.Get(f"norm_1pe__{ch:02d}")
    if not graph:
        raise KeyError(f"Unable to find norm_1pe__{ch:02d} in norm_1pe.root")
    norm1_pe.append(graph.Clone(f"norm_1p_ch{ch}"))
f.Close()

DC_values = []
with Path("DC_results.txt").open() as f:
    for line in f:
            DC_values.append(tuple(map(float, line.split()[:4])))


print("DC values:", DC_values)

norm1_pe_hist = []
template_integrals = []
for ch, graph in enumerate(norm1_pe):
    x_values = graph.GetX()
    y_values = graph.GetY()
    n_points = graph.GetN()

    template_hist = ROOT.TH1F(
        f"norm_1pe_hist_ch{ch:02d}",
        f"1pe template channel {ch:02d}",
        n_points,
        x_values[0],
        x_values[n_points - 1],
    )
    template_hist.SetDirectory(0)

    for point_index in range(1, n_points):
        repeats = int(max(0.0, y_values[point_index] * 100.0))
        for _ in range(repeats):
            template_hist.Fill(x_values[point_index])

    template_hist.Scale(1.0 / 100.0)
    norm1_pe_hist.append(template_hist)

    max_bin = template_hist.GetMaximumBin()
    left_bin = max(1, max_bin - 5 * N_SAMPLES // 320)
    right_bin = min(template_hist.GetNbinsX(), max_bin + 9 * N_SAMPLES // 320)
    template_integral = 0.0
    for bin_index in range(left_bin, right_bin + 1):
        template_integral += template_hist.GetBinContent(bin_index)

    if template_integral <= 0.0:
        raise RuntimeError(f"Template integral is zero for channel {ch}")

    template_integrals.append(template_integral)

digitized_waveforms = []
n_events = 1
n_channels = min(main_hist3D.GetNbinsY(), len(norm1_pe_hist), len(DC_values))

for evt in range(n_events):
    event_bin = 1 
    event_waveforms = []

    for ch in range(n_channels):
        print("[INFO] Channel", ch, "Event", evt)
        waveform_hist = ROOT.TH1F(
            f"wf_evt{evt:02d}_ch{ch:02d}",
            "Digitized waveform without dark count or cross-talk",
            N_SAMPLES,
            0.0,
            TIME_WINDOW_NS,
        )
        waveform_hist.SetDirectory(0)
        waveform_hist.GetXaxis().SetTitle("Time [ns]")
        waveform_hist.GetYaxis().SetTitle("Amplitude [mV]")
        waveform_hist.SetLineWidth(2)

        template = norm1_pe_hist[ch]
        template_scale = DC_values[ch][0] / template_integrals[ch]

        for time_bin in range(N_SAMPLES):
            photons = main_hist3D.GetBinContent(event_bin, ch + 1, time_bin + 1)
            if photons <= 0.0:
                continue

            visible_template_len = min(template.GetNbinsX(), N_SAMPLES - TIME0_BINS - time_bin)
            if visible_template_len <= 0:
                continue

            for template_index in range(visible_template_len):
                destination_bin = TIME0_BINS + time_bin + template_index + 1
                waveform_hist.AddBinContent(
                    destination_bin,
                    photons * template_scale * template.GetBinContent(template_index + 1),
                )

        event_waveforms.append(waveform_hist)

    digitized_waveforms.append(event_waveforms)

print(f"Digitized {len(digitized_waveforms)} events for {n_channels} channels without dark count or cross-talk.")

output_file = ROOT.TFile.Open("digitized_waveforms.root", "RECREATE")
if not output_file or output_file.IsZombie():
    raise OSError("Unable to create output ROOT file: digitized_waveforms.root")

for event_waveforms in digitized_waveforms:
    for waveform_hist in event_waveforms:
        output_file.cd()
        written = waveform_hist.Write(waveform_hist.GetName(), ROOT.TObject.kOverwrite)
        if written <= 0:
            output_file.Close()
            raise OSError(f"Unable to write histogram '{waveform_hist.GetName()}'")

output_file.Close()
print("Saved all digitized channels to digitized_waveforms.root")

plot_wom_channels([1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], 'plasma', 'Example WOM Channel Values')
plot_wom_channels([1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], 'plasma', 'Example WOM Channel Values 2')
plt.show()
