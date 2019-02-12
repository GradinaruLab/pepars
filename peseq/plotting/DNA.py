import seaborn
from plotly import graph_objs

from . import plotting


def plot_sequence_count_histogram(sequence_counts):

    seaborn.distplot(numpy.log10(sequence_counts))
    plt.show()


def plot_nucleotide_prevalence_bar_chart(
        nucleotide_counts,
        output_file_path=None,
        interactive=False):
    """
    Generates a stacked bar chart of the prevalence of different nucleotides.

    :param nucleotide_counts: A dictionary of arrays - an entry for each of
        A, C, G, T, N (optional). Each value in the array represents the count
        of that nucleotide in that position in the sequence
    :param output_file_path: The path to save the file, optional
    :param interactive: Whether this is an interactive plot (i.e. in a notebook)
    :return: Nothing
    """

    sequence_length = len(nucleotide_counts[
                              list(nucleotide_counts.keys())[0]])

    position_labels = list(range(sequence_length))

    traces = []

    for nucleotide in sorted(nucleotide_counts):

        trace = graph_objs.Bar(
            x=position_labels,
            y=nucleotide_counts[nucleotide],
            name=nucleotide
        )

        traces.append(trace)

    layout = graph_objs.Layout(
        barmode="stack"
    )

    figure = graph_objs.Figure(data=traces, layout=layout)

    return plotting.generate_plotly_plot(figure,
                                         output_file_path=output_file_path,
                                         interactive=interactive)

def plot_significance_z_scores(
        z_scores,
        output_file_path=None,
        interactive=False,
        colorscale="RdBu",
        min_value=None,
        max_value=None):

    if min_value is None:
        min_value = -max(z_scores.max().max(), -z_scores.min().min())
    if max_value is None:
        max_value = max(z_scores.max().max(), -z_scores.min().min())

    figure_traces = []

    heatmap = graph_objs.Heatmap(
        x=z_scores.index,
        y=z_scores.columns,
        z=z_scores.values.T,
        zmin=min_value,
        zmax=max_value,
        zauto=False,
        colorscale=colorscale
    )

    figure_traces.append(heatmap)

    layout = graph_objs.Layout(
        title="Amino Acid Presence",
        xaxis=dict(title="Amino Acid"),
        yaxis=dict(title="Position")
    )

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    return generate_plotly_plot(figure, output_file_path=output_file_path,
                                interactive=interactive)
