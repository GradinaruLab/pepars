from plotly import graph_objs

from . import plotting


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

    for nucleotide in nucleotide_counts:

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
