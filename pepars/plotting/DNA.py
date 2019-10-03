import seaborn
import numpy
from plotly import graph_objs

from . import plotting
from ..utils import DNA as DNA_utils
from ..analysis import DNA as DNA_analysis


def plot_sequence_count_histogram(sequence_counts):

    seaborn.distplot(numpy.log10(sequence_counts))
    plt.show()


def plot_quality_score_distribution(
        quality_score_distribution,
        output_file_path=None,
        interactive=False,
        sample_name=None):

    num_reads = quality_score_distribution.sum(axis=1)[0]

    traces = []

    for position in range(quality_score_distribution.shape[0]):

        min_value = None
        q1_value = None
        median_value = None
        q3_value = None
        max_value = None
        cumulative_count = 0

        for quality_score in range(quality_score_distribution.shape[1]):

            quality_score_count = quality_score_distribution[position][
                quality_score]

            cumulative_count += quality_score_count

            if min_value is None and quality_score_count > 0:
                min_value = quality_score

            if q1_value is None and cumulative_count >= 0.25 * num_reads:
                q1_value = quality_score

            if median_value is None and cumulative_count >= 0.5 * num_reads:
                median_value = quality_score

            if q3_value is None and cumulative_count >= 0.75 * num_reads:
                q3_value = quality_score

            if quality_score_count > 0:
                max_value = quality_score

        y_values = [min_value, q1_value, median_value, median_value, q3_value,
                    max_value]

        trace = graph_objs.Box(
            y=y_values,
            name="%i" % (position + 1)
        )

        traces.append(trace)

    if sample_name is None:
        title = "Quality Score Distribution"
    else:
        title = "%s Quality Score Distribution" % sample_name

    layout = graph_objs.Layout(
        title=title,
        xaxis=dict(title="Position"),
        yaxis=dict(title="Quality score")
    )

    figure = graph_objs.Figure(data=traces, layout=layout)

    return plotting.generate_plotly_plot(figure,
                                         output_file_path=output_file_path,
                                         interactive=interactive)


def plot_nucleotide_prevalence_bar_chart(
        nucleotide_counts,
        output_file_path=None,
        interactive=False,
        sample_name=None):
    """
    Generates a stacked bar chart of the prevalence of different nucleotides.

    :param nucleotide_counts: A dictionary of arrays - an entry for each of
        A, C, G, T, N (optional). Each value in the array represents the count
        of that nucleotide in that position in the sequence
    :param output_file_path: The path to save the file, optional
    :param interactive: Whether this is an interactive plot (i.e. in a notebook)
    :param sample_name: The name of the sample to include in the plot title
    :return: Nothing
    """

    sequence_length = len(nucleotide_counts[
                              list(nucleotide_counts.keys())[0]])

    position_labels = list(range(sequence_length))

    traces = []

    for nucleotide in DNA_utils.get_nucleotides():

        trace = graph_objs.Bar(
            x=position_labels,
            y=nucleotide_counts[nucleotide],
            name=nucleotide
        )

        traces.append(trace)

    if "N" in nucleotide_counts:

        trace = graph_objs.Bar(
            x=position_labels,
            y=nucleotide_counts["N"],
            name="N"
        )

        traces.append(trace)

    if sample_name is None:
        title = "Nucleotide Distribution"
    else:
        title = "%s Nucleotide Distribution" % sample_name

    layout = graph_objs.Layout(
        barmode="stack",
        title=title
    )

    figure = graph_objs.Figure(data=traces, layout=layout)

    return plotting.generate_plotly_plot(figure,
                                         output_file_path=output_file_path,
                                         interactive=interactive)


def plot_amino_acid_bias(amino_acid_sequence_counts,
                         template_sequence=None,
                         amino_acid_probabilities=None,
                         sample_name=None,
                         ignore_counts=False,
                         biggest_value=None,
                         **kwargs):

    amino_acids = DNA_utils.get_amino_acids()

    if template_sequence is not None:
        unbiased_amino_acid_probabilities = \
            DNA_analysis.get_amino_acid_probabilities_from_template(
                template_sequence, allow_stop_codon=False)
    else:
        unbiased_amino_acid_probabilities = amino_acid_probabilities

    sequence_length = unbiased_amino_acid_probabilities.shape[0]

    amino_acid_position_counts = numpy.zeros(
        unbiased_amino_acid_probabilities.shape)

    for sequence, count in amino_acid_sequence_counts.items():

        if ignore_counts:
            count = 1

        for amino_acid_index, amino_acid in enumerate(sequence):
            amino_acid_position_counts[
                amino_acid_index, DNA_utils.AMINO_ACID_INDEX_MAP[
                    amino_acid]] += count

    sample_amino_acid_probabilities = \
        amino_acid_position_counts / \
        amino_acid_position_counts.sum(axis=1)[:, None]

    sample_amino_acid_probabilities_min = sample_amino_acid_probabilities[
        sample_amino_acid_probabilities > 0].min() / 2

    sample_amino_acid_probabilities[sample_amino_acid_probabilities == 0.0] = \
        sample_amino_acid_probabilities_min

    unbiased_amino_acid_probabilities_min = unbiased_amino_acid_probabilities[
        unbiased_amino_acid_probabilities > 0].min() / 2

    unbiased_amino_acid_probabilities[
        unbiased_amino_acid_probabilities == 0.0] = \
        unbiased_amino_acid_probabilities_min

    amino_acid_biases = numpy.log2(
        sample_amino_acid_probabilities / unbiased_amino_acid_probabilities)

    if biggest_value is None:
        biggest_value = max(abs(amino_acid_biases.max()),
                            abs(amino_acid_biases.min()))

    trace = graph_objs.Heatmap(
        z=amino_acid_biases,
        zmin=-biggest_value,
        zmax=biggest_value,
        zauto=False,
        x=amino_acids,
        y=[str(i) for i in range(1, sequence_length + 1)],
        colorbar=dict(
            title="Log2 Bias"
        ),
        colorscale=[
            [0.0, 'rgb(49,54,149)'],
            [0.1111111111111111, 'rgb(69,117,180)'],
            [0.2222222222222222, 'rgb(116,173,209)'],
            [0.3333333333333333, 'rgb(171,217,233)'],
            [0.4444444444444444, 'rgb(224,243,248)'],
            [0.5, 'rgb(255, 255, 255)'],
            [0.5555555555555556, 'rgb(254,224,144)'],
            [0.6666666666666666, 'rgb(253,174,97)'],
            [0.7777777777777778, 'rgb(244,109,67)'],
            [0.8888888888888888, 'rgb(215,48,39)'],
            [1.0, 'rgb(165,0,38)']
        ]
    )

    if sample_name is None:
        title = "Amino Acid Bias"
    else:
        title = "%s Amino Acid Bias" % sample_name

    layout = graph_objs.Layout(
        title=title
    )

    data = [trace]

    figure = graph_objs.Figure(data=data, layout=layout)

    plotting.generate_plotly_plot(figure, **kwargs)


OUTLINE_WIDTH = 2.5
OUTLINE_COLOR = "black"


def plot_signficant_amino_acid_biases(
        amino_acid_biases,
        p_values=None,
        biggest_value=None,
        sample_name=None,
        p_value_threshold=None,
        invert_outline=False,
        **kwargs
):
    amino_acid_biases = amino_acid_biases.transpose()

    if p_values is not None:
        p_values = p_values.transpose()

    amino_acids = DNA_utils.get_amino_acids()
    sequence_length = amino_acid_biases.shape[0]

    if biggest_value is None:
        biggest_value = max(abs(numpy.nanmax(amino_acid_biases)),
                            abs(numpy.nanmin(amino_acid_biases)))

    if invert_outline:
        xgap = 0
        ygap = 0
    else:
        xgap = 5
        ygap = 5

    traces = []

    if not invert_outline:

        for i in range(amino_acid_biases.shape[0]):
            for j in range(amino_acid_biases.shape[1]):

                if p_values[i, j] > p_value_threshold:
                    amino_acid_biases[i, j] = numpy.nan
                    continue

                if numpy.isnan(p_values[i, j]):
                    amino_acid_biases[i, j] = numpy.nan
                    continue

                outline_trace = graph_objs.Scatter(
                    x=[j + 1 - 0.43, j + 1 - 0.43],
                    y=[i + 1 - 0.43, i + 1 + 0.43],
                    mode="lines",
                    line={
                        "width": OUTLINE_WIDTH,
                        "color": OUTLINE_COLOR
                    },
                    showlegend=False,
                    hoverinfo="skip"
                )

                traces.append(outline_trace)

                outline_trace = graph_objs.Scatter(
                    x=[j + 1 + 0.43, j + 1 + 0.43],
                    y=[i + 1 - 0.43, i + 1 + 0.43],
                    mode="lines",
                    line={
                        "width": OUTLINE_WIDTH,
                        "color": OUTLINE_COLOR
                    },
                    showlegend=False,
                    hoverinfo="skip"
                )

                traces.append(outline_trace)

                outline_trace = graph_objs.Scatter(
                    x=[j + 1 - 0.43, j + 1 + 0.43],
                    y=[i + 1 - 0.43, i + 1 - 0.43],
                    mode="lines",
                    line={
                        "width": OUTLINE_WIDTH,
                        "color": OUTLINE_COLOR
                    },
                    showlegend=False,
                    hoverinfo="skip"
                )

                traces.append(outline_trace)

                outline_trace = graph_objs.Scatter(
                    x=[j + 1 + 0.43, j + 1 - 0.43],
                    y=[i + 1 + 0.43, i + 1 + 0.43],
                    mode="lines",
                    line={
                        "width": OUTLINE_WIDTH,
                        "color": OUTLINE_COLOR
                    },
                    showlegend=False,
                    hoverinfo="skip"
                )

                traces.append(outline_trace)
    else:

        crosshatch_x_values = []
        crosshatch_y_values = []

        for i in range(amino_acid_biases.shape[0]):
            for j in range(amino_acid_biases.shape[1]):

                if not numpy.isnan(p_values[i, j]) and p_values[
                        i, j] < p_value_threshold:
                    continue

                amino_acid_biases[i, j] = 0

                crosshatch_x_values.append(j + 1)
                crosshatch_y_values.append(i + 1)

        scatter = graph_objs.Scatter(
            x=crosshatch_x_values,
            y=crosshatch_y_values,
            mode="markers",
            marker={
                "symbol": "x-thin-open",
                "size": 35,
                "color": "gray"
            },
            hoverinfo="skip",
            showlegend=False
        )

        traces = [scatter]

    trace = graph_objs.Heatmap(
        z=amino_acid_biases,
        zmin=-biggest_value,
        zmax=biggest_value,
        zauto=False,
        x=[i + 1 for i in range(len(amino_acids))],
        xgap=xgap,
        ygap=ygap,
        y=[str(i) for i in range(1, sequence_length + 1)],
        colorbar=dict(
            title="Bias Level"
        ),
        colorscale=[
            [0.0, 'rgb(49,54,149)'],
            [0.1111111111111111, 'rgb(69,117,180)'],
            [0.2222222222222222, 'rgb(116,173,209)'],
            [0.3333333333333333, 'rgb(171,217,233)'],
            [0.4444444444444444, 'rgb(224,243,248)'],
            [0.5, 'rgb(255, 255, 255)'],
            [0.5555555555555556, 'rgb(254,224,144)'],
            [0.6666666666666666, 'rgb(253,174,97)'],
            [0.7777777777777778, 'rgb(244,109,67)'],
            [0.8888888888888888, 'rgb(215,48,39)'],
            [1.0, 'rgb(165,0,38)']
        ]
    )

    if sample_name is None:
        title = "Amino Acid Bias"
    else:
        title = "%s Amino Acid Bias" % sample_name

    layout = graph_objs.Layout(
        title=title,
        xaxis={
            "ticktext": amino_acids,
            "tickvals": [i + 1 for i in range(len(amino_acids))],
            "title": "Amino Acid",
            "range": [0.5, len(amino_acids) + 0.5],
            "gridcolor": "rgba(0,0,0,0)"
        },
        yaxis={
            "title": "Position",
            "range": [0.5, sequence_length + 0.5],
            "gridcolor": "rgba(0,0,0,0)"
        },
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)"
    )

    data = [trace]
    data.extend(traces)

    figure = graph_objs.Figure(data=data, layout=layout)

    plotting.generate_plotly_plot(figure, **kwargs)

    return figure
