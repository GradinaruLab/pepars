import numpy

from plotly import offline as plotly_offline
from plotly import graph_objs


def init_notebook_mode():
    plotly_offline.init_notebook_mode()


def generate_plotly_plot(figure, output_file_path=None, interactive=False):

    if output_file_path is not None:
        if len(output_file_path) < 5 or output_file_path[-5:] != ".html":
            output_file_path += ".html"

    if interactive:
        plotly_offline.iplot(figure)

        if output_file_path:
            plotly_offline.iplot(figure, filename=output_file_path)

    elif output_file_path:
        plotly_offline.plot(figure, filename=output_file_path, auto_open=False)


def plot_histogram(values,
                   interactive=False,
                   output_file_path=None,
                   title="Histogram",
                   x_axis_title=None,
                   y_axis_title="Count",
                   num_bins=None,
                   log_scale=False,
                   x_range=None,
                   y_range=None):

    figure_traces = []

    figure_parameters = {
        "x": values
    }

    if num_bins is not None:
        figure_parameters["nbinsx"] = num_bins

    histogram = graph_objs.Histogram(figure_parameters)

    figure_traces.append(histogram)

    layout_parameters = {
        "title": title,
        "xaxis": {
        },
        "yaxis": {
            "title": y_axis_title
        },
        "hovermode": "closest"
    }

    if x_axis_title is not None:
        layout_parameters["xaxis"]["title"] = x_axis_title

    if x_range is not None:
        layout_parameters["xaxis"]["range"] = x_range
    if y_range is not None:
        layout_parameters["yaxis"]["range"] = y_range

    if log_scale:
        layout_parameters["yaxis"]["type"] = "log"

    layout = graph_objs.Layout(layout_parameters)

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    generate_plotly_plot(figure,
                         output_file_path=output_file_path,
                         interactive=interactive)


def plot_count_distribution(sample_counts, **kwargs):

    if isinstance(sample_counts, list):
        sample_counts = {"": sample_counts}

    figure_traces = []

    for sample, counts in sample_counts.items():

        counts = sorted(counts)

        values, counts = numpy.unique(counts, return_counts=True)

        cumulative_counts = numpy.cumsum(counts)

        percentiles = cumulative_counts / sum(counts)

        scatter = graph_objs.Scatter(
            x=values,
            y=percentiles,
            name=sample
        )

        figure_traces.append(scatter)

    if len(sample_counts) == 1:
        title = "%s Sequence Count Distribution" % list(sample_counts)[0]
    else:
        title = "Sequence Count Distribution"

    layout = graph_objs.Layout(
        hovermode="closest",
        xaxis=dict(
            type="log",
            title="Sequence Count"
        ),
        yaxis=dict(
            range=[0, 1],
            title="Probability"
        ),
        title=title
    )

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    return generate_plotly_plot(figure, **kwargs)


def plot_scatter(x_values,
                 y_values,
                 output_file_path=None,
                 interactive=False,
                 text_labels=None):

    figure_traces = []

    if text_labels is None:
        scatter = graph_objs.Scatter(
            x=list(x_values),
            y=list(y_values),
            mode="markers"
        )
    else:
        scatter = graph_objs.Scatter(
            x=list(x_values),
            y=list(y_values),
            mode="markers",
            hoverinfo="text",
            text=text_labels
        )

    figure_traces.append(scatter)

    layout = graph_objs.Layout(
        hovermode="closest"
    )

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    return generate_plotly_plot(figure, output_file_path=output_file_path,
                                interactive=interactive)


def plot_bar_chart(values,
                   condition_names,
                   group_names=None,
                   errors=None,
                   interactive=False,
                   output_file_path=None,
                   title="Bar Chart"):
    """
    Values can be either a 1-dimensional array or 2-dimensional array.
    If 1-dimensional, the values should correspond to different conditions and
    will be placed side by side on the bar chart.
    If 2-dimensional, the first dimension corresponds to different conditions,
    the 2nd dimension to different groups
    """

    values = numpy.array(values)
    if errors is not None:
        errors = numpy.array(errors)
        if errors.shape != values.shape:
            raise ValueError("Error and values must be of the same size")

    if condition_names is None or len(condition_names) != values.shape[0]:
        raise ValueError("Must specify a condition name for each value")

    figure_traces = []

    if len(values.shape) == 1 or values.shape[1] == 1:
        values = numpy.reshape(values, newshape=(values.shape[0], 1))
    else:
        if group_names is None:
            raise ValueError("Must specify group names for 2D bar charts")

    for condition_index, condition_values in enumerate(values):

        figure_parameters = {
            "y": values[condition_index, :],
            "name": condition_names[condition_index]
        }

        if group_names is not None:
            figure_parameters["x"] = group_names

        if errors is not None:
            figure_parameters["error_y"] = {
                "array": errors[condition_index, :],
                "type": "data",
                "visible": True
            }

        trace = graph_objs.Bar(figure_parameters)

        figure_traces.append(trace)

    layout_parameters = {
        "barmode": "group",
        "title": title,
        "xaxis": {
        },
        "hovermode": "closest"
    }

    if group_names is None:
        layout_parameters["xaxis"]["visible"] = False

    layout = graph_objs.Layout(layout_parameters)

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    generate_plotly_plot(figure,
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

    return generate_plotly_plot(figure,
                                output_file_path=output_file_path,
                                interactive=interactive)