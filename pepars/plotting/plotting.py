import numpy

from plotly import offline as plotly_offline
from plotly import graph_objs

INTERACTIVE_MODE = False


def init_notebook_mode():
    plotly_offline.init_notebook_mode()
    global INTERACTIVE_MODE
    INTERACTIVE_MODE = True


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
                   y_range=None,
                   x_axis_log_scale=False,
                   y_axis_log_scale=False):

    is_integer_values = True

    # There has got to be a better way to do this...
    # Sample 10% of the values to see if they're all integers
    num_values_to_sample = int(numpy.ceil(numpy.sqrt(len(values))))
    for i in range(num_values_to_sample):
        random_index = numpy.random.randint(0, len(values))
        if values[random_index] != int(values[random_index]):
            is_integer_values = False
            break

    figure_traces = []

    if num_bins is None:
        # Estimate the number of bins using the Rice Rule
        num_bins = int(2 * numpy.ceil(numpy.power(len(values), 1 / 3)))

        if is_integer_values and num_bins > max(values):
            num_bins = int(numpy.ceil(max(values)))

    if is_integer_values:

        min_value = min(values)
        max_value = max(values)

        value_range = max_value - min_value
        bin_size = numpy.ceil(value_range / num_bins)
        min_bin = numpy.floor(min_value / bin_size) * bin_size
        max_bin = numpy.ceil(max_value / bin_size) * bin_size
        bins = numpy.arange(min_bin, max_bin+bin_size, bin_size)
    else:
        bins=num_bins

    y, x = numpy.histogram(values, bins=bins)

    histogram = graph_objs.Bar(
        x=x,
        y=y,
    )

    figure_traces.append(histogram)

    layout_parameters = {
        "title": title,
        "xaxis": {
        },
        "yaxis": {
            "title": y_axis_title
        },
        "hovermode": "closest",
        "bargap": 0
    }

    if x_axis_title is not None:
        layout_parameters["xaxis"]["title"] = x_axis_title

    if x_range is not None:
        layout_parameters["xaxis"]["range"] = x_range
    if y_range is not None:
        layout_parameters["yaxis"]["range"] = y_range

    if log_scale:
        layout_parameters["yaxis"]["type"] = "log"

    if x_axis_log_scale:
        layout["xaxis"]["type"] = "log"

    if y_axis_log_scale:
        layout["yaxis"]["type"] = "log"

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
                 text_labels=None,
                 trace_names=None,
                 title=None,
                 x_range=None,
                 y_range=None,
                 x_axis_title=None,
                 y_axis_title=None,
                 x_axis_log_scale=False,
                 y_axis_log_scale=False,
                 x_intersect=None,
                 y_intersect=None,
                 mode="markers"):

    figure_traces = []

    # If the x values have a 2nd dimension, we have multiple traces
    try:
        len(x_values[0])

        for series_index in range(len(x_values)):

            if trace_names is not None:
                trace_name = trace_names[series_index]
            else:
                trace_name = "%i" % (series_index + 1)

            if text_labels is None:
                scatter = graph_objs.Scatter(
                    x=list(x_values[series_index]),
                    y=list(y_values[series_index]),
                    mode=mode,
                    name=trace_name
                )
            else:
                scatter = graph_objs.Scatter(
                    x=list(x_values[series_index]),
                    y=list(y_values[series_index]),
                    mode=mode,
                    hoverinfo="text",
                    text=text_labels,
                    name=trace_name
                )
            figure_traces.append(scatter)

    except TypeError:

        if text_labels is None:
            scatter = graph_objs.Scatter(
                x=list(x_values),
                y=list(y_values),
                mode=mode
            )
        else:
            scatter = graph_objs.Scatter(
                x=list(x_values),
                y=list(y_values),
                mode=mode,
                hoverinfo="text",
                text=text_labels
            )

        figure_traces.append(scatter)

    layout = graph_objs.Layout(
        hovermode="closest"
    )

    layout["xaxis"] = {}
    layout["yaxis"] = {}

    if x_range:
        layout["xaxis"]["range"] = x_range

    if y_range:
        layout["yaxis"]["range"] = y_range

    if x_axis_title:
        layout["xaxis"]["title"] = x_axis_title

    if y_axis_title:
        layout["yaxis"]["title"] = y_axis_title

    if title:
        layout["title"] = title

    if x_axis_log_scale:
        layout["xaxis"]["type"] = "log"

    if y_axis_log_scale:
        layout["yaxis"]["type"] = "log"

    if y_intersect:
        extra_trace = graph_objs.Scatter(
            x=[0, max(x_values)],
            y=[y_intersect, y_intersect],
            mode="markers+lines"
        )

        figure_traces.append(extra_trace)

    if x_intersect:
        extra_trace = graph_objs.Scatter(
            x=[x_intersect, x_intersect],
            y=[0, max(y_values)],
            mode="markers+lines"
        )

        figure_traces.append(extra_trace)

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    return generate_plotly_plot(figure, output_file_path=output_file_path,
                                interactive=interactive)


def plot_bar_chart(values,
                   condition_names,
                   group_names=None,
                   errors=None,
                   interactive=False,
                   output_file_path=None,
                   title="Bar Chart",
                   y_axis_title=None):
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

    if y_axis_title:
        layout_parameters["yaxis"] = {}
        layout_parameters["yaxis"]["title"] = y_axis_title

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
