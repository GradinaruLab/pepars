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
                   log_scale=False):

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

    if log_scale:
        layout_parameters["yaxis"]["type"] = "log"

    layout = graph_objs.Layout(layout_parameters)

    figure = graph_objs.Figure(data=figure_traces, layout=layout)

    generate_plotly_plot(figure,
                         output_file_path=output_file_path,
                         interactive=interactive)
