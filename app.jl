using Printf
using Random

using DifferentialEquations

using Dash
using PlotlyJS

println("Initializing Dash app...")

# define a system to be simulated
# negative values in solution when parameter values are altered
function μ_func(du, u, p, t)
    SC, P, D = u
    λ, a₁, a₂, a₃, a₄, a₅, a₆, d₁, d₂, d₃ = p

    # stem
    du[1] = λ + a₁ * SC - a₃ * SC - d₁ * SC
    # progenitor
    du[2] = a₂ * SC + a₄ * P + 2 * a₃ * SC - a₆ * P - d₂ * P
    # differentiated
    du[3] = a₅ * P + 2a₆ * P - d₃ * D
end

function σ_func(du, u, p, t)
  du[1] = 10
  du[2] = 10
  du[3] = 10
end

# TODO: make parameter values and ODE solution local in user session
cell_types = ["Stem cells", "Progenitor cells", "Red blood cells"]
num_cell_types = length(cell_types)
ss_params = [20.0, 0.33, 0.3, 0.35, 0.33, 0.33, 0.33, 0.0, 0.4, 0.5]
ss_params_min = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.24]
ss_param_names = ["lambda", "a1", "a2", "a3", "a4", "a5", "a6", "d1", "d2",
                  "d3"]
num_ss_params = length(ss_params)
ss_param2idx = Dict([ss_param_names[i] => i for i = 1:num_ss_params])
time_span = (0.0, 100.0)
init_condition = [1000.0, 100.0, 100.0]

# solve the system given current parameter values
function get_ss_solution()
    ss_sde_problem = SDEProblem(μ_func, σ_func, init_condition, time_span,
                                ss_params)
    solve(ss_sde_problem, SOSRI())
end

ss_solution = get_ss_solution()

# determine hero status based on number of cells
function get_hero_status(num_curr_cells)
    if num_curr_cells[1] >= 1000 && num_curr_cells[3] >= 6000
        return "happy"
    end

    return "sad"
end

# create and launch a Dash app
app = dash(update_title="")
app.title = "SuperHero Strawberry Demo"
time_min = time_span[1]
time_max = time_span[2]
num_time_slider_ticks = 6
slider_time_points = collect(LinRange(time_min, time_max, num_time_slider_ticks))
slider_markers = Dict([Int(i) => @sprintf("%d", i) for i in slider_time_points])
line_plot_time_points = collect(LinRange(time_min, time_max, 100))
num_d3_increments = 0
num_d3_decrements = 0
cell_colors = ["#4063D8", "#DAA520", "#CB3C33"]

# define page layout
app.layout = html_div() do
    html_h1("Welcome to stem cell superhero!"),
    html_p("Drag the slider to change current time. Can you make our strawberry superhero happy with the right balance of stem and blood cells?"),
    html_div(
        (
            dcc_graph(id="blood-cells"),
            dcc_graph(id="hero-status"),
            dcc_graph(id="bar-plot"),
            dcc_graph(id="line-plot"),
        ),
        id="output-panels",
    ),
    html_div(
        (
            html_div(
                (
                    # TODO: display actual parameter values
                    html_div(
                        (
                            html_span("Change stem cell death rate: "),
                            html_button(
                                id="increase-d1",
                                children="+",
                                n_clicks=0,
                            ),
                            html_button(
                                id="decrease-d1",
                                children="-",
                                n_clicks=0,
                            ),
                        ),
                        className="param-control-group",
                    ),
                    html_div(
                        (
                            html_span("Change red blood cell death rate: "),
                            html_button(
                                id="increase-d3",
                                children="+",
                                n_clicks=0,
                            ),
                            html_button(
                                id="decrease-d3",
                                children="-",
                                n_clicks=0,
                            ),
                        ),
                        className="param-control-group",
                    ),
                ),
                id="param-controls",
            ),
            dcc_slider(
                id="time-slider",
                included=true,
                min=time_min,
                max=time_max,
                step=0.1,
                marks=slider_markers,
                value=time_span[1],
                updatemode="drag",
            ),
        ),
        id="input-widgets",
    )
end

# define callbacks to user inputs
# update all panels given current time on time slider
callback!(
    app,
    Output("blood-cells", "figure"),
    Output("hero-status", "figure"),
    Output("bar-plot", "figure"),
    Output("line-plot", "figure"),
    Input("time-slider", "value"),
) do curr_time
    num_curr_cells = ss_solution(curr_time)
    num_cells = ss_solution(line_plot_time_points)
    # yaxis_range = [0, maximum(num_cells)]
    yaxis_range = [0, max(10000, maximum(num_cells))]

    # render blood cells according to ODE values at current time
    # TODO: add legend (color -> cell type)
    indiv_cell_colors = cat(
        [fill(cell_colors[i], trunc(Int, num_curr_cells[i]))
            for i=1:num_cell_types]...,
        dims=1)
    shuffle!(indiv_cell_colors)
    blood_cell_data = scatter(
        y=rand(length(indiv_cell_colors)),
        marker=attr(color=indiv_cell_colors),
        mode="markers",
    )
    blood_cell_layout = Layout(
        xaxis=attr(
            visible=false,
        ),
        yaxis=attr(
            visible=false,
        ),
    )

    # render hero status according to ODE values at current time
    hero_status = get_hero_status(num_curr_cells)
    hero_status_data = scatter(
        x=[0, 1],
        y=[0, 1],
        mode="markers",
        marker_opacity=0,
    )
    hero_status_layout = Layout(
        template=templates.plotly_white,
        xaxis=attr(
            visible=false,
            range=[0, 1],
        ),
        yaxis=attr(
            visible=false,
            range=[0, 1],
        ),
        images=[
            attr(
                source="assets/strawberry-superhero-$(hero_status).png",
                x=0.5,
                y=1.0,
                sizex=1.0,
                sizey=1.0,
                xanchor="center",
                yanchor="center",
            )
        ]
    )

    # render ODE values at on current time as bar plot
    # create bars for all cell types
    bar_plot_data = bar(
        x=cell_types,
        y=num_curr_cells,
        marker_color=cell_colors,
    )

    # define plot layout
    bar_plot_layout = Layout(
        xaxis_title="Cell type",
        yaxis_range=yaxis_range,
        yaxis_title="Number of cells",
        font=attr(size=16),
    )

    # render ODE values as line plot for the entire duration of time
    # create traces for all cell types
    line_plot_data = [
        scatter(
            x=line_plot_time_points,
            y=num_cells[i, :],
            name=cell_types[i],
            line=attr(color=cell_colors[i]),
        ) for i = 1:num_cell_types
    ]

    # create a vertical line marking current time
    curr_time_line = [vline(curr_time)]

    # define plot layout
    line_plot_layout = Layout(
        xaxis_title="Time",
        yaxis_range=yaxis_range,
        yaxis_title="Number of cells",
        shapes=curr_time_line,
        showlegend=true,
	linewidth=10,
        legend=attr(x=1, xanchor="right", y=1, bgcolor="#FFFFFF00"),
        font=attr(size=16),
    )

    return (
        Plot(blood_cell_data, blood_cell_layout),
        Plot(hero_status_data, hero_status_layout),
        Plot(bar_plot_data, bar_plot_layout),
        Plot(line_plot_data, line_plot_layout),
    )
end

# change parameter value from button input
callback!(
    app,
    Output("time-slider", "value"),
    Input("increase-d1", "n_clicks"),
    Input("decrease-d1", "n_clicks"),
    Input("increase-d3", "n_clicks"),
    Input("decrease-d3", "n_clicks"),
    State("time-slider", "value"),
) do _, _, _, _, curr_time
    ctx = callback_context()
    if length(ctx.triggered) == 0
        return curr_time
    end

    global ss_solution

    button_id = split(ctx.triggered[1].prop_id, ".")[1]
    action, param = split(button_id, "-")
    param_idx = ss_param2idx[param]
    param_old = ss_params[param_idx]

    if action == "increase"
        ss_params[param_idx] += 0.05
    else
        ss_params[param_idx] = max(ss_params_min[param_idx],
                                   ss_params[param_idx] - 0.05)
    end

    # update SDE solution only when the parameter has a new value
    if ss_params[param_idx] != param_old
        ss_solution = get_ss_solution()
    end

    curr_time
end


println("Launching Dash app...")
run_server(app, "0.0.0.0", 8050, debug=true)
using Printf
using Random

using DifferentialEquations

using Dash
using PlotlyJS

println("Initializing Dash app...")

# define a system to be simulated
# negative values in solution when parameter values are altered
function μ_func(du, u, p, t)
    SC, P, D = u
    λ, a₁, a₂, a₃, a₄, a₅, a₆, d₁, d₂, d₃ = p

    # stem
    du[1] = λ + a₁ * SC - a₃ * SC - d₁ * SC
    # progenitor
    du[2] = a₂ * SC + a₄ * P + 2 * a₃ * SC - a₆ * P - d₂ * P
    # differentiated
    du[3] = a₅ * P + 2a₆ * P - d₃ * D
end

function σ_func(du, u, p, t)
  du[1] = 10
  du[2] = 10
  du[3] = 10
end

# TODO: make parameter values and ODE solution local in user session
cell_types = ["Stem", "Progenitor", "Differentiated"]
num_cell_types = length(cell_types)
ss_params = [20.0, 0.33, 0.3, 0.35, 0.33, 0.33, 0.33, 0.0, 0.4, 0.5]
ss_params_min = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.24]
ss_param_names = ["lambda", "a1", "a2", "a3", "a4", "a5", "a6", "d1", "d2",
                  "d3"]
num_ss_params = length(ss_params)
ss_param2idx = Dict([ss_param_names[i] => i for i = 1:num_ss_params])
time_span = (0.0, 100.0)
init_condition = [1000.0, 100.0, 100.0]

# solve the system given current parameter values
function get_ss_solution()
    ss_sde_problem = SDEProblem(μ_func, σ_func, init_condition, time_span,
                                ss_params)
    solve(ss_sde_problem, SOSRI())
end

ss_solution = get_ss_solution()

# determine hero status based on number of cells
function get_hero_status(num_curr_cells)
    if num_curr_cells[1] >= 1000 && num_curr_cells[3] >= 6000
        return "happy"
    end

    return "sad"
end

# create and launch a Dash app
app = dash(update_title="")
app.title = "SuperHero Strawberry Demo"
time_min = time_span[1]
time_max = time_span[2]
num_time_slider_ticks = 6
slider_time_points = collect(LinRange(time_min, time_max, num_time_slider_ticks))
slider_markers = Dict([Int(i) => @sprintf("%d", i) for i in slider_time_points])
line_plot_time_points = collect(LinRange(time_min, time_max, 100))
num_d3_increments = 0
num_d3_decrements = 0
cell_colors = ["#4063D8", "#DAA520", "#CB3C33"]

# define page layout
app.layout = html_div() do
    html_h1("SuperHero Strawberry Demo"),
    html_p("Drag the slider at the bottom to change current time"),
    html_div(
        (
            dcc_graph(id="blood-cells"),
            dcc_graph(id="hero-status"),
            dcc_graph(id="bar-plot"),
            dcc_graph(id="line-plot"),
        ),
        id="output-panels",
    ),
    html_div(
        (
            html_div(
                (
                    # TODO: display actual parameter values
                    html_div(
                        (
                            html_span("Change death rate of stem cells: "),
                            html_button(
                                id="increase-d1",
                                children="+",
                                n_clicks=0,
                            ),
                            html_button(
                                id="decrease-d1",
                                children="-",
                                n_clicks=0,
                            ),
                        ),
                        className="param-control-group",
                    ),
                    html_div(
                        (
                            html_span("Change death rate of differentiated cells: "),
                            html_button(
                                id="increase-d3",
                                children="+",
                                n_clicks=0,
                            ),
                            html_button(
                                id="decrease-d3",
                                children="-",
                                n_clicks=0,
                            ),
                        ),
                        className="param-control-group",
                    ),
                ),
                id="param-controls",
            ),
            dcc_slider(
                id="time-slider",
                included=true,
                min=time_min,
                max=time_max,
                step=0.1,
                marks=slider_markers,
                value=time_span[1],
                updatemode="drag",
            ),
        ),
        id="input-widgets",
    )
end

# define callbacks to user inputs
# update all panels given current time on time slider
callback!(
    app,
    Output("blood-cells", "figure"),
    Output("hero-status", "figure"),
    Output("bar-plot", "figure"),
    Output("line-plot", "figure"),
    Input("time-slider", "value"),
) do curr_time
    num_curr_cells = ss_solution(curr_time)
    num_cells = ss_solution(line_plot_time_points)
    # yaxis_range = [0, maximum(num_cells)]
    yaxis_range = [0, max(10000, maximum(num_cells))]

    # render blood cells according to ODE values at current time
    # TODO: add legend (color -> cell type)
    indiv_cell_colors = cat(
        [fill(cell_colors[i], trunc(Int, num_curr_cells[i]))
            for i=1:num_cell_types]...,
        dims=1)
    shuffle!(indiv_cell_colors)
    blood_cell_data = scatter(
        y=rand(length(indiv_cell_colors)),
        marker=attr(color=indiv_cell_colors),
        mode="markers",
    )
    blood_cell_layout = Layout(
        xaxis=attr(
            visible=false,
        ),
        yaxis=attr(
            visible=false,
        ),
    )

    # render hero status according to ODE values at current time
    hero_status = get_hero_status(num_curr_cells)
    hero_status_data = scatter(
        x=[0, 1],
        y=[0, 1],
        mode="markers",
        marker_opacity=0,
    )
    hero_status_layout = Layout(
        template=templates.plotly_white,
        xaxis=attr(
            visible=false,
            range=[0, 1],
        ),
        yaxis=attr(
            visible=false,
            range=[0, 1],
        ),
        images=[
            attr(
                source="assets/strawberry-superhero-$(hero_status).png",
                x=0.5,
                y=1.0,
                sizex=1.0,
                sizey=1.0,
                xanchor="center",
                yanchor="center",
            )
        ]
    )

    # render ODE values at on current time as bar plot
    # create bars for all cell types
    bar_plot_data = bar(
        x=cell_types,
        y=num_curr_cells,
        marker_color=cell_colors,
    )

    # define plot layout
    bar_plot_layout = Layout(
        xaxis_title="Cell type",
        yaxis_range=yaxis_range,
        yaxis_title="Number of cells",
        font=attr(size=16),
    )

    # render ODE values as line plot for the entire duration of time
    # create traces for all cell types
    line_plot_data = [
        scatter(
            x=line_plot_time_points,
            y=num_cells[i, :],
            name=cell_types[i],
            line=attr(color=cell_colors[i]),
        ) for i = 1:num_cell_types
    ]

    # create a vertical line marking current time
    curr_time_line = [vline(curr_time)]

    # define plot layout
    line_plot_layout = Layout(
        xaxis_title="Time",
        yaxis_range=yaxis_range,
        yaxis_title="Number of cells",
        shapes=curr_time_line,
        showlegend=true,
        legend=attr(x=1, xanchor="right", y=1, bgcolor="#FFFFFF00"),
        font=attr(size=16),
    )

    return (
        Plot(blood_cell_data, blood_cell_layout),
        Plot(hero_status_data, hero_status_layout),
        Plot(bar_plot_data, bar_plot_layout),
        Plot(line_plot_data, line_plot_layout),
    )
end

# change parameter value from button input
callback!(
    app,
    Output("time-slider", "value"),
    Input("increase-d1", "n_clicks"),
    Input("decrease-d1", "n_clicks"),
    Input("increase-d3", "n_clicks"),
    Input("decrease-d3", "n_clicks"),
    State("time-slider", "value"),
) do _, _, _, _, curr_time
    ctx = callback_context()
    if length(ctx.triggered) == 0
        return curr_time
    end

    global ss_solution

    button_id = split(ctx.triggered[1].prop_id, ".")[1]
    action, param = split(button_id, "-")
    param_idx = ss_param2idx[param]
    param_old = ss_params[param_idx]

    if action == "increase"
        ss_params[param_idx] += 0.05
    else
        ss_params[param_idx] = max(ss_params_min[param_idx],
                                   ss_params[param_idx] - 0.05)
    end

    # update SDE solution only when the parameter has a new value
    if ss_params[param_idx] != param_old
        ss_solution = get_ss_solution()
    end

    curr_time
end


println("Launching Dash app...")
run_server(app, "0.0.0.0", 8050, debug=true)
