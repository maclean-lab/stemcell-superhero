using Printf
using Random

using ArgParse
using JSON3

using DifferentialEquations
using Dash
using PlotlyJS

"Parse command line arguments."
function get_command_line_args()
    arg_parse_settings = ArgParseSettings()

    @add_arg_table arg_parse_settings begin
        "--path_prefix"
            help = "path prefix for this web app, e.g. \"/scsh/\" when deployed
                    on macleanlab.usc.edu/scsh"
            arg_type = String
            default = "/"
        "--debug"
            help = "enable Dash dev tools"
            action = :store_true
    end

    return parse_args(arg_parse_settings)
end

# define a system to be simulated
function μ_func!(du, u, p, t)
    SC, P, D = u
    λ, a₁, a₂, a₃, a₄, a₅, a₆, d₁, d₂, d₃ = p

    # stem
    du[1] = λ + a₁ * SC - a₃ * SC - d₁ * SC
    # progenitor
    du[2] = a₂ * SC + a₄ * P + 2 * a₃ * SC - a₆ * P - d₂ * P
    # differentiated
    du[3] = a₅ * P + 2a₆ * P - d₃ * D
end

function σ_func!(du, u, p, t)
    du[1] = 10
    du[2] = 10
    du[3] = 10
end

cell_types = ["Stem cells", "Progenitor cells", "Red blood cells"]
cell_colors = ["#4063D8", "#DAA520", "#CB3C33"]
num_cell_types = length(cell_types)
ss_params = ["lambda", "a1", "a2", "a3", "a4", "a5", "a6", "d1", "d2", "d3"]
ss_param_vals_init = [25.0, 0.33, 0.3, 0.35, 0.33, 0.33, 0.33, 0.01, 0.4, 0.5]
ss_param_mins = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.24]
ss_param_indices = Dict(pn => i for (i, pn) in enumerate(ss_params))
ss_param_increments = Dict("d1" => 0.01, "d3" => 0.05)
time_min = 0.0
time_max = 100.0
time_span = (time_min, time_max)
time_step = 0.1
time_points = collect(time_min:time_step:time_max)
time_point_indices = Dict(t => i for (i, t) in enumerate(time_points))
num_cells_init = [1000.0, 100.0, 100.0]

"Solve the system given current parameter values."
function get_ss_solution(params)
    # solve the system with given parameters
    ss_prob = SDEProblem(μ_func!, σ_func!, num_cells_init, time_span, params;
                        seed=0)
    ss_sol = Array(solve(ss_prob, SOSRI(), saveat=time_step))

    # return solution as a dictionary
    return Dict(cell_types[i] => ss_sol[i, :] for i in 1:num_cell_types)
end

"Determine hero status based on number of cells."
function get_hero_status(num_cells, curr_time)
    if curr_time > 10
        if num_cells[3] >= 2000 || num_cells[3] <= 1000
            return "sad"
        elseif 0.01 <= num_cells[1] / num_cells[3] <= 0.15
            return "happy"
        else
            return "neutral"
        end
    end

    return "neutral"
end




"Get prompt for displaying parameter value in parameter control widget."
function get_param_control_prompt(param_full_name, param_value)
    return @sprintf("%s: %.2f", param_full_name, param_value)
end

"Create widgets for playing with parameters and time."
function get_input_widgets()
    param_control_widgets = html_div(
        children=[
            html_div(
                children=[
                    html_span(
                        children=get_param_control_prompt(
                            "Stem cell death rate",
                            ss_param_vals_init[ss_param_indices["d1"]]
                        ),
                        id="d1-value",
                    ),
                    html_button(
                        children="-",
                        id="decrease-d1",
                        n_clicks=0,
                    ),
                    html_button(
                        children="+",
                        id="increase-d1",
                        n_clicks=0,
                    ),
                ],
                className="param-control-widget",
            ),
            html_div(
                children=[
                    html_span(
                        children=get_param_control_prompt(
                            "Red blood cell death rate",
                            ss_param_vals_init[ss_param_indices["d3"]]
                        ),
                        id="d3-value",
                    ),
                    html_button(
                        id="decrease-d3",
                        children="-",
                        n_clicks=0,
                    ),
                    html_button(
                        id="increase-d3",
                        children="+",
                        n_clicks=0,
                    ),
                ],
                className="param-control-widget",
            ),
        ],
        id="param-controls",
    )

    # create time slider widget
    time_slider_step_size = 20
    time_slider_ticks = time_min:time_slider_step_size:time_max
    time_slider_markers = Dict(
        Int(i) => @sprintf("%d", i) for i in time_slider_ticks
    )

    time_slider = dcc_slider(
        id="time-slider", included=true, min=time_min, max=time_max,
        step=time_step, marks=time_slider_markers, value=time_min,
        updatemode="drag",
    )

    return param_control_widgets, time_slider
end

"Add a callback that updates all output panels."
function add_output_update_callback(app, assets_dir)
    callback!(
        app,
        Output("blood-cells", "figure"),
        Output("hero-status", "figure"),
        Output("bar-plot", "figure"),
        Output("line-plot", "figure"),
        Input("time-slider", "value"),
        Input("ss-solution", "data")
    ) do curr_time, ss_solution
        curr_time_idx = time_point_indices[curr_time]
        num_curr_cells = [ss_solution[ct][curr_time_idx] for ct in cell_types]
        yaxis_max = max(
            5000,
            maximum(maximum(cell_nums) for (_, cell_nums) in ss_solution)
        )
        yaxis_range = [0, yaxis_max]

        # render blood cells according to cell numbers at current time
        # for performance we only render 1/10 cells
        indiv_cell_colors = cat(
            [fill(cell_colors[i], trunc(Int, num_curr_cells[i] / 10))
            for i=1:num_cell_types]...,
            dims=1)
        num_bone_cells = length(indiv_cell_colors)
        shuffle!(Random.seed!(0), indiv_cell_colors)
        bone_x0 = -1.0
        bone_y0 = 0.05
        bone_x1 = Float64(num_bone_cells) + 1.0
        bone_y1 = 0.65
        blood_cell_data = scatter(
            y=rand(num_bone_cells) .* 0.5 .+ 0.1,
            marker=attr(color=indiv_cell_colors),
            mode="markers",
        )
        blood_cell_layout = Layout(
            template=templates.plotly_white,
            xaxis=attr(
                visible=false,
            ),
            yaxis=attr(
                visible=false,
                range=[0.0, 1.0],
            ),
            shapes=[
                line(x0=bone_x0, y0=bone_y0, x1=bone_x1, y1=bone_y0;
                     line_color="#B8A487", line_width=8),
                line(x0=bone_x0, y0=bone_y1, x1=bone_x1, y1=bone_y1;
                     line_color="#B8A487", line_width=8),
                rect(x0=bone_x0, y0=bone_y0, x1=bone_x1, y1=bone_y1;
                     fillcolor="#E6D5B3", layer="below", line_width=0),
            ],
            images=[
                attr(
                    source=joinpath(assets_dir, "bone_marrow.png"),
                    x=0.5,
                    y=1.0,
                    sizex=1.0,
                    sizey=1.0,
                    xanchor="center",
                    yanchor="center",
                    layer="below",
                )
            ],
        )

        # render hero status according to cell numbers at current time
        hero_status = get_hero_status(num_curr_cells, curr_time)
        hero_status_data = scatter(
            x=[0.12, 0.12, 0.12],
            y=[0.65, 0.5, 0.35],
            mode="markers",
            marker=attr(color=cell_colors, size=20),
        )
        hero_status_layout = Layout(
            template=templates.plotly_white,
            xaxis=attr(
                visible=false,
                range=[0.0, 1.0],
            ),
            yaxis=attr(
                visible=false,
                range=[0.0, 1.0],
            ),
            images=[
                attr(
                    source=joinpath(assets_dir, "hero_$(hero_status).png"),
                    x=0.65,
                    y=1.0,
                    sizex=1.0,
                    sizey=1.0,
                    xanchor="center",
                    yanchor="center",
                    layer="below",
                )
            ],
            annotations=[
                attr(
                    x=0.15,
                    y=0.65,
                    xanchor="left",
                    text="Stem cells",
                    font=attr(size=16, color=:black),
                    showarrow=false,
                ),
                attr(
                    x=0.15,
                    y=0.5,
                    xanchor="left",
                    text="Progenitor cells",
                    font=attr(size=16, color=:black),
                    showarrow=false,
                ),
                attr(
                    x=0.15,
                    y=0.35,
                    xanchor="left",
                    text="Red blood cells",
                    font=attr(size=16, color=:black),
                    showarrow=false,
                ),
            ],
        )

        # render cell numbers at on current time as bar plot
        bar_plot_data = bar(
            x=cell_types,
            y=num_curr_cells,
            marker_color=cell_colors,
        )
        bar_plot_layout = Layout(
            xaxis_title="Cell type",
            yaxis_range=yaxis_range,
            yaxis_title="Number of cells",
            font=attr(size=16),
        )

        # render cell numbers as line plot for the entire duration of time
        line_plot_data = [
            scatter(
                x=time_points,
                y=ss_solution[cell_types[i]],
                name=cell_types[i],
                line=attr(color=cell_colors[i]),
            ) for i = 1:num_cell_types
        ]
        curr_time_line = [vline(curr_time)]  # vertical line at current time
        line_plot_layout = Layout(
            xaxis=attr(dtick=20),
            xaxis_title="Time",
            yaxis_range=yaxis_range,
            yaxis_title="Number of cells",
            shapes=curr_time_line,
            showlegend=false,
            linewidth=10,
            font=attr(size=16),
        )

        return (
            Plot(blood_cell_data, blood_cell_layout),
            Plot(hero_status_data, hero_status_layout),
            Plot(bar_plot_data, bar_plot_layout),
            Plot(line_plot_data, line_plot_layout),
        )
    end
end

"Add a callback that changes parameter value from button input."
function add_param_control_callback(app)
    callback!(
        app,
        Output("ss-params", "data"),
        Output("ss-solution", "data"),
        Output("d1-value", "children"),
        Output("d3-value", "children"),
        Input("increase-d1", "n_clicks"),
        Input("decrease-d1", "n_clicks"),
        Input("increase-d3", "n_clicks"),
        Input("decrease-d3", "n_clicks"),
        State("ss-params", "data"),
        prevent_initial_call=true
    ) do _, _, _, _, param_data
        param_vals = Float64.(copy(param_data))
        ctx = callback_context()
        button_id = split(ctx.triggered[1].prop_id, ".")[1]
        action, param = split(button_id, "-")
        param_idx = ss_param_indices[param]
        param_old = param_vals[param_idx]

        if action == "increase"
            param_vals[param_idx] += ss_param_increments[param]
        else
            param_vals[param_idx] = max(
                ss_param_mins[param_idx],
                param_vals[param_idx] - ss_param_increments[param]
            )
        end

        if param_vals[param_idx] == param_old
            throw(PreventUpdate())
        end

        return (
            param_vals,
            get_ss_solution(param_vals),
            get_param_control_prompt(
                "Stem cell death rate",
                param_vals[ss_param_indices["d1"]]
            ),
            get_param_control_prompt(
                "Red blood cell death rate",
                param_vals[ss_param_indices["d3"]]
            ),
        )
    end
end

"Create and launch the Stem Cell Superhero app."
function main()
    println("Initializing Dash app...")
    flush(stdout)

    # get command line arguments
    command_line_args = get_command_line_args()
    path_prefix = command_line_args["path_prefix"]
    assets_dir = "$(path_prefix)assets"

    # create and launch a Dash app
    app = dash(update_title="", requests_pathname_prefix=path_prefix)
    app.title = "Stem Cell Superhero"
    plot_config = Dict("staticPlot"=>true)

    # define page layout
    app.layout = html_div() do
        dcc_store(id="ss-params", data=ss_param_vals_init),
        dcc_store(id="ss-solution", data=get_ss_solution(ss_param_vals_init)),
        html_h1(children="Welcome to stem cell superhero!"),
        html_p(children="Drag the slider to change current time. Can you make
                         our strawberry superhero happy with the right balance
                         of stem and blood cells?"),
        html_div(
            children=[
                dcc_graph(id="blood-cells", config=plot_config),
                dcc_graph(id="hero-status", config=plot_config),
                dcc_graph(id="bar-plot", config=plot_config),
                dcc_graph(id="line-plot", config=plot_config),
            ],
            id="output-panels",
        ),
        html_div(children=get_input_widgets(), id="input-widgets")
    end

    # add callbacks
    add_output_update_callback(app, assets_dir)
    add_param_control_callback(app)

    # launch Dash app
    println("Launching Dash app...")
    flush(stdout)
    run_server(app, "0.0.0.0", 8080, debug=command_line_args["debug"])
end

main()