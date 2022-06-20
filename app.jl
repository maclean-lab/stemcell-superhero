using Dash
using PlotlyJS

using Catalyst
using DifferentialEquations

# define a system to be simulated
println("Initializing Dash app...")
ss_system_spec = @reaction_network begin
    λ, 0 --> SC # production of SC
    a₁, SC --> 2 * SC
    a₂, SC --> SC + P
    a₃, SC --> 2 * P
    a₄, P --> 2 * P
    a₅, P --> P + D
    a₆, P --> 2 * D
    d₁, SC --> 0 # death rate of SC
    d₂, P --> 0 # death rate of P
    d₃, D --> 0 # death rate of D
end λ a₁ a₂ a₃ a₄ a₅ a₆ d₁ d₂ d₃

species_names = ["Stem", "Progenitor", "Differentiated"]
num_species = size(species_names, 1)
ss_params = [0.25, 0.3, 0.3, 0.4, 0.33, 0.33, 0.33, 0.0, 0.1, 0.2]
time_span = (0.0, 200.0)
init_condition = [100.0, 0.0, 0.0] # [SC,P,D]

ss_ode_problem = ODEProblem(ss_system_spec, init_condition, time_span,
                            ss_params)
ss_solution = solve(ss_ode_problem, Tsit5())

# create and launch a Dash app
app = dash(update_title = "")
app.title = "SuperHero Strawberry Demo"
time_min = time_span[1]
time_max = time_span[2] / 2
num_time_slider_ticks = 5
slider_time_points = collect(LinRange(time_min, time_max, num_time_slider_ticks))
slider_markers = Dict([Int(i) => string(i) for i in slider_time_points])
line_plot_time_points = collect(LinRange(time_min, time_max, 100))
num_d3_increments = 0
num_d3_decrements = 0

# define page layout
app.layout = html_div() do
    html_h1("SuperHero Strawberry Demo"),
    html_p("Drag the slider at the bottom to change current time"),
    html_div(
        (
            dcc_graph(id = "blood-cells"),
            dcc_graph(id = "hero-status"),
            dcc_graph(id = "bar-plot"),
            dcc_graph(id = "line-plot"),
        ),
        id = "output-panels",
    ),
    html_div(
        (
            dcc_slider(
                id = "time-slider",
                included = true,
                min = time_min,
                max = time_max,
                step = 0.1,
                marks = slider_markers,
                value = time_span[1],
                updatemode = "drag",
            ),
            html_span("Change death rate of differentiated cells: "),
            html_button(
                id = "increase-d3",
                children = "+",
                n_clicks = 0,
            ),
            html_button(
                id = "decrease-d3",
                children = "-",
                n_clicks = 0,
            ),
        ),
        id = "input-widgets",
    )
end

# define callbacks to user inputs
callback!(
    app,
    Output("blood-cells", "figure"),
    Output("hero-status", "figure"),
    Output("bar-plot", "figure"),
    Output("line-plot", "figure"),
    Input("time-slider", "value"),
    Input("increase-d3", "n_clicks"),
    Input("decrease-d3", "n_clicks"),
) do curr_time, num_increase_clicks, num_decrease_clicks
    global num_d3_increments
    global num_d3_decrements
    global ss_solution

    if num_increase_clicks != num_d3_increments || num_decrease_clicks != num_d3_decrements
        if num_increase_clicks != num_d3_increments
            num_d3_increments = num_increase_clicks
            ss_params[end] += 0.05
        else
            num_d3_decrements = num_decrease_clicks
            ss_params[end] -= 0.05
        end

        ss_ode_problem = ODEProblem(ss_system_spec, init_condition, time_span,
                                    ss_params)
        ss_solution = solve(ss_ode_problem, Tsit5())
    end
    ss_solution_data = ss_solution(line_plot_time_points)
    yaxis_range = [0, maximum(ss_solution_data)]

    # TODO: render blood cells according to ODE values at current time
    blood_cell_data = scatter(
        x = [0],
        y = [0],
        text = ["Placeholder for blood cells"],
        textfont_size = 24,
        mode = "text",
    )
    blood_cell_layout = Layout(
        title = string("Current time: ", curr_time),
    )

    # TODO: render hero status according to ODE values at current time
    hero_status_data = scatter(
        x = [0],
        y = [0],
        text = ["Placeholder for hero status"],
        textfont_size = 24,
        mode = "text",
    )

    # render ODE values at on current time as bar plot
    # create bars for all species
    bar_plot_data = bar(
        x = species_names,
        y = ss_solution(curr_time),
        marker_color = 1:num_species,
    )

    # define plot layout
    bar_plot_layout = Layout(
        xaxis_title = "Species",
        yaxis_range = yaxis_range,
        yaxis_title = "Value",
    )

    # render ODE values as line plot for the entire duration of time
    # create traces for all species
    line_plot_data = [
        scatter(
            x = line_plot_time_points,
            y = ss_solution_data[i, :],
            name = species_names[i],
        ) for i = 1:num_species
    ]

    # create a vertical line marking current time
    curr_time_line = [vline(curr_time)]

    # define plot layout
    line_plot_layout = Layout(
        xaxis_title = "Time",
        yaxis_range = yaxis_range,
        yaxis_title = "Value",
        shapes = curr_time_line,
        showlegend = true,
        legend = attr(x = 1, xanchor = "right", y = 1, bgcolor = "#FFFFFF00"),
    )

    return (
        Plot(blood_cell_data, blood_cell_layout),
        Plot(hero_status_data),
        Plot(bar_plot_data, bar_plot_layout),
        Plot(line_plot_data, line_plot_layout),
    )
end

println("Lauching Dash app...")
run_server(app, "0.0.0.0", debug = true)
