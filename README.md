# MacLean Lab Outreach Project
Interactive cell differentiation model demo

## Set up the web app for the first time
1. Install Julia 1.7
2. Clone this repository
3. Install all packages specified in `Project.toml`
    - In a terminal, `cd` into the cloned repository
    - Launch Julia REPL using the following command
        ```
        julia --project
        ```
    - Press `]` to switch to Julia's package management mode
    - Enter `instantiate` to install on pacakges
    - Exit Julia REPL once all packages are installed

## Run the web app
### Run on the server
```
nohup julia --project app.jl --path_prefix="/scsh/" > app.log &
```
The app will run in the background and keep running after user logs out.
Any output is saved to `app.log`.

### Run locally for development
```
julia --project app.jl --debug
```
Dash dev tools are enabled with the `--debug` flag.
