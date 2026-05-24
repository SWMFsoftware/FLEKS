using Batsrus
using PyPlot

println("====================================================")
println("SWMF/FLEKS 2D Fast Wave Propagation Plotting Tool")
println("====================================================")

# Paths to the test case output directories
const PIC_DIR = joinpath(@__DIR__, "fast_wave_pic", "output")
const HYBRID_DIR = joinpath(@__DIR__, "fast_wave_hybrid", "output")

# Function to extract cycle number from filename
function get_cycle_from_filename(filename::String)
    # Filenames are like: z=0_fluid_region1_1_t00000000_n00000000.out
    m = match(r"_n(\d+)\.out$", filename)
    if m === nothing
        return -1
    else
        return parse(Int, m.captures[1])
    end
end

# Function to load a 1D slice from a 2D dataset
function get_1d_slice(file_path::String, var_name::String)
    data = load(file_path)
    
    # Robustly handle access (property vs dictionary)
    x_raw = hasproperty(data, :x) ? data.x : data["x"]
    
    sym = Symbol(var_name)
    var_raw = hasproperty(data, sym) ? getproperty(data, sym) : data[var_name]
    
    # Convert to standard Julia Array to bypass DimensionalData's custom indexing quirks/bugs
    x = collect(x_raw)
    var = collect(var_raw)
    
    # SWMF grids:
    # x has coordinates (nx, ny, ndim) in 2D, or (nx, ny, nz, ndim) in 3D.
    # var has shape (nx, ny) in 2D, or (nx, ny, nz) in 3D.
    var_ndims = length(size(var))
    
    if var_ndims == 1
        # 1D dataset
        # size(x) is (nx, 1) where x[:, 1] is x coordinates
        return x[:, 1], var
    elseif var_ndims == 2
        # 2D dataset
        # size(x) is (nx, ny, 2) where x[:, :, 1] is x coordinates
        nx, ny = size(var)
        mid_y = div(ny, 2) + 1
        return x[:, mid_y, 1], var[:, mid_y]
    elseif var_ndims == 3
        # 3D dataset
        # size(x) is (nx, ny, nz, 3) where x[:, :, :, 1] is x coordinates
        nx, ny, nz = size(var)
        mid_y = div(ny, 2) + 1
        mid_z = div(nz, 2) + 1
        return x[:, mid_y, mid_z, 1], var[:, mid_y, mid_z]
    else
        error("Unsupported variable dimensions: $var_ndims")
    end
end

function main()
    # Check if directories exist
    if !isdir(PIC_DIR) || !isdir(HYBRID_DIR)
        println("ERROR: Output directories not found!")
        println("Please run 'python3 tests/test_standalone/validate_tests.py' first to generate and save output files.")
        return
    end

    # Get and sort out files
    pic_files = filter(f -> endswith(f, ".out"), readdir(PIC_DIR))
    hybrid_files = filter(f -> endswith(f, ".out"), readdir(HYBRID_DIR))

    # Parse and sort by cycle number
    pic_cycles = Dict(get_cycle_from_filename(f) => joinpath(PIC_DIR, f) for f in pic_files if get_cycle_from_filename(f) >= 0)
    hybrid_cycles = Dict(get_cycle_from_filename(f) => joinpath(HYBRID_DIR, f) for f in hybrid_files if get_cycle_from_filename(f) >= 0)

    # Find common cycles
    common_cycles = sort(collect(intersect(keys(pic_cycles), keys(hybrid_cycles))))

    if isempty(common_cycles)
        println("ERROR: No common output cycles found between PIC and Hybrid PIC simulations!")
        return
    end

    println("Found $(length(common_cycles)) common output cycles: ", common_cycles)

    # We want to select a few key cycles to plot.
    # E.g. t = 0.0 (cycle 0), t = 320.0 (cycle 1600), t = 640.0 (cycle 3200)
    # Let's plot 3 time steps: Initial (t=0.0), Middle (t=320.0), and Final (t=640.0).
    selected_cycles = Int[]
    if 0 in common_cycles
        push!(selected_cycles, 0)
    end
    
    # Try to find a cycle near the middle (320.0s / 1600 cycles)
    mid_target = 1600
    if !isempty(common_cycles)
        mid_cycle = common_cycles[argmin(abs.(common_cycles .- mid_target))]
        if !(mid_cycle in selected_cycles)
            push!(selected_cycles, mid_cycle)
        end
    end
    
    # Try to find a cycle near the end (640.0s / 3200 cycles)
    end_target = 3200
    if !isempty(common_cycles)
        end_cycle = common_cycles[argmin(abs.(common_cycles .- end_target))]
        if !(end_cycle in selected_cycles)
            push!(selected_cycles, end_cycle)
        end
    end

    # Sort selected cycles
    selected_cycles = sort(selected_cycles)
    println("Selected cycles for plotting: ", selected_cycles)

    # Setup the multi-panel figure
    # 3 rows (density rho0, velocity ux0, magnetic field by), columns for each selected cycle/time
    num_times = length(selected_cycles)
    fig, axes = plt.subplots(3, num_times, figsize=(5 * num_times, 9), sharex="row", sharey="row")
    
    # Ensure axes is always a 2D array even if num_times == 1
    if num_times == 1
        axes = reshape(axes, (3, 1))
    end

    # Plot parameters
    dt = 0.2
    
    for (col_idx, cycle) in enumerate(selected_cycles)
        t_phys = cycle * dt
        
        # Load files
        pic_file = pic_cycles[cycle]
        hybrid_file = hybrid_cycles[cycle]
        
        # Plot Density (Row 1)
        ax_rho = axes[1, col_idx]
        x_pic, rho_pic = get_1d_slice(pic_file, "rhos0")
        x_hyb, rho_hyb = get_1d_slice(hybrid_file, "rhos0")
        
        ax_rho.plot(x_pic, rho_pic, label="Full PIC", color="#4361ee", linewidth=2.0)
        ax_rho.plot(x_hyb, rho_hyb, label="Hybrid PIC", color="#f72585", linewidth=2.0, linestyle="--")
        ax_rho.set_title("t = $(t_phys) s (Cycle $cycle)")
        if col_idx == 1
            ax_rho.set_ylabel("Ion Density (rhos0)", fontsize=11, fontweight="bold")
        end
        ax_rho.grid(true, linestyle=":", alpha=0.6)
        
        # Plot Velocity (Row 2)
        ax_ux = axes[2, col_idx]
        _, ux_pic = get_1d_slice(pic_file, "uxs0")
        _, ux_hyb = get_1d_slice(hybrid_file, "uxs0")
        
        ax_ux.plot(x_pic, ux_pic, color="#4361ee", linewidth=2.0)
        ax_ux.plot(x_hyb, ux_hyb, color="#f72585", linewidth=2.0, linestyle="--")
        if col_idx == 1
            ax_ux.set_ylabel("Ion Velocity (uxs0)", fontsize=11, fontweight="bold")
        end
        ax_ux.grid(true, linestyle=":", alpha=0.6)

        # Plot Magnetic Field (Row 3)
        ax_by = axes[3, col_idx]
        _, by_pic = get_1d_slice(pic_file, "by")
        _, by_hyb = get_1d_slice(hybrid_file, "by")
        
        ax_by.plot(x_pic, by_pic, color="#4361ee", linewidth=2.0)
        ax_by.plot(x_hyb, by_hyb, color="#f72585", linewidth=2.0, linestyle="--")
        if col_idx == 1
            ax_by.set_ylabel("Magnetic Field (by)", fontsize=11, fontweight="bold")
        end
        ax_by.set_xlabel("x", fontsize=11, fontweight="bold")
        ax_by.grid(true, linestyle=":", alpha=0.6)
    end

    # Add legend to the first panel
    axes[1, 1].legend(loc="upper right", frameon=true, shadow=false, fancybox=true)

    # Tight layout and save
    plt.tight_layout()
    output_png = joinpath(@__DIR__, "wave_propagation_comparison.png")
    plt.savefig(output_png, dpi=300)
    plt.close()
    
    println("====================================================")
    println("SUCCESS: Generated comparison plot at:")
    println(output_png)
    println("====================================================")
end

main()
