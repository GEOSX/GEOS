#!/opt/cray/pe/python/3.11.5/bin/python

import matplotlib.pyplot as plt
import argparse
import glob
import sys
import re
import os

def parse_log_file(file_path):
    """Parse the log file from GEOS and extract required statistics."""
    try:
        print(f"Parsing {file_path = }")
        with open(file_path, 'r') as file:
            data = file.read()

    except FileNotFoundError:
        print(f"Error: File not found - {file_path}", file=sys.stderr)

    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}", file=sys.stderr)

    # Regex patterns to extract the required information
    patterns = {
        'non_linear_iters': re.compile(r'number of successful nonlinear iterations: (\d+).*number of discarded nonlinear iterations: (\d+)', re.S),
        'linear_iters': re.compile(r'number of successful linear iterations: (\d+).*number of discarded linear iterations: (\d+)', re.S),
        'applysol_time': re.compile(r'apply solution time = (\d+\.\d+) s'),
        'assemble_time': re.compile(r'assemble time = (\d+\.\d+) s'),
        'convcheck_time': re.compile(r'convergence check time = (\d+\.\d+) s'),
        'update_state_time': re.compile(r'update state time = (\d+\.\d+) s'),
        'ls_create_time': re.compile(r'linear solver create time = (\d+.\d+) s'),
        'ls_setup_time': re.compile(r'linear solver setup time = (\d+.\d+) s'),
        'ls_solve_time': re.compile(r'linear solver solve time = (\d+.\d+) s'),
        'run_time': re.compile(r'run time\s+.*\((\d+\.\d+) s\)')
    }

    # Dictionary to hold the results
    stats = {}

    # Extract and process all patterns
    for key, pattern in patterns.items():
        if key in ['non_linear_iters', 'linear_iters']:
            xmatch = pattern.search(data)
            stats[key] = int(xmatch.group(1)) + int(xmatch.group(2)) if xmatch else 0

        else:
            xmatch = pattern.search(data)
            stats[key] = float(xmatch.group(1)) if xmatch else 0

    # Calculate GEOS time (run time excluding hypre time)
    stats['geos_time'] = stats['run_time'] - (stats['ls_create_time'] + stats['ls_setup_time'] + stats['ls_solve_time'])

    # Calculate averages
    averages = {
        'avg_geos_time': 'geos_time',
        'avg_ls_create_time': 'ls_create_time',
        'avg_ls_setup_time': 'ls_setup_time',
        'avg_ls_solve_time': 'ls_solve_time',
        'avg_ls_iters': 'linear_iters',
        'avg_run_time': 'run_time'
    }

    if stats['non_linear_iters'] <= 0:
        print(f"Invalid number of non-linear iterations!", file=sys.stderr)
        return None

    for avg_key, total_key in averages.items():
        stats[avg_key] = stats[total_key] / stats['non_linear_iters']

    # Store level
    xmatch = re.search(r'level(\d+)', file_path)
    stats['level'] = int(xmatch.group(1)) if xmatch else -1

    return stats

def print_results(model, results):
    """Print results in markdown table format"""
    if not results or not all(results):
        print(f"Invalid results!", file=sys.stderr)
        return

    # Setup header
    headers = ["Level", "GEOS (s)", "Matrix creation (s)", "Hypre setup (s)", "Hypre solve (s)", "Total (s)", "Hypre iters"]

    # Gather all rows data including headers
    data_rows = [
        headers,
        *[
            [
                "{:02d}".format(result['level']),
                "{:.4f}".format(result['avg_geos_time']),
                "{:.4f}".format(result['avg_ls_create_time']),
                "{:.4f}".format(result['avg_ls_setup_time']),
                "{:.4f}".format(result['avg_ls_solve_time']),
                "{:.4f}".format(result['avg_run_time']),
                "{:.1f}".format(result['avg_ls_iters']),
            ]
            for result in results
        ]
    ]

    # Calculate maximum column widths
    max_widths = [max(len(row[i]) for row in data_rows) for i in range(len(headers))]

    # Print header row with adjusted widths
    header_row = "| " + " | ".join(header.ljust(max_widths[i]) for i, header in enumerate(headers)) + " |"
    print(f"\nSummary table for {model}:\n")
    print(header_row)

    # Adjust the divider line based on maximum widths
    divider_line = "|" + "|".join(["-" * (max_widths[i] + 2) for i in range(len(headers))]) + "|"
    print(divider_line)

    # Print each data row with adjusted widths
    for row in data_rows[1:]:  # Skip the header in data_rows for output
        formatted_row = "| " + " | ".join(item.ljust(max_widths[i]) for i, item in enumerate(row)) + " |"
        print(formatted_row)

def plot_results(model, model_data, results, fs=16, ms=10):
    if not results:
        return

    plt.rcParams.update({'font.size': fs})
    fig, ax = plt.subplots(figsize=(10, 6))

    # Extracting ranks and global DOFs
    levels = [result['level'] for result in results]
    ranks = [entry['ranks'] for entry in model_data if entry['level'] in levels]
    global_dofs = [entry['global_dofs'] for entry in model_data if entry['level'] in levels]

    # Setting up data for plotting
    geos_times = [result['avg_geos_time'] for result in results]
    ls_create_times = [result['avg_ls_create_time'] for result in results]
    ls_setup_times = [result['avg_ls_setup_time'] for result in results]
    ls_solve_times = [result['avg_ls_solve_time'] for result in results]
    total_times = [result['avg_run_time'] for result in results]

    # Plotting each metric
    ax.plot(ranks, geos_times, 'o-', markersize=ms, label='GEOS')
    ax.plot(ranks, ls_create_times, 's-', markersize=ms, label='Matrix Creation')
    ax.plot(ranks, ls_setup_times, '^-', markersize=ms, label='Hypre Setup')
    ax.plot(ranks, ls_solve_times, 'v-', markersize=ms, label='Hypre Solve')
    ax.plot(ranks, total_times, 'd-', markersize=ms, label='Total')

    # Adding labels and legend
    ax.set_title(f"Weak scaling - {model}", fontsize=fs+4)
    ax.set_xlabel('Number of GPUs (Global DOFs)', labelpad=18, fontsize=fs+1)
    ax.set_ylabel('Time (s)', labelpad=18, fontsize=fs+1)
    ax.set_xscale('log', base=2)
    ax.set_xticks(ranks)
    ax.set_xticklabels([f"{r:,} ({d})" for r, d in zip(ranks, global_dofs)])
    ax.legend(fontsize=fs)

    # Show grid lines
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')  # Add gridlines with custom settings

    # Update margins
    plt.tight_layout()

    # Save figure to file
    figname = f"weakscal_{model}.png"
    print(f"Saving figure to file {figname}...")
    plt.savefig(figname, format='png', dpi=450)

    # Show the plot
    plt.show()

def find_newest_file(model, level):
    """Finds the newest output file under the directory relative to a physics model and for certain problem level (mesh refinement)."""
    target_dir = os.path.join(model, f"level{level:02d}")
    files = glob.glob(os.path.join(target_dir, '*.out'))
    if not files:
        return target_dir, None

    # Extract the job ID numbers and sort the files by these IDs
    suffixes = {'mechanics': 'mechanics',
                'singlePhaseFlow': 'singlePhaseFlow',
                'compositionalMultiphaseFlow': 'compflow'}
    files_with_ids = [(file, int(re.search(rf"(\d+)-{suffixes[model]}", file).group(1))) for file in files]

    # Sort by the job ID number
    return target_dir, max(files_with_ids, key=lambda x: x[1])[0]

def main():
    """Process multiple input files and store results using command line arguments."""
    valid_models = ['mechanics', 'singlePhaseFlow', 'compositionalMultiphaseFlow']

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', required=True, choices=valid_models, help="Physics model")
    parser.add_argument('-l', '--levels', nargs='+', type=int, help="Levels to process")
    args = parser.parse_args()

    # Setup problem info
    problems = {
        "mechanics": [
            {"level": 2, "ranks": 8, "global_dofs": "20.7M"},
            {"level": 3, "ranks": 64, "global_dofs": "162M"},
            {"level": 4, "ranks": 512, "global_dofs": "1.3B"},
            {"level": 5, "ranks": 4096, "global_dofs": "10.2B"},
            {"level": 6, "ranks": 32768, "global_dofs": "81.3B"}
        ],
        "singlePhaseFlow": [
            {"level": 2, "ranks": 4, "global_dofs": "6.6M"},
            {"level": 3, "ranks": 32, "global_dofs": "52.8M"},
            {"level": 4, "ranks": 256, "global_dofs": "422M"},
            {"level": 5, "ranks": 2048, "global_dofs": "3.4B"},
            {"level": 6, "ranks": 16384, "global_dofs": "27B"}
        ],
        "compositionalMultiphaseFlow": [
            {"level": 2, "ranks": 4, "global_dofs": "19.8M"},
            {"level": 3, "ranks": 32, "global_dofs": "153M"},
            {"level": 4, "ranks": 256, "global_dofs": "1.2B"},
            {"level": 5, "ranks": 2048, "global_dofs": "9.8B"}
        ]
    }

    # Set default levels if not passed in
    if args.levels is None:
        args.levels = [entry['level'] for entry in problems[args.model]]

    # Read log files
    results = []
    file_paths = []
    for level in args.levels:
        target_dir, file_path = find_newest_file(args.model, level)
        if file_path:
            file_paths.append(file_path)
        else:
            print(f"No suitable files found in {target_dir}")

    # Parse statistics
    results = [parse_log_file(fp) for fp in file_paths]

    # Print summary table
    print_results(args.model, results)

    # Plot results
    plot_results(args.model, problems[args.model], results)

if __name__ == "__main__":
    main()
