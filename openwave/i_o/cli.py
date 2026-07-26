"""
OpenWave CLI - Command Line Interface for running Xperiments.

This module provides the command-line entry point for OpenWave,
allowing users to interactively select and run Xperiments.
"""

import argparse
import os
import sys
import subprocess
import webbrowser
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.text import Text

# Conditional import for simple_term_menu (not available on Windows)
try:
    from simple_term_menu import TerminalMenu

    HAS_INTERACTIVE_MENU = True
except (ImportError, NotImplementedError):
    HAS_INTERACTIVE_MENU = False

console = Console()


def _get_version():
    """Get the package version string."""
    try:
        from openwave import __version__

        return __version__
    except ImportError:
        from importlib.metadata import version

        return version("OPENWAVE")


# Hardcoded tutorial entry
TUTORIAL_URL = "https://github.com/openwave-labs/openwave/blob/main/TUTORIAL.md"
TUTORIAL_ENTRY = ("OpenWave TUTORIAL", TUTORIAL_URL)


def get_experiments_list():
    """
    Get a list of available Xperiment launchers from the xperiments directory.
    Each collection has exactly one _launcher.py file.
    Collections are subdirectories under xperiments/.
    Collections can be Xperiment methods or groups.


    Returns:
        list: List of tuples containing (display_name, file_path)
    """
    # Get the xperiments directory path
    # Navigate from i_o module to parent package, then to xperiments
    package_dir = Path(__file__).parent.parent
    xperiments_dir = package_dir / "xperiments"

    if not xperiments_dir.exists():
        print(f"Error: Xperiments directory not found at {xperiments_dir}")
        sys.exit(1)

    experiments = []

    # Find all _launcher.py files in collection subdirectories
    for launcher_path in xperiments_dir.rglob("_launcher.py"):
        # Skip files in directories that start with underscore
        if any(part.startswith("_") for part in launcher_path.parent.parts):
            continue

        # Get collection directory
        collection_dir = launcher_path.parent
        collection_name = collection_dir.name

        # Skip if directly in xperiments root (no collection)
        if collection_dir == xperiments_dir:
            continue

        # Get display name from collection's __init__.py docstring
        display_name = None
        init_path = collection_dir / "__init__.py"

        if init_path.exists():
            try:
                with open(init_path, "r") as f:
                    lines = f.readlines()
                    # Look for first non-empty line in docstring
                    if len(lines) > 0 and '"""' in lines[0]:
                        for line in lines[1:]:
                            if '"""' in line:
                                break
                            stripped = line.strip()
                            if stripped:
                                display_name = stripped
                                break
            except Exception:
                pass

        # Fallback to formatted collection name
        if not display_name:
            display_name = (
                collection_name.replace("__", ": ").replace("_", " ").replace("-", " ").title()
            )

        experiments.append((display_name, str(launcher_path)))

    # Sort by display name (which starts with A/, B/, C/ etc.)
    experiments.sort(key=lambda x: x[0])

    # Insert hardcoded tutorial entry at the beginning
    experiments.insert(0, TUTORIAL_ENTRY)

    return experiments


def show_menu_simple(experiments):
    """
    Display a simple numbered menu for Xperiment selection.
    Fallback for systems where interactive menu is not available.

    Args:
        experiments: List of tuples containing (display_name, file_path)

    Returns:
        tuple: (display_name, file_path) of the selected xperiment
    """
    pkg_version = _get_version()

    console.print()
    console.print(
        Panel(
            Text(f"OPENWAVE  v{pkg_version}", justify="center", style="bold"),
            subtitle="Available XPERIMENTS",
            border_style="green",
            width=64,
        )
    )
    console.print()

    # Display numbered list of experiments with spacing
    for idx, (display_name, _) in enumerate(experiments, 1):
        console.print(f"  [bold]{idx}.[/bold] {display_name}")
        # console.print()  # Spacing between items

    console.print(f"  [dim]{len(experiments) + 1}. EXIT[/dim]")
    console.print()

    while True:
        try:
            choice = input("\nSelect an Xperiment (enter number): ").strip()
            choice_num = int(choice)

            if choice_num == len(experiments) + 1:
                print("Exiting...")
                sys.exit(0)

            if 1 <= choice_num <= len(experiments):
                return experiments[choice_num - 1]
            else:
                print(
                    f"Invalid choice. Please enter a number between 1 and {len(experiments) + 1}"
                )
        except ValueError:
            print("Invalid input. Please enter a number.")
        except KeyboardInterrupt:
            print("\n\nExiting...")
            sys.exit(0)


def show_menu_interactive(experiments):
    """
    Display an interactive menu using arrow keys for Xperiment selection.
    Only available on Unix-like systems (Linux, macOS).

    Args:
        experiments: List of tuples containing (display_name, file_path)

    Returns:
        tuple: (display_name, file_path) of the selected xperiment
    """
    if not HAS_INTERACTIVE_MENU:
        # Fallback to simple menu if interactive menu not available
        return show_menu_simple(experiments)

    pkg_version = _get_version()

    # Build plain-text banner — TerminalMenu miscounts ANSI escapes as visible chars
    w = 64
    heading = f"OPENWAVE  v{pkg_version}"
    subtitle = "↑/↓ navigate  ·  ENTER to select XPERIMENT"
    title_str = (
        "\n"
        f"  ╭{'─' * (w - 2)}╮\n"
        f"  │{heading:^{w - 2}}│\n"
        f"  ╰{'─' * (w - 2)}╯\n"
        f"  {subtitle:^{w}}\n"
    )

    # Build menu options with blank lines between items for readability
    menu_options = []
    for display_name, _ in experiments:
        menu_options.append(f"  • {display_name}")
        # menu_options.append(None)  # Spacing between items

    menu_options.append("  ─── EXIT ───")
    exit_idx = len(menu_options) - 1

    terminal_menu = TerminalMenu(
        menu_options,
        title=title_str,
        menu_cursor="▶ ",
        menu_cursor_style=("fg_green", "bold"),
        menu_highlight_style=("bg_green", "fg_black"),
        cycle_cursor=True,
        skip_empty_entries=True,
    )

    choice_idx = terminal_menu.show()

    if choice_idx is None or choice_idx == exit_idx:
        print("Exiting...")
        sys.exit(0)

    return experiments[choice_idx]


def run_experiment(display_name, file_path):
    """
    Run the selected Xperiment file or open a URL.

    Args:
        display_name: Display name of the xperiment
        file_path: Path to the xperiment Python file, or a URL

    Returns:
        int: The return code from the experiment process
    """
    # Handle tutorial URL specially
    if file_path == TUTORIAL_URL:
        console.print()
        console.print(
            Panel("Opening the tutorial in your default browser...", border_style="cyan", width=64)
        )
        console.print()
        webbrowser.open(TUTORIAL_URL)
        console.print("[green]Done![/green]")
        return 0

    console.print()
    console.print(
        Panel(
            f"[bold]{display_name}[/bold]",
            title="Running XPERIMENT",
            border_style="cyan",
            width=64,
        )
    )
    console.print()

    try:
        # Run the xperiment using subprocess
        result = subprocess.run(
            [sys.executable, file_path],
            env=os.environ.copy(),
        )
        return result.returncode
    except KeyboardInterrupt:
        print("\n\nExperiment interrupted by user.")
        return 0
    except Exception as e:
        print(f"\nError running XPERIMENT: {e}")
        return 1


def main():
    """
    Main entry point for the OpenWave CLI.

    This function is called when running 'openwave -x' from the command line.
    Runs the selected xperiment and exits when it closes.
    """
    # Get list of available experiments
    experiments = get_experiments_list()

    if not experiments:
        print("No Xperiments found in the xperiments directory.")
        sys.exit(1)

    # Show interactive menu and get user selection
    display_name, file_path = show_menu_interactive(experiments)

    # Run the selected xperiment
    returncode = run_experiment(display_name, file_path)

    # Exit after xperiment closes
    style = "green" if returncode == 0 else "red"
    console.print()
    console.print(
        Panel(
            f"XPERIMENT closed (exit code: {returncode})",
            border_style=style,
            width=64,
        )
    )
    console.print()

    sys.exit(returncode)


def cli_main():
    """
    Main entry point for the 'openwave' command.

    Handles command-line arguments and routes to appropriate functionality.
    """
    parser = argparse.ArgumentParser(
        description="OpenWave - Subatomic Simulator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-x",
        "--xperiments",
        action="store_true",
        help="launch the xperiments selector",
    )

    args = parser.parse_args()

    if args.xperiments:
        main()
    else:
        parser.print_help()
        sys.exit(0)


if __name__ == "__main__":
    cli_main()
