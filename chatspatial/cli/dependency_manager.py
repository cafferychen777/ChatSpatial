"""
ChatSpatial Dependency Management CLI

A user-friendly command-line interface for managing dependencies.
This tool makes dependency management from a pain point into a competitive advantage.

Commands:
    chatspatial deps check       - Check dependency status
    chatspatial deps doctor      - Diagnose and fix common issues
    chatspatial deps install     - Install dependencies by level
    chatspatial deps features    - Show available features
    chatspatial deps safe-mode   - Run in safe mode (core deps only)

Design Philosophy:
- Zero configuration for basic usage
- Progressive disclosure of complexity
- Clear, actionable error messages
- Never leave users stranded

Author: Linus-approved dependency CLI
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import click
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.prompt import Confirm
from rich.syntax import Syntax
from rich.table import Table

# Handle rich import gracefully
try:
    from rich import box
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.prompt import Confirm
    from rich.syntax import Syntax
    from rich.table import Table

    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False

    # Fallback console class
    class Console:
        def print(self, *args, **kwargs):
            print(*args)

        def rule(self, *args, **kwargs):
            print("-" * 50)


console = Console()

# Import our dependency management system
try:
    from ..utils.dependency_fallbacks import (enable_safe_mode,
                                              list_available_fallbacks)
    from ..utils.smart_import import (DependencyLevel, check_environment,
                                      get_smart_importer)
except ImportError:
    # Handle case where we're running outside the package
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from chatspatial.utils.dependency_fallbacks import (
        list_available_fallbacks)
    from chatspatial.utils.smart_import import (DependencyLevel,
                                                check_environment,
                                                get_smart_importer)


@click.group()
def deps():
    """ChatSpatial Dependency Management"""
    pass


@deps.command()
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "json", "compact"]),
    default="table",
    help="Output format",
)
@click.option(
    "--level",
    type=click.Choice(["core", "recommended", "advanced", "experimental", "all"]),
    default="all",
    help="Show specific dependency level",
)
@click.option("--missing-only", is_flag=True, help="Show only missing dependencies")
def check(output_format: str, level: str, missing_only: bool):
    """Check dependency status"""

    if RICH_AVAILABLE:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Checking dependencies...", total=None)
            env_report = check_environment()
    else:
        print("Checking dependencies...")
        env_report = check_environment()

    dep_status = env_report["dependency_status"]

    if output_format == "json":
        click.echo(json.dumps(env_report, indent=2))
        return
    elif output_format == "compact":
        _print_compact_status(dep_status, level, missing_only)
        return

    # Rich table format (default)
    _print_detailed_status(env_report, level, missing_only)


def _print_compact_status(
    dep_status: Dict[str, Dict[str, Any]], level: str, missing_only: bool
):
    """Print compact status without rich formatting"""
    levels_to_show = (
        [level]
        if level != "all"
        else ["core", "recommended", "advanced", "experimental"]
    )

    for level_name in levels_to_show:
        if level_name not in dep_status:
            continue

        deps = dep_status[level_name]
        print(f"\n{level_name.upper()} Dependencies:")
        print("-" * (len(level_name) + 15))

        for name, info in deps.items():
            if missing_only and info["available"]:
                continue

            status = "✅" if info["available"] else "❌"
            print(f"{status} {name:<20} - {info['description']}")

            if not info["available"]:
                print(f"   Install: {info['install_cmd']}")


def _print_detailed_status(env_report: Dict[str, Any], level: str, missing_only: bool):
    """Print detailed status with rich formatting"""

    if not RICH_AVAILABLE:
        _print_compact_status(env_report["dependency_status"], level, missing_only)
        return

    dep_status = env_report["dependency_status"]
    stats = env_report["statistics"]

    # Header
    console.rule("[bold blue]ChatSpatial Dependency Status")

    # Statistics panel
    stats_text = "\n".join(
        [
            f"Core: {stats['core']}",
            f"Recommended: {stats['recommended']}",
            f"Advanced: {stats['advanced']}",
            f"Experimental: {stats['experimental']}",
        ]
    )

    console.print(Panel(stats_text, title="[bold]Summary", border_style="green"))

    # Critical issues
    if env_report["critical_issues"]:
        issues_text = "\n".join(env_report["critical_issues"])
        console.print(
            Panel(issues_text, title="[bold red]Critical Issues", border_style="red")
        )

    # Detailed tables
    levels_to_show = (
        [level]
        if level != "all"
        else ["core", "recommended", "advanced", "experimental"]
    )

    for level_name in levels_to_show:
        if level_name not in dep_status:
            continue

        deps = dep_status[level_name]

        if missing_only:
            deps = {k: v for k, v in deps.items() if not v["available"]}

        if not deps:
            continue

        table = Table(
            title=f"{level_name.capitalize()} Dependencies",
            box=box.ROUNDED,
            show_header=True,
            header_style="bold magenta",
        )

        table.add_column("Package", style="cyan", no_wrap=True)
        table.add_column("Status", justify="center")
        table.add_column("Description", style="dim")
        table.add_column("Install Command", style="green")

        for name, info in deps.items():
            status = "✅ Available" if info["available"] else "❌ Missing"
            install_cmd = info["install_cmd"] if not info["available"] else ""

            table.add_row(name, status, info["description"], install_cmd)

        console.print(table)


@deps.command()
@click.option("--fix", is_flag=True, help="Automatically fix detected issues")
@click.option("--verbose", is_flag=True, help="Show detailed diagnostic information")
def doctor(fix: bool, verbose: bool):
    """Diagnose and fix common dependency issues"""

    console.rule("[bold green]ChatSpatial Doctor")

    if RICH_AVAILABLE:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running diagnostics...", total=None)
            issues = _diagnose_issues()
    else:
        print("Running diagnostics...")
        issues = _diagnose_issues()

    if not issues:
        console.print("🎉 [green]No issues detected! Your environment looks good.")
        return

    console.print(f"[yellow]Found {len(issues)} issues:[/yellow]")

    for i, issue in enumerate(issues, 1):
        console.print(f"\n[bold]{i}. {issue['title']}[/bold]")
        console.print(f"   {issue['description']}")

        if issue.get("solution"):
            console.print(f"   [green]Solution: {issue['solution']}[/green]")

        if issue.get("command"):
            console.print(f"   [blue]Command: {issue['command']}[/blue]")

            if fix and issue.get("auto_fixable", False):
                if Confirm.ask(f"Apply fix for: {issue['title']}?"):
                    _apply_fix(issue)

    if not fix:
        console.print(
            "\n[dim]Run with --fix to automatically resolve fixable issues[/dim]"
        )


def _diagnose_issues() -> List[Dict[str, Any]]:
    """Diagnose common dependency issues"""
    issues = []

    env_report = check_environment()

    # Check for critical missing dependencies
    for name, info in env_report["dependency_status"]["core"].items():
        if not info["available"]:
            issues.append(
                {
                    "title": f"Missing core dependency: {name}",
                    "description": f"Core functionality will be broken without {name}",
                    "solution": f"Install {name}",
                    "command": info["install_cmd"],
                    "severity": "critical",
                    "auto_fixable": True,
                }
            )

    # Check for common version conflicts
    issues.extend(_check_version_conflicts())

    # Check for platform-specific issues
    issues.extend(_check_platform_issues())

    # Check for environment issues
    issues.extend(_check_environment_issues())

    return issues


def _check_version_conflicts() -> List[Dict[str, Any]]:
    """Check for known version conflicts"""
    issues = []

    # Check for dask version conflict with squidpy
    try:
        import dask

        dask_version = tuple(map(int, dask.__version__.split(".")[:2]))
        if dask_version >= (2025, 1):
            issues.append(
                {
                    "title": "Dask version conflict",
                    "description": f"Dask {dask.__version__} may cause squidpy import issues",
                    "solution": "Downgrade dask to 2024.x.x",
                    "command": "pip install 'dask<2025'",
                    "severity": "warning",
                    "auto_fixable": True,
                }
            )
    except ImportError:
        pass

    # Check for TensorFlow/PyTorch conflicts
    try:
        import tensorflow as tf
        import torch

        # Add specific version conflict checks here
    except ImportError:
        pass

    return issues


def _check_platform_issues() -> List[Dict[str, Any]]:
    """Check for platform-specific issues"""
    issues = []

    # Check for R installation for rpy2
    importer = get_smart_importer()
    rpy2_result = importer.smart_import("rpy2", required=False)

    if not rpy2_result.success and "rpy2" in str(rpy2_result.error).lower():
        issues.append(
            {
                "title": "R not found for rpy2",
                "description": "rpy2 requires R to be installed on your system",
                "solution": "Install R from https://cran.r-project.org/",
                "severity": "info",
                "auto_fixable": False,
            }
        )

    return issues


def _check_environment_issues() -> List[Dict[str, Any]]:
    """Check for environment setup issues"""
    issues = []

    # Check Python version
    if sys.version_info < (3, 8):
        issues.append(
            {
                "title": "Python version too old",
                "description": f"Python {sys.version} is not supported. Need Python 3.8+",
                "solution": "Upgrade Python to 3.8 or newer",
                "severity": "critical",
                "auto_fixable": False,
            }
        )

    # Check for pip
    try:
        import pip
    except ImportError:
        issues.append(
            {
                "title": "pip not available",
                "description": "pip is required for installing dependencies",
                "solution": "Install pip",
                "severity": "critical",
                "auto_fixable": False,
            }
        )

    return issues


def _apply_fix(issue: Dict[str, Any]) -> None:
    """Apply an automatic fix"""
    if not issue.get("command"):
        console.print("❌ [red]No command available for this fix[/red]")
        return

    try:
        import subprocess

        result = subprocess.run(
            issue["command"].split(), check=True, capture_output=True, text=True
        )
        console.print("✅ [green]Fix applied successfully[/green]")
    except subprocess.CalledProcessError as e:
        console.print(f"❌ [red]Fix failed: {e}[/red]")


@deps.command()
@click.argument(
    "level", type=click.Choice(["core", "recommended", "advanced", "experimental"])
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be installed without installing"
)
@click.option(
    "--force", is_flag=True, help="Force installation even if already installed"
)
def install(level: str, dry_run: bool, force: bool):
    """Install dependencies by level"""

    console.rule(f"[bold blue]Installing {level} dependencies")

    importer = get_smart_importer()

    if dry_run:
        console.print("[dim]Dry run - showing what would be installed:[/dim]")

    try:
        dependency_level = DependencyLevel(level)
        results = importer.install_dependencies(dependency_level, dry_run=dry_run)

        if not dry_run:
            # Show results
            successful = sum(1 for success in results.values() if success)
            total = len(results)

            if successful == total:
                console.print(
                    f"🎉 [green]Successfully installed all {total} packages![/green]"
                )
            else:
                console.print(
                    f"⚠️  [yellow]Installed {successful}/{total} packages[/yellow]"
                )

                # Show failed installations
                failed = [name for name, success in results.items() if not success]
                if failed:
                    console.print(f"[red]Failed to install: {', '.join(failed)}[/red]")

    except Exception as e:
        console.print(f"❌ [red]Installation failed: {e}[/red]")


@deps.command()
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "json"]),
    default="table",
    help="Output format",
)
def features(output_format: str):
    """Show available features based on installed dependencies"""

    importer = get_smart_importer()
    features = importer.get_available_features()
    fallbacks = list_available_fallbacks()

    if output_format == "json":
        click.echo(json.dumps({"features": features, "fallbacks": fallbacks}, indent=2))
        return

    console.rule("[bold blue]Available Features")

    if RICH_AVAILABLE:
        table = Table(
            title="Feature Matrix",
            box=box.ROUNDED,
            show_header=True,
            header_style="bold magenta",
        )

        table.add_column("Category", style="cyan", no_wrap=True)
        table.add_column("Available Methods", style="green")
        table.add_column("Fallbacks", style="yellow")

        for category, methods in features.items():
            methods_str = ", ".join(methods) if methods else "None"
            category_fallbacks = [f for f in fallbacks if category in f]
            fallbacks_str = (
                ", ".join(category_fallbacks) if category_fallbacks else "None"
            )

            table.add_row(category, methods_str, fallbacks_str)

        console.print(table)
    else:
        for category, methods in features.items():
            print(f"{category}: {', '.join(methods) if methods else 'None'}")


@deps.command()
@click.option(
    "--config-only", is_flag=True, help="Only generate safe mode configuration"
)
def safe_mode(config_only: bool):
    """Enable safe mode (core dependencies only)"""

    console.rule("[bold yellow]ChatSpatial Safe Mode")

    console.print(
        """
    🛡️  Safe Mode enables basic functionality using only core dependencies.
    
    Available in Safe Mode:
    • Basic preprocessing (scanpy)
    • Clustering and UMAP
    • Basic visualization
    • Differential expression analysis
    • Simple spatial plots (if coordinates available)
    
    Not Available:
    • Advanced spatial analysis methods
    • Cell-cell communication analysis
    • Spatial domain identification
    • Advanced deconvolution methods
    """
    )

    if config_only:
        config = {
            "mode": "safe",
            "enabled_features": ["preprocessing", "clustering", "basic_visualization"],
            "disabled_features": [
                "spatial_analysis",
                "cell_communication",
                "deconvolution",
            ],
        }

        config_path = Path.cwd() / "chatspatial_safe_config.json"
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)

        console.print(f"✅ Safe mode configuration saved to {config_path}")
    else:
        console.print("To use safe mode in your analysis:")

        if RICH_AVAILABLE:
            code = """
# Enable safe mode
from chatspatial.utils.dependency_fallbacks import enable_safe_mode

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Enable safe mode
safe_mode = enable_safe_mode(adata)

# Run basic analysis
safe_mode.basic_preprocessing()
safe_mode.basic_visualization()

# Check what's available
capabilities = safe_mode.get_capabilities()
print(capabilities)
"""

            syntax = Syntax(code, "python", theme="monokai", line_numbers=True)
            console.print(Panel(syntax, title="Example Usage", border_style="green"))
        else:
            print(
                """
Example usage:

from chatspatial.utils.dependency_fallbacks import enable_safe_mode

adata = sc.read_h5ad('your_data.h5ad')
safe_mode = enable_safe_mode(adata)
safe_mode.basic_preprocessing()
safe_mode.basic_visualization()
"""
            )


if __name__ == "__main__":
    deps()
