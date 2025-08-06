#!/usr/bin/env python3

import subprocess
from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt
import os

EXECUTABLE = "./main"  # o il path al tuo binario C++
console = Console()

def clear():
    os.system("clear" if os.name == "posix" else "cls")

def run_solver(equation: str) -> str:
    try:
        result = subprocess.run([EXECUTABLE, equation], capture_output=True, text=True, timeout=3)
        return result.stdout.strip() if result.returncode == 0 else result.stderr.strip()
    except Exception as e:
        return f"Error during execution: {e}"

def main():
    while True:
        clear()
        console.print(Panel("[bold cyan]🧮 EQUATION SOLVER (C++)[/bold cyan]", border_style="blue"))
        equation = Prompt.ask("[bold green]Insert an equation[/bold green] (Enter to exit)")
        if not equation.strip():
            console.print("\n[bold red]Exit👋[/bold red]")
            break
        output = run_solver(equation)
        console.print(Panel(output, title=equation, border_style="green"))
        input("\nPremi invio per continuare...")

if __name__ == "__main__":
    main()
