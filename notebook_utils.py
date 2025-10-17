# notebook_generator.py
import json
import datetime
from pathlib import Path

def make_notebook(csv_path: Path):
    nb = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# Auto Benchmark Notebook\n",
                    f"_Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}_\n",
                    "\nRun the cells below to visualize your benchmark results.\n"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "import pandas as pd, matplotlib.pyplot as plt\n",
                    f"df = pd.read_csv(r'{csv_path.resolve()}')\n",
                    "df"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "plt.figure(figsize=(6,4))\n",
                    "df_plot = df.sort_values('avg_s') if 'avg_s' in df.columns else df\n",
                    "plt.bar(df_plot['script'], df_plot['avg_s'])\n",
                    "plt.ylabel('Average runtime [s]')\n",
                    "plt.title('Average runtime per script')\n",
                    "plt.xticks(rotation=20, ha='right')\n",
                    "plt.grid(axis='y', linestyle='--', alpha=0.5)\n",
                    "plt.tight_layout()\n",
                    "plt.show()\n"
                ]
            }
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"}
        },
        "nbformat": 4,
        "nbformat_minor": 5
    }
    nb_path = csv_path.parent / "Benchmark_Auto.ipynb"
    with open(nb_path, "w", encoding="utf-8") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
    print(f"\033[1;33m[auto-notebook]\033[0m Created {nb_path.name} for Jupyter visualization.")
