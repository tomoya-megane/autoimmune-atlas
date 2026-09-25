"""Autoimmune Atlas の Dash アプリを起動する。"""

import argparse
import json
from pathlib import Path

from autoimmune_atlas.snapshot import load_snapshot
from autoimmune_atlas.ui.application import create_app

BASE_DIR = Path(__file__).resolve().parent

try:
    SNAPSHOT = load_snapshot()
    LOAD_ERROR = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    SNAPSHOT = None
    LOAD_ERROR = error

app = create_app(
    SNAPSHOT,
    LOAD_ERROR,
    assets_folder=BASE_DIR / "assets",
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Autoimmune Atlas")
    parser.add_argument(
        "--dev", action="store_true", help="Automatically reload when code changes"
    )
    args = parser.parse_args()
    app.run(host="127.0.0.1", port=8050, debug=args.dev)
