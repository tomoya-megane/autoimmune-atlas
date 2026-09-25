"""Autoimmune Atlas の Dash アプリを起動する。"""

import argparse
import json
from pathlib import Path

from backend.snapshot import load_snapshot
from backend.ui.application import create_app

BASE_DIR = Path(__file__).resolve().parent

try:
    snapshot_data = load_snapshot()
    load_error = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    snapshot_data = None
    load_error = error

app = create_app(
    snapshot_data,
    load_error,
    assets_folder=BASE_DIR / "assets",
)


if __name__ == "__main__":

    class Arguments(argparse.Namespace):
        dev: bool = False

    parser = argparse.ArgumentParser(description="Run Autoimmune Atlas")
    parser.add_argument(
        "--dev", action="store_true", help="Automatically reload when code changes"
    )
    args = parser.parse_args(namespace=Arguments())
    app.run(  # pyright: ignore[reportUnknownMemberType] -- Dash の可変長引数が Unknown。
        host="127.0.0.1", port=8050, debug=args.dev
    )
