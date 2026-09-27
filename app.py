"""Autoimmune Atlas の Dash アプリを起動する。"""

import argparse
import json
import sys
from pathlib import Path

from autoimmune_atlas.genetics import load_genetics
from autoimmune_atlas.snapshot import load_snapshot
from autoimmune_atlas.ui.application import create_app

BASE_DIR = Path(__file__).resolve().parent

try:
    snapshot_data = load_snapshot()
    load_error = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    snapshot_data = None
    load_error = error

# genetics.json が読めなくても、薬剤ページは止めない。
try:
    genetics_data = load_genetics()
    genetics_error = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    print(f"genetics.json を読めません: {error}", file=sys.stderr)
    genetics_data = None
    genetics_error = error

app = create_app(
    snapshot_data,
    load_error,
    assets_folder=BASE_DIR / "assets",
    genetics=genetics_data,
    genetics_error=genetics_error,
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
