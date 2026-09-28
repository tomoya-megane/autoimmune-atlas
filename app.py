"""Autoimmune Atlas の Dash アプリを起動する。"""

import argparse
import json
from pathlib import Path

from autoimmune_atlas.snapshot import load_snapshot
from autoimmune_atlas.ui.application import create_app

BASE_DIR = Path(__file__).resolve().parent

try:
    snapshot_data = load_snapshot()
    load_error = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    snapshot_data = None
    load_error = error

app = create_app(snapshot_data, load_error, assets_folder=BASE_DIR / "assets")


if __name__ == "__main__":

    class Arguments(argparse.Namespace):
        dev: bool = False
        host: str = "127.0.0.1"
        port: int = 8050

    parser = argparse.ArgumentParser(description="Run Autoimmune Atlas")
    parser.add_argument(
        "--dev", action="store_true", help="Automatically reload when code changes"
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Address to listen on. Use 0.0.0.0 to accept other machines",
    )
    parser.add_argument("--port", type=int, default=8050, help="Port to listen on")
    args = parser.parse_args(namespace=Arguments())
    if args.dev and args.host not in ("127.0.0.1", "localhost"):
        # 開発モードは Werkzeug のデバッガを公開するので、外から繋がる host では起動しない。
        parser.error("--dev は --host 127.0.0.1 でしか使えない")
    app.run(  # pyright: ignore[reportUnknownMemberType] -- Dash の可変長引数が Unknown。
        host=args.host, port=args.port, debug=args.dev
    )
