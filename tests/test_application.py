"""アプリの入口と assets のパスが作業ディレクトリに左右されないことを検証する。"""

from __future__ import annotations

import os
import runpy
import tempfile
import unittest
from pathlib import Path
from typing import cast
from unittest.mock import patch

from autoimmune_atlas import refresh, snapshot
from tests.ui_fixture import (
    ASSETS_PATH,
    DashApplication,
)


class ApplicationPathTests(unittest.TestCase):
    def test_data_paths_keep_the_repository_data_directory(self):
        expected = ASSETS_PATH.parent / "data" / "snapshot.json"
        self.assertEqual(snapshot.DATA_PATH, expected)
        self.assertEqual(refresh.DATA_PATH, expected)

    def test_entrypoint_and_assets_do_not_depend_on_working_directory(self):
        entrypoint = ASSETS_PATH.parent / "app.py"
        previous = Path.cwd()
        try:
            with tempfile.TemporaryDirectory() as directory:
                os.chdir(directory)
                with patch(
                    "autoimmune_atlas.snapshot.load_snapshot", return_value=None
                ):
                    namespace = cast(
                        dict[str, object],
                        runpy.run_path(str(entrypoint), run_name="app_path_test"),
                    )
            dash_app_value = namespace["app"]
            assert hasattr(dash_app_value, "server") and hasattr(
                dash_app_value, "config"
            )
            dash_app = cast(DashApplication, dash_app_value)
            assets_folder = dash_app.config["assets_folder"]
            assert isinstance(assets_folder, str)
            self.assertEqual(Path(assets_folder), ASSETS_PATH)
            client = dash_app.server.test_client()
            for asset in ("/assets/style.css", "/assets/icons/info.svg"):
                response = client.get(asset)
                try:
                    self.assertEqual(response.status_code, 200)
                finally:
                    response.close()
        finally:
            os.chdir(previous)
