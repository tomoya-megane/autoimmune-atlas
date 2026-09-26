"""Dash アプリを組み立てる。"""

from pathlib import Path

from dash import Dash

from backend.models import Snapshot
from backend.ui.callbacks import register_callbacks
from backend.ui.layout import dashboard_layout, unavailable_layout


def create_app(
    snapshot: Snapshot | None,
    error: Exception | None = None,
    *,
    assets_folder: str | Path,
) -> Dash:
    """保存済みデータまたは明示的な fixture から Dash アプリを作る。"""
    # 元記録の操作部は、選択した詳細を描画するときに追加する。
    application = Dash(
        __name__,
        title="Autoimmune Atlas",
        suppress_callback_exceptions=True,
        assets_folder=str(assets_folder),
    )
    application.layout = (
        unavailable_layout(error) if snapshot is None else dashboard_layout(snapshot)
    )
    if snapshot is not None:
        register_callbacks(application, snapshot)
    return application
