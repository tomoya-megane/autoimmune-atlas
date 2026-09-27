"""Dash アプリを組み立てる。"""

from pathlib import Path

from dash import Dash, Input, Output, dcc, html

from autoimmune_atlas.models import Snapshot
from autoimmune_atlas.ui.callbacks import register_callbacks
from autoimmune_atlas.ui.genetics_callbacks import register_genetics_callbacks
from autoimmune_atlas.ui.genetics_layout import genetics_page, genetics_unavailable_page
from autoimmune_atlas.ui.layout import dashboard_layout, unavailable_layout

NO_SNAPSHOT = "The snapshot is not available, so the genetics page cannot be shown."


def create_app(
    snapshot: Snapshot | None,
    error: Exception | None = None,
    *,
    assets_folder: str | Path,
) -> Dash:
    """保存済みデータまたは明示的な fixture から、2 ページの Dash アプリを作る。"""
    application = Dash(
        __name__,
        title="Autoimmune Atlas",
        suppress_callback_exceptions=True,
        assets_folder=str(assets_folder),
    )
    if snapshot is None:
        drug_page = unavailable_layout(error)
        gene_page = genetics_unavailable_page(NO_SNAPSHOT)
    else:
        drug_page = dashboard_layout(snapshot)
        gene_page = genetics_page(snapshot)
        register_genetics_callbacks(application, snapshot)
    # 初期表示は薬剤ページにし、URL の callback が正しいほうを出す。
    gene_page.hidden = True  # pyright: ignore[reportAttributeAccessIssue] - Dash の型定義に hidden がない。
    application.layout = html.Div(
        [dcc.Location(id="url", refresh=False), drug_page, gene_page]
    )
    if snapshot is not None:
        register_callbacks(application, snapshot)

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("drug-page", "hidden"),
        Output("genetics-page", "hidden"),
        Input("url", "pathname"),
    )
    def switch_page(pathname: str | None) -> tuple[bool, bool]:
        """/genetics なら遺伝子ページを、それ以外は薬剤ページを出す。"""
        on_genetics = pathname == "/genetics"
        return on_genetics, not on_genetics

    return application
