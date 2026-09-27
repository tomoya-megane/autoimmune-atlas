// ページの hidden が切り替わったら、Plotly に幅を測り直させる。
// hidden の中で描かれた図は幅 0 のままで、窓の大きさが変わるまで描き直されない。
// dcc.Location の変化は popstate では拾えないので、hidden 属性の変化を監視する。
function watchPages() {
    const pages = ["drug-page", "genetics-page"].map((id) => document.getElementById(id)).filter(Boolean);
    if (pages.length < 2) return false;
    const observer = new MutationObserver(() => window.dispatchEvent(new Event("resize")));
    pages.forEach((page) => observer.observe(page, {attributes: true, attributeFilter: ["hidden"]}));
    return true;
}

const bootstrap = new MutationObserver(() => {
    if (watchPages()) bootstrap.disconnect();
});
bootstrap.observe(document.body, {childList: true, subtree: true});
