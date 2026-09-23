// 説明をボタンに接して置き、画面の端では横位置を収める。
function showHelp(event) {
    const tip = event.target.closest?.(".info-tip");
    if (!tip || (event.type === "pointerenter" && event.target !== tip)) return;
    if (event.type !== "reposition") tip.classList.remove("dismissed");
    const button = tip.querySelector(".info-button").getBoundingClientRect();
    const panel = tip.querySelector(".info-content");
    const box = panel.getBoundingClientRect();
    panel.style.left = `${Math.max(12, Math.min(button.left, window.innerWidth - box.width - 12))}px`;
    panel.style.top = `${button.bottom + box.height <= window.innerHeight - 12 ? button.bottom : Math.max(12, button.top - box.height)}px`;
}

document.addEventListener("pointerenter", showHelp, true);
document.addEventListener("focusin", showHelp, true);
function repositionVisibleHelp() {
    document.querySelectorAll(".info-tip:hover, .info-tip:focus-within").forEach((tip) => showHelp({type: "reposition", target: tip}));
}
document.addEventListener("scroll", repositionVisibleHelp, true);
window.addEventListener("resize", repositionVisibleHelp);
document.addEventListener("keydown", (event) => {
    if (event.key !== "Escape") return;
    document.querySelectorAll(".info-tip:hover, .info-tip:focus-within").forEach((tip) => tip.classList.add("dismissed"));
});
