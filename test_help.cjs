// node test_help.cjs で、ブラウザーや追加ライブラリなしにイベント処理を確認する。
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const handlers = {};
const classes = new Set();
let buttonBox = {left: 920, top: 500, bottom: 522};
let panelBox = {width: 340, height: 160};
const panel = {style: {}, getBoundingClientRect: () => panelBox};
const button = {getBoundingClientRect: () => buttonBox, closest: () => tip};
const tip = {
    closest: () => tip,
    classList: {add: name => classes.add(name), remove: name => classes.delete(name)},
    querySelector: selector => selector === ".info-button" ? button : panel,
};
const context = {
    document: {addEventListener: (name, handler) => handlers[name] = handler, querySelectorAll: () => [tip]},
    window: {innerWidth: 1000, innerHeight: 600, addEventListener: (name, handler) => handlers[name] = handler},
};
vm.runInNewContext(fs.readFileSync(path.join(__dirname, "assets/help.js"), "utf8"), context);

// 右端かつ下端のボタンでは、説明を左へ寄せて上側に置く。
handlers.pointerenter({type: "pointerenter", target: tip});
assert.equal(panel.style.left, "648px");
assert.equal(panel.style.top, "340px");

// Escape の後は、内側へのポインター移動や画面移動だけで再表示しない。
handlers.keydown({key: "Escape"});
handlers.pointerenter({type: "pointerenter", target: button});
assert(classes.has("dismissed"));
buttonBox = {left: 300, top: 100, bottom: 122};
panelBox = {width: 351, height: 140};
context.window.innerWidth = 375;
handlers.resize();
handlers.scroll();
assert(classes.has("dismissed"));
assert.equal(panel.style.left, "12px");
assert.equal(panel.style.top, "122px");

// 新たなフォーカスや外側からのポインター移動なら説明を再表示する。
handlers.focusin({type: "focusin", target: button});
assert(!classes.has("dismissed"));
handlers.keydown({key: "Escape"});
handlers.pointerenter({type: "pointerenter", target: tip});
assert(!classes.has("dismissed"));
console.log("Tooltip interaction checks passed (simulated DOM).");
