const fs = require('fs');
require('mathjax-full/js/util/asyncLoad/node.js');
const { mathjax } = require('mathjax-full/js/mathjax.js');
const { TeX } = require('mathjax-full/js/input/tex.js');
const { SVG } = require('mathjax-full/js/output/svg.js');
const { liteAdaptor } = require('mathjax-full/js/adaptors/liteAdaptor.js');
const { RegisterHTMLHandler } = require('mathjax-full/js/handlers/html.js');
const { AllPackages } = require('mathjax-full/js/input/tex/AllPackages.js');

const adaptor = liteAdaptor({ fontSize: 16 });
RegisterHTMLHandler(adaptor);

let html = fs.readFileSync('navier-stokes-working-notes.html', 'utf8');
html = html.replace(/<script>\s*window\.MathJax[\s\S]*?<\/script>\s*/, '');
html = html.replace(/<script src="https:\/\/cdnjs[^"]*"[^>]*><\/script>\s*/, '');

(async () => {
  // Pre-parse with retries so lazy entity tables load first
  let parsed;
  await mathjax.handleRetriesFor(() => { parsed = adaptor.parse(html, 'text/html'); });
  const tex = new TeX({ packages: AllPackages,
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$']],
    processEscapes: true });
  const svg = new SVG({ fontCache: 'none' });
  let doc;
  await mathjax.handleRetriesFor(() => {
    doc = mathjax.document(parsed, { InputJax: tex, OutputJax: svg });
    doc.render();
  });
  let out = adaptor.doctype(doc.document) + "\n" + adaptor.outerHTML(adaptor.root(doc.document));
  fs.writeFileSync('navier-stokes-working-notes-rendered.html', out);
  console.log('rendered bytes:', out.length);
  console.log('math items:', (out.match(/<mjx-container/g) || []).length);
  console.log('leftover $$ count:', (out.match(/\$\$/g) || []).length);
  const err = (out.match(/data-mjx-error="([^"]*)"/g) || []);
  console.log('tex errors:', err.length, err.slice(0,5));
})().catch(e => { console.error('FAILED:', e.message); process.exit(1); });
