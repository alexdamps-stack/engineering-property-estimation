// Re-author the upstream grouped vectors as standalone, editable SVG assets.
import { chromium } from 'playwright';
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
const root = new URL('../../docs/design-reference/iso-10628/', import.meta.url);
const fontData = readFileSync(new URL('sources/font/osifont.ttf', root)).toString('base64');
const fontStyle = `<style>@font-face{font-family:osifont;src:url(data:font/ttf;base64,${fontData}) format('truetype');font-weight:normal;font-style:normal}</style>`;
const save = (path, data) => writeFileSync(new URL(path, root), data);
mkdirSync(new URL('symbols/', root), { recursive: true });
const browser = await chromium.launch({ channel: 'chrome', headless: true });
const page = await browser.newPage();
const symbols = [];
const sheets = [];
try {
  for (let sheet = 1; sheet <= 7; sheet++) {
    await page.goto(new URL(`sources/sheet-${sheet}.svg`, root).href);
    await page.evaluate(async data => {
      const font = new FontFace('osifont', `url(data:font/ttf;base64,${data})`);
      await font.load(); document.fonts.add(font); await document.fonts.ready;
    }, fontData);
    const result = await page.evaluate(() => {
      const svg = document.documentElement;
      const serializer = new XMLSerializer();
      const defs = svg.querySelector('defs');
      const entries = [];
      for (const label of svg.querySelectorAll('text')) {
        const text = label.textContent.replace(/\s+/g, ' ').trim();
        const match = text.match(/^REG#:\s*(\S+?)\s*DESC:\s*(.+)$/);
        if (!match) continue;
        // Two upstream groups occur out of sequence; resolve by their stable IDs.
        const overrides = { X8179: 'g35204', C2005: 'g4862' };
        const drawing = overrides[match[1]] ? svg.querySelector('#'+overrides[match[1]]) : label.nextElementSibling;
        // Each sheet also has a legend example without a following drawing group.
        if (drawing?.localName !== 'g') continue;
        const box = drawing.getBBox();
        if (box.width + box.height <= 0) throw new Error(`Empty symbol: ${text}`);
        const ancestors = [];
        for (let p = drawing.parentElement; p !== svg; p = p.parentElement) {
          const shell = p.cloneNode(false);
          shell.removeAttribute('id');
          ancestors.unshift(shell);
        }
        const clone = drawing.cloneNode(true);
        let wrapped = clone;
        for (let i = ancestors.length - 1; i >= 0; i--) {
          ancestors[i].append(wrapped);
          wrapped = ancestors[i];
        }
        // Bounds in document user coordinates, including all nested transforms.
        const m = svg.getScreenCTM().inverse().multiply(drawing.getScreenCTM());
        const corners = [[box.x,box.y],[box.x+box.width,box.y],[box.x,box.y+box.height],[box.x+box.width,box.y+box.height]].map(([x,y]) => new DOMPoint(x,y).matrixTransform(m));
        const x = Math.min(...corners.map(p=>p.x));
        const y = Math.min(...corners.map(p=>p.y));
        const w = Math.max(...corners.map(p=>p.x))-x;
        const h = Math.max(...corners.map(p=>p.y))-y;
        // Generous padding also contains stroke caps, joins and arrow markers.
        const pad = Math.max(2, Math.max(w,h)*0.08);
        const viewBox = [x-pad,y-pad,w+2*pad,h+2*pad].map(n=>Number(n.toFixed(6)));
        entries.push({ registration: match[1], description: match[2], sourceGroup: drawing.id, sourceLabel: label.id,
          sourceLayer: drawing.parentElement.id, viewBox,
          defs: defs ? serializer.serializeToString(defs) : '', geometry: serializer.serializeToString(wrapped) });
      }
      return { title: svg.querySelector('title')?.textContent, entries };
    });
    sheets.push({ sheet, title: result.title, count: result.entries.length });
    for (const entry of result.entries) {
      const { defs, geometry, ...meta } = entry;
      const slug = entry.description.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g,'').slice(0,100).replace(/-$/,'');
      const id = `s${sheet}-${entry.registration.toLowerCase()}-${slug}`;
      const esc = s => s.replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('"','&quot;');
      const out = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="${entry.viewBox.join(' ')}" width="240" height="240" role="img" aria-labelledby="asset-title">\n<title id="asset-title">${esc(entry.registration+' — '+entry.description)}</title>\n<desc>Adapted from Steven Baltakatei Sandoval's ISO 10628-2:2012 symbol sheets. CC BY-SA 4.0; https://creativecommons.org/licenses/by-sa/4.0/. Standalone vector extraction; geometry retained. Source sheet ${sheet}, group ${entry.sourceGroup}.</desc>\n${defs}\n${geometry}\n</svg>\n`;
      save(`symbols/${id}.svg`, geometry.includes('<text') ? out.replace(defs, defs+fontStyle) : out);
      symbols.push({ id, sheet, ...meta, file: `symbols/${id}.svg` });
    }
    console.log(`Sheet ${sheet}: ${result.entries.length} symbols`);
  }
  if (symbols.length !== 296) throw new Error(`Expected 296 symbols; found ${symbols.length}`);
  if (new Set(symbols.map(s=>s.id)).size !== symbols.length) throw new Error('Duplicate asset IDs');
  save('manifest.json', JSON.stringify({ formatVersion: 1, attribution: 'Steven Baltakatei Sandoval (Baltakatei)', license: 'CC-BY-SA-4.0', source: 'https://commons.wikimedia.org/wiki/File:ISO_10628-2_2012_Symbols.pdf', method: 'Grouped source SVG geometry re-authored as standalone SVG with padded bounds; labels kept in metadata.', sheets, symbols }, null, 2)+'\n');
  const data = JSON.stringify({ sheets, symbols }).replaceAll('<','\\u003c');
  save('catalog.html', `<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>ISO 10628 symbol reference</title>
<style>body{font:15px system-ui;margin:32px;background:#f4f6f8;color:#1b2838}h1{margin-bottom:8px}p{max-width:950px;line-height:1.6}header{position:sticky;top:0;background:#f4f6f8;padding:12px 0;z-index:1}input,select{font:inherit;padding:10px;border:1px solid #aebdca;border-radius:5px}input{width:min(520px,70%)}#grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:16px}article{background:white;border:1px solid #d4dde5;border-radius:8px;padding:16px}img{width:100%;height:140px;object-fit:contain}b{display:block;margin:12px 0 5px}small{color:#56677a}a{color:#175a93}#count{margin:12px 0} .name{min-height:44px;line-height:1.4}</style>
<h1>Process symbol reference</h1><p>296 editable symbols from the seven ISO 10628-2:2012 reference sheets drawn by Steven Baltakatei Sandoval. Adapted into individual SVGs under <a href="https://creativecommons.org/licenses/by-sa/4.0/">CC BY-SA 4.0</a>. <a href="README.md">AI design guidance</a> · <a href="sources/original.pdf">Original PDF</a> · <a href="manifest.json">Manifest</a></p>
<header><input id="search" type="search" aria-label="Search symbols" placeholder="Search registration or description, e.g. pump"><select id="sheet" aria-label="Sheet"><option value="">All sheets</option></select><div id="count" aria-live="polite"></div></header><main id="grid"></main>
<script>const data=${data}; const search=document.querySelector('#search'),sheet=document.querySelector('#sheet'),grid=document.querySelector('#grid');for(const s of data.sheets){const o=document.createElement('option');o.value=s.sheet;o.textContent='Sheet '+s.sheet+' · '+s.count+' symbols';sheet.append(o)}function render(){grid.replaceChildren();const q=search.value.toLowerCase();const found=data.symbols.filter(s=>(!sheet.value||s.sheet==sheet.value)&&(s.registration+' '+s.description).toLowerCase().includes(q));document.querySelector('#count').textContent=found.length+' / 296 symbols';for(const s of found){const card=document.createElement('article'),a=document.createElement('a'),img=document.createElement('img'),b=document.createElement('b'),name=document.createElement('div'),small=document.createElement('small');a.href=s.file;img.src=s.file;img.alt=s.description;img.loading='lazy';a.append(img);b.textContent=s.registration;name.className='name';name.textContent=s.description;small.textContent='Sheet '+s.sheet+' · '+s.sourceGroup;card.append(a,b,name,small);grid.append(card)}}search.addEventListener('input',render);sheet.addEventListener('change',render);render();</script></html>`);
  console.log(`Generated ${symbols.length} SVG assets and catalog at ${fileURLToPath(root)}`);
} finally { await browser.close(); }
