import assert from 'node:assert/strict';
import { chromium } from 'playwright';
import { readFileSync, writeFileSync, mkdirSync, readdirSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
const root = new URL('../../docs/design-reference/iso-10628/', import.meta.url);
const manifest = JSON.parse(readFileSync(new URL('manifest.json', root)));
assert.equal(manifest.symbols.length, 296);
assert.equal(readdirSync(new URL('symbols/', root)).length, 296);
const provenance = JSON.parse(readFileSync(new URL('sources/provenance.json', root)));
for (const s of provenance.files) assert.equal(createHash('sha256').update(readFileSync(new URL(`sources/${s.file}`, root))).digest('hex'), s.sha256);
const fontProvenance = JSON.parse(readFileSync(new URL('sources/font/provenance.json', root)));
for (const s of fontProvenance.files) assert.equal(createHash('sha256').update(readFileSync(new URL(`sources/font/${s.file}`, root))).digest('hex'), s.sha256);
const browser = await chromium.launch({ channel: 'chrome', headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const results = [];
try {
  for (const {sheet} of manifest.sheets) {
    const entries = manifest.symbols.filter(s=>s.sheet===sheet);
    await page.goto(new URL(`sources/sheet-${sheet}.svg`,root).href);
    const sourceGroups = await page.evaluate(layers => layers.flatMap(id=>Array.from(document.getElementById(id).children).filter(e=>e.localName==='g').map(e=>e.id)), [...new Set(entries.map(s=>s.sourceLayer))]);
    assert.deepEqual([...sourceGroups].sort(),entries.map(s=>s.sourceGroup).sort(),`Unmapped source groups on sheet ${sheet}`);
  }
  for (const s of manifest.symbols) {
    await page.goto(new URL(s.file, root).href);
    await page.evaluate(() => document.fonts.ready);
    const result = await page.evaluate(async () => {
      const svg = document.documentElement;
      if (svg.localName !== 'svg' || document.querySelector('parsererror')) throw new Error('Invalid SVG');
      const ids = new Set(Array.from(svg.querySelectorAll('[id]')).map(e=>e.id));
      const refs = Array.from(new XMLSerializer().serializeToString(svg).matchAll(/(?:url\(#([^)]*)\)|(?:xlink:)?href="#([^"]*)")/g)).map(m=>m[1]||m[2]);
      const missingRefs = refs.filter(id=>!ids.has(id));
      const img = new Image();
      const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)], {type:'image/svg+xml'}));
      img.src = url; await img.decode();
      const c = document.createElementNS('http://www.w3.org/1999/xhtml', 'canvas');
      c.width=512;c.height=512;const ctx=c.getContext('2d');ctx.drawImage(img,0,0,512,512);
      const pixels=ctx.getImageData(0,0,512,512).data;
      let ink=0,edge=0;
      for(let y=0;y<512;y++)for(let x=0;x<512;x++){
        const i=(y*512+x)*4;
        if(pixels[i+3]>16 && Math.min(pixels[i],pixels[i+1],pixels[i+2])<200){ink++;if(x<2||y<2||x>509||y>509)edge++;}
      }
      URL.revokeObjectURL(url);
      return { missingRefs, ink, edge };
    });
    assert.deepEqual(result.missingRefs, [], `Unresolved SVG dependencies: ${s.id}`);
    assert.ok(result.ink > 0, `Blank rendering: ${s.id}`);
    assert.equal(result.edge, 0, `Clipped drawing: ${s.id}`);
    results.push({id:s.id,...result});
  }
  await page.goto(new URL('catalog.html', root).href);
  assert.equal(await page.locator('article').count(),296);
  await page.locator('#search').fill('pump');
  const pumpCount = await page.locator('article').count();
  assert.ok(pumpCount>0 && pumpCount<296);
  await page.locator('#search').fill('this-symbol-does-not-exist');
  assert.equal(await page.locator('article').count(),0);
  await page.locator('#search').fill('');
  mkdirSync(new URL('previews/',root),{recursive:true});
  for(const sheet of manifest.sheets){
    await page.locator('#sheet').selectOption(String(sheet.sheet));
    assert.equal(await page.locator('article').count(),sheet.count);
    // Trigger lazy images before taking a full-page screenshot.
    await page.locator('img').evaluateAll(imgs=>imgs.forEach(i=>i.loading='eager'));
    await page.locator('img').evaluateAll(imgs=>Promise.all(imgs.map(i=>i.decode())));
    await page.screenshot({path:fileURLToPath(new URL(`previews/sheet-${sheet.sheet}.png`,root)),fullPage:true});
  }
  writeFileSync(new URL('verification.json',root),JSON.stringify({checked:296,sourceHashes:'passed',xmlAndReferences:'passed',nonblankRendering:'passed',unclippedRendering:'passed',catalogSearch:'passed',sheetCounts:manifest.sheets.map(({sheet,count})=>({sheet,count})),pumpSearchCount:pumpCount,renderer:'Playwright / installed Google Chrome',results},null,2)+'\n');
  console.log('Verified 296 rendered SVGs, source hashes, dependency references, edge padding, search and all seven sheet filters.');
} finally {await browser.close();}
