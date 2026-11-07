import { createHash } from 'node:crypto';
import { execFileSync } from 'node:child_process';
import { mkdirSync, writeFileSync, readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
const dir = fileURLToPath(new URL('../../docs/design-reference/iso-10628/sources/', import.meta.url));
mkdirSync(dir, { recursive: true });
const files = [];
for (let i = 0; i <= 7; i++) {
  const name = i ? `ISO_10628-2_2012_Symbols_Sheet_${i}.svg` : 'ISO_10628-2_2012_Symbols.pdf';
  const hash = createHash('md5').update(name).digest('hex');
  const url = `https://upload.wikimedia.org/wikipedia/commons/${hash[0]}/${hash.slice(0, 2)}/${name}`;
  const local = i ? `sheet-${i}.svg` : 'original.pdf';
  if (!existsSync(`${dir}/${local}`)) {
    execFileSync('curl.exe', ['--ipv4', '--silent', '--show-error', '--location', '--fail', '--max-time', '60', url, '-o', `${dir}/${local}`]);
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  files.push({ file: local, url, sha256: createHash('sha256').update(readFileSync(`${dir}/${local}`)).digest('hex') });
  console.log(local);
}
writeFileSync(`${dir}/provenance.json`, JSON.stringify({ retrieved: '2026-11-07', author: 'Steven Baltakatei Sandoval', license: 'CC-BY-SA-4.0', source: 'https://commons.wikimedia.org/wiki/File:ISO_10628-2_2012_Symbols.pdf', files }, null, 2) + '\n');
const fontDir = `${dir}/font`;
mkdirSync(fontDir, { recursive: true });
const fontFiles = [
  ['osifont.ttf', 'https://raw.githubusercontent.com/hikikomori82/osifont/master/osifont.ttf'],
  ['MANUAL', 'https://raw.githubusercontent.com/hikikomori82/osifont/master/MANUAL'],
  ['README.md', 'https://raw.githubusercontent.com/hikikomori82/osifont/master/README.md'],
  ['GPL-3.0.txt', 'https://raw.githubusercontent.com/spdx/license-list-data/main/text/GPL-3.0-or-later.txt'],
];
const fontSources = [];
for (const [file, url] of fontFiles) {
  if (!existsSync(`${fontDir}/${file}`)) execFileSync('curl.exe', ['--silent','--show-error','--location','--fail','--max-time','60',url,'-o',`${fontDir}/${file}`]);
  fontSources.push({file,url,sha256:createHash('sha256').update(readFileSync(`${fontDir}/${file}`)).digest('hex')});
}
writeFileSync(`${fontDir}/provenance.json`,JSON.stringify({retrieved:'2026-11-07',source:'https://github.com/hikikomori82/osifont',license:'GPL-3.0 with upstream font exception (see MANUAL)',files:fontSources},null,2)+'\n');
