import { createHash } from 'node:crypto';
import { execFileSync } from 'node:child_process';
import { mkdirSync, writeFileSync, readFileSync, existsSync, renameSync, rmSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

// Every source is checked against the SHA-256 recorded when the library was
// built (2026-11-07). A mismatch fails instead of recording whatever arrived.
// GitHub sources are pinned to commits; Wikimedia Commons has no immutable URL
// for a file's current revision, so its files rely on the pinned hashes alone.
const sha256 = path => createHash('sha256').update(readFileSync(path)).digest('hex');

function fetchVerified(url, path, expected, curlArgs) {
  if (!existsSync(path)) {
    const partial = `${path}.partial`;
    rmSync(partial, { force: true });
    execFileSync('curl.exe', [...curlArgs, '--silent', '--show-error', '--location', '--fail', '--max-time', '60', url, '-o', partial]);
    const actual = sha256(partial);
    if (actual !== expected) {
      rmSync(partial, { force: true });
      throw new Error(`${url}: SHA-256 ${actual} does not match pinned ${expected}; review upstream change explicitly`);
    }
    renameSync(partial, path);
    return true;
  }
  const actual = sha256(path);
  if (actual !== expected) throw new Error(`${path}: SHA-256 ${actual} does not match pinned ${expected}`);
  return false;
}

const dir = fileURLToPath(new URL('../../docs/design-reference/iso-10628/sources/', import.meta.url));
mkdirSync(dir, { recursive: true });
const sheetSha256 = [
  '047d4f3561bedf3132ab311533edbbb44de9d2a713a0e4fa5c00685c07b67e7d', // original.pdf
  'b0fbbde870ffc72f6e83255a89f35bcae2b8a04b7d6017a3e730fafb666fdaa0',
  '6b45e50237fbc91031889ff541b6cc78dabe0f4e676bc52aa8e4f31d07ff9d63',
  '6ebb6d69078f0d2673182a05354c65af408de1d57af219a3d067338e6486a7ee',
  '43474d6ecb760ea427572e5d4f64bb95c61728f78f01ea5bb017a2af2b445e6e',
  '420aca0b09b9dd3fe87ddb1fa7f9ae8e57e5059b901b9b780fd0dc7452b71744',
  '89c08ba5e644e548a0b4889194849fe60bf3af49a2a702a2932a64b4fcd1fb3b',
  'e36ebfeb57ae5205a4f161b7b30bb371451a944e648c075cb1be527664fae8a1',
];
const files = [];
for (let i = 0; i <= 7; i++) {
  const name = i ? `ISO_10628-2_2012_Symbols_Sheet_${i}.svg` : 'ISO_10628-2_2012_Symbols.pdf';
  const hash = createHash('md5').update(name).digest('hex');
  const url = `https://upload.wikimedia.org/wikipedia/commons/${hash[0]}/${hash.slice(0, 2)}/${name}`;
  const local = i ? `sheet-${i}.svg` : 'original.pdf';
  if (fetchVerified(url, `${dir}/${local}`, sheetSha256[i], ['--ipv4'])) {
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  files.push({ file: local, url, sha256: sheetSha256[i] });
  console.log(local);
}
writeFileSync(`${dir}/provenance.json`, JSON.stringify({ retrieved: '2026-11-07', author: 'Steven Baltakatei Sandoval', license: 'CC-BY-SA-4.0', source: 'https://commons.wikimedia.org/wiki/File:ISO_10628-2_2012_Symbols.pdf', files }, null, 2) + '\n');

const fontDir = `${dir}/font`;
mkdirSync(fontDir, { recursive: true });
// osifont master and spdx/license-list-data main as of the 2026-11-07 retrieval.
const OSIFONT_COMMIT = '2e9aa86a8a09b044e08c00f5a4a2505dc3ca9f6e';
const SPDX_COMMIT = '31846481de7ea17d92f2f46fad76abe58b2e9bb0';
const fontFiles = [
  ['osifont.ttf', `https://raw.githubusercontent.com/hikikomori82/osifont/${OSIFONT_COMMIT}/osifont.ttf`, '31e457a464b27ad0e3137bf957f0f4044ed9a3678df91eea2aea55c97c677208'],
  ['MANUAL', `https://raw.githubusercontent.com/hikikomori82/osifont/${OSIFONT_COMMIT}/MANUAL`, '3b22dd431554618d5cb31e48892111448dc07025a2e4a5ef1ee32f922d608b43'],
  ['README.md', `https://raw.githubusercontent.com/hikikomori82/osifont/${OSIFONT_COMMIT}/README.md`, '7ba462231fd835682f2594db1207ef2eed8dc4b0f25a160dfe8350313de526c2'],
  ['GPL-3.0.txt', `https://raw.githubusercontent.com/spdx/license-list-data/${SPDX_COMMIT}/text/GPL-3.0-or-later.txt`, 'fb981668c18a279e285fc4d83fba1e836cc84dd4daa73c9697d3cfd2d8aca6e0'],
];
const fontSources = [];
for (const [file, url, expected] of fontFiles) {
  fetchVerified(url, `${fontDir}/${file}`, expected, []);
  fontSources.push({ file, url, sha256: expected });
}
writeFileSync(`${fontDir}/provenance.json`, JSON.stringify({ retrieved: '2026-11-07', source: 'https://github.com/hikikomori82/osifont', license: 'GPL-3.0 with upstream font exception (see MANUAL)', files: fontSources }, null, 2) + '\n');
