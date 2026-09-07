const path = require('path');
const fs = require('fs');
const [,, jsPath, inputPath, ...args] = process.argv;
const factoryName = 'create_' + path.basename(jsPath, '.js').replace(/-/g, '_');
const mod = require(path.resolve(jsPath));
const factory = mod[factoryName] || mod.default || mod;
let out = '';
factory({ print: (s) => { out += s + '\n'; }, printErr: (s) => process.stderr.write(s + '\n') }).then((m) => {
  m.FS.writeFile('/input.fas', fs.readFileSync(inputPath));
  const code = m.callMain([...args, '/input.fas']);
  process.stdout.write(out);
  process.exit(code);
}).catch((e) => { console.error(e); process.exit(1); });
