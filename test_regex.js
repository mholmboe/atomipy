const pdb = `CRYST1 1
MODEL 1
ATOM 1
ENDMDL
CRYST1 2
MODEL 2
ATOM 2
ENDMDL`;

const crystMatch = pdb.match(/^CRYST1.*$/m);
const crystStr = crystMatch[0];
console.log("Found crystStr:", crystStr);
const processed = pdb.replace(/^(MODEL\s+\d+)\r?\n(?!CRYST1)/gm, `$1\n${crystStr}\n`);
console.log(processed);
