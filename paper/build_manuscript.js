const fs = require('fs');
const D = require('docx');
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, ImageRun, AlignmentType,
  HeadingLevel, WidthType, BorderStyle, ShadingType, Footer, PageNumber, LevelFormat, PageBreak } = D;

const stats = JSON.parse(fs.readFileSync('data/stats.json'));
const mols = JSON.parse(fs.readFileSync('data/mols.json'));
const FONT = 'Times New Roman';
const W = 9026; // A4 content width (DXA) with 1" margins

// ---- citations: [#k,...] cites REFS keys (k = position in REFS, 1-based; ranges a–b allowed).
// Numbers are assigned in order of first citation and the reference list is printed in that order.
const citeOrder = [];
function cite(spec) {
  const keys = [];
  for (const part of spec.split(',')) {
    const [a, b] = part.split('–').map(Number);
    for (let k = a; k <= (b || a); k++) keys.push(k);
  }
  const nums = keys.map(k => { if (!citeOrder.includes(k)) citeOrder.push(k); return citeOrder.indexOf(k) + 1; }).sort((x, y) => x - y);
  const out = [];
  for (let i = 0; i < nums.length;) {
    let j = i;
    while (j + 1 < nums.length && nums[j + 1] === nums[j] + 1) j++;
    out.push(j - i >= 2 ? `${nums[i]}–${nums[j]}` : nums.slice(i, j + 1).join(','));
    i = j + 1;
  }
  return out.join(',');
}

// ---- inline markup: _{sub} ^{sup} **bold** *italic* [[flag]] (yellow highlight)
function runs(text, base = {}) {
  text = text.replace(/\\\*/g, '\u0001'); // \* is a literal asterisk, not italics
  const out = [];
  const re = /(\[#[^\]]+\]|\_\{[^}]*\}|\^\{[^}]*\}|\*\*[^*]+\*\*|\*[^*]+\*|\[\[[^\]]*\]\])/g;
  let last = 0, m;
  const push = (t, o = {}) => { if (t) out.push(new TextRun({ text: t.replace(/\u0001/g, '*'), font: FONT, size: 22, ...base, ...o })); };
  while ((m = re.exec(text))) {
    push(text.slice(last, m.index));
    const t = m[0];
    if (t.startsWith('[#')) push('[' + cite(t.slice(2, -1)) + ']', { superScript: true });
    else if (t.startsWith('_{')) push(t.slice(2, -1), { subScript: true });
    else if (t.startsWith('^{')) push(t.slice(2, -1), { superScript: true });
    else if (t.startsWith('**')) push(t.slice(2, -2), { bold: true });
    else if (t.startsWith('[[')) push('[' + t.slice(2, -2) + ']', { shading: { type: ShadingType.CLEAR, fill: 'FFFF00', color: 'auto' } });
    else push(t.slice(1, -1), { italics: true });
    last = m.index + t.length;
  }
  push(text.slice(last));
  return out;
}
const P = (text, o = {}) => new Paragraph({ children: runs(text, o.run), spacing: { after: 120, line: 300 }, alignment: o.align || AlignmentType.JUSTIFIED, ...o.para });
const H1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text: t, font: FONT, size: 26, bold: true })], spacing: { before: 280, after: 120 }, keepNext: true });
const H2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun({ text: t, font: FONT, size: 22, bold: true, italics: true })], spacing: { before: 200, after: 100 }, keepNext: true });
const CAP = (t) => new Paragraph({ children: runs(t, { size: 20 }), spacing: { before: 60, after: 200, line: 260 }, alignment: AlignmentType.JUSTIFIED });
const TCAP = (t) => new Paragraph({ children: runs(t, { size: 20 }), spacing: { before: 200, after: 80, line: 260 }, alignment: AlignmentType.JUSTIFIED, keepNext: true });

// ---- tables
const thin = { style: BorderStyle.SINGLE, size: 4, color: '000000' };
const none = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
function cell(text, w, o = {}) {
  const keep = o.keep !== false;
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    borders: { top: o.top ? thin : none, bottom: o.bottom ? thin : none, left: none, right: none },
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    shading: o.fill ? { type: ShadingType.CLEAR, fill: o.fill, color: 'auto' } : undefined,
    children: [new Paragraph({ keepNext: keep, alignment: o.align || AlignmentType.CENTER, children: runs(text, { size: 20, bold: !!o.bold }) })],
  });
}
function table(cols, header, rows) {
  const total = cols.reduce((a, b) => a + b, 0);
  const hdr = new TableRow({ tableHeader: true, children: header.map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: i < 2 ? AlignmentType.LEFT : AlignmentType.CENTER })) });
  const body = rows.map((r, ri) => new TableRow({ cantSplit: true, children: r.map((c, i) => cell(c, cols[i], { align: i < 2 ? AlignmentType.LEFT : AlignmentType.CENTER, bottom: ri === rows.length - 1, keep: ri < rows.length - 1 })) }));
  return new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] });
}
const f3 = (x) => x.toFixed(3);
function statRows(lab, level = 'pooled') {
  const rows = [];
  for (const [ds, name] of [['unstrained', 'Unstrained'], ['strained', 'Strained'], ['combined', 'Combined']]) {
    for (const [k, i] of [['ID', 0], ['HY', 1], ['HO', 2]]) {
      const s = stats[`${lab}|${ds}|${level}|${k}`];
      rows.push([i === 0 ? name : '', k, String(s.n), f3(s.rmse), f3(s.mae), f3(s.me), f3(s.std)]);
    }
  }
  return rows;
}
const T_HDR = ['Dataset', 'Scheme', 'N', 'RMSE', 'MAE', 'ME', 'STD'];
const T_COLS = [1700, 1300, 1000, 1256, 1256, 1256, 1258];

// ---- figures
function pngSize(f) { const b = fs.readFileSync('figures/' + f); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
const ALT = {
  'fig1_schemes.png': 'Isodesmic, hypohomodesmotic and homodesmotic reactions generated for butan-2-one, drawn as structures with stoichiometric coefficients.',
  'fig2_distributions.png': 'Violin and box plots of post-refinement errors for isodesmic, hypohomodesmotic and homodesmotic reactions, for unstrained and strained molecules.',
  'fig3_rmse.png': 'Bar chart of RMSE before and after single-point refinement for the three schemes, with bootstrap intervals, for unstrained, strained and all molecules.',
  'fig4_geometry.png': 'Parity plots of homodesmotic errors for two geometry levels before and after refinement, and cumulative distributions of the spread across geometries.',
  'fig5_strained.png': 'Dot plot of post-refinement errors for the 22 strained molecules, with structure drawings, for the three schemes.',
  'fig6_scalability.png': 'Number of distinct reference species and total reference molecules versus target size for polyenes, aromatic hydrocarbons and C60, and the homodesmotic reaction of C60.',
  'figS2_wildcard_libraries.png': 'The wildcard fragment libraries of the isogyric, isodesmic and hypohomodesmotic levels drawn as Kekule structures, with X for any atom.',
  'figS3_wildcard_library_homodesmotic.png': 'The wildcard fragment library of the homodesmotic level drawn as Kekule structures, with X for any atom and Y for any non-carbon atom.',
  'figS1_distributions_preSPR.png': 'Violin and box plots of pre-refinement errors for the three schemes, for unstrained and strained molecules.',
};
function fig(file, wpx) {
  const [w, h] = pngSize(file);
  return new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, spacing: { before: 120 },
    children: [new ImageRun({ type: 'png', data: fs.readFileSync('figures/' + file), transformation: { width: wpx, height: Math.round(wpx * h / w) },
      altText: { title: file, description: ALT[file] || file, name: file } })] });
}

// ---- references (single list, numbered in order of first citation)
const REFS = [
 'Wu, H.; Doner, A. C.; Pang, H. W.; Green, W. H. Toward Accurate Quantum Mechanical Thermochemistry: (2) Optimal Methods for Enthalpy Calculations from Comprehensive Benchmarks of 284 Model Chemistries. *J. Phys. Chem. A* **2025**, *129* (17), 3949–3962.',
 'Ruscic, B.; Bross, D. H. Thermochemistry. In *Computer Aided Chemical Engineering*, Vol. 45; Faravelli, T., Manenti, F., Ranzi, E., Eds.; Elsevier, 2019; pp 3–114.',
 'Elm, J.; Ayoubi, D.; Engsvang, M.; Jensen, A. B.; Knattrup, Y.; Kubečka, J.; Bready, C. J.; Fowler, V. R.; Harold, S. E.; Longsworth, O. M.; et al. Quantum Chemical Modeling of Organic Enhanced Atmospheric Nucleation: A Critical Review. *WIREs Comput. Mol. Sci.* **2023**, *13* (5), e1662.',
 'Khan, S. S.; Yu, X.; Wade, J. R.; Malmgren, R. D.; Broadbelt, L. J. Thermochemistry of Radicals and Molecules Relevant to Atmospheric Chemistry: Determination of Group Additivity Values Using G3//B3LYP Theory. *J. Phys. Chem. A* **2009**, *113* (17), 5176–5194.',
 'Curtarolo, S.; Hart, G. L.; Nardelli, M. B.; Mingo, N.; Sanvito, S.; Levy, O. The High-Throughput Highway to Computational Materials Design. *Nat. Mater.* **2013**, *12* (3), 191–201.',
 'Friedrich, R.; Usanmaz, D.; Oses, C.; Supka, A.; Fornari, M.; Buongiorno Nardelli, M.; Toher, C.; Curtarolo, S. Coordination Corrected Ab Initio Formation Enthalpies. *npj Comput. Mater.* **2019**, *5*, 59.',
 'Curtiss, L. A.; Raghavachari, K.; Trucks, G. W.; Pople, J. A. Gaussian-2 Theory for Molecular Energies of First- and Second-Row Compounds. *J. Chem. Phys.* **1991**, *94* (11), 7221–7230.',
 'Pople, J. A.; Head-Gordon, M.; Fox, D. J.; Raghavachari, K.; Curtiss, L. A. Gaussian-1 Theory: A General Procedure for Prediction of Molecular Energies. *J. Chem. Phys.* **1989**, *90* (10), 5622–5629.',
 'Petersson, G. A.; Tensfeldt, T. G.; Montgomery, J. A. A Complete Basis Set Model Chemistry. III. The Complete Basis Set-Quadratic Configuration Interaction Family of Methods. *J. Chem. Phys.* **1991**, *94* (9), 6091–6101.',
 'Martin, J. M. L.; de Oliveira, G. Towards Standard Methods for Benchmark Quality Ab Initio Thermochemistry—W1 and W2 Theory. *J. Chem. Phys.* **1999**, *111* (5), 1843–1856.',
 'Karton, A.; Rabinovich, E.; Martin, J. M.; Ruscic, B. W4 Theory for Computational Thermochemistry: In Pursuit of Confident Sub-kJ/mol Predictions. *J. Chem. Phys.* **2006**, *125* (14), 144108.',
 'Curtiss, L. A.; Raghavachari, K.; Redfern, P. C.; Rassolov, V.; Pople, J. A. Gaussian-3 (G3) Theory for Molecules Containing First and Second-Row Atoms. *J. Chem. Phys.* **1998**, *109* (18), 7764–7776.',
 'Curtiss, L. A.; Redfern, P. C.; Raghavachari, K. Gaussian-4 Theory. *J. Chem. Phys.* **2007**, *126* (8), 084108.',
 'Boese, A. D.; Oren, M.; Atasoylu, O.; Martin, J. M.; Kállay, M.; Gauss, J. W3 Theory: Robust Computational Thermochemistry in the kJ/mol Accuracy Range. *J. Chem. Phys.* **2004**, *120* (9), 4129–4141.',
 'Tajti, A.; Szalay, P. G.; Császár, A. G.; Kállay, M.; Gauss, J.; Valeev, E. F.; Flowers, B. A.; Vázquez, J.; Stanton, J. F. HEAT: High Accuracy Extrapolated Ab Initio Thermochemistry. *J. Chem. Phys.* **2004**, *121* (23), 11599–11613.',
 'Wheeler, S. E.; Houk, K. N.; Schleyer, P. v. R.; Allen, W. D. A Hierarchy of Homodesmotic Reactions for Thermochemistry. *J. Am. Chem. Soc.* **2009**, *131* (7), 2547–2560.',
 'Hehre, W. J.; Ditchfield, R.; Radom, L.; Pople, J. A. Molecular Orbital Theory of the Electronic Structure of Organic Compounds. V. Molecular Theory of Bond Separation. *J. Am. Chem. Soc.* **1970**, *92* (16), 4796–4801.',
 'Hess, B. A., Jr.; Schaad, L. J. Ab Initio Calculation of Resonance Energies. Benzene and Cyclobutadiene. *J. Am. Chem. Soc.* **1983**, *105* (26), 7500–7505.',
 'Nyulászi, L.; Várnai, P.; Veszprémi, T. About the Aromaticity of Five-Membered Heterocycles. *J. Mol. Struct. THEOCHEM* **1995**, *358*, 55–61.',
 'Vianello, R.; Liebman, J. F.; Maksić, Z. B. In Search of Ultrastrong Brønsted Neutral Organic Superacids: A DFT Study on Some Cyclopentadiene Derivatives. *Chem. Eur. J.* **2004**, *10* (22), 5751–5760.',
 'Barić, D.; Kovačević, B.; Maksić, Z. B.; Müller, T. A Novel Approach in Analyzing Aromaticity by Homo- and Isostructural Reactions: An Ab Initio Study of Fluorobenzenes. *J. Phys. Chem. A* **2005**, *109* (46), 10594–10606.',
 'Chesnut, D. B.; Davis, K. M. Resonance Revisited: A Consideration of the Calculation of Cyclic Conjugation Energies. *J. Comput. Chem.* **1997**, *18* (4), 584–593.',
 'El-Nahas, A. M.; Bozzelli, J. W.; Simmie, J. M.; Navarro, M. V.; Black, G.; Curran, H. J. Thermochemistry of Acetonyl and Related Radicals. *J. Phys. Chem. A* **2006**, *110*, 13618–13623.',
 'George, P.; Bock, C. W.; Trachtman, M. The Matching of Structural Elements in Reactions for Evaluating Stabilization Energies for Benzene and Monosilabenzene. *Theor. Chim. Acta* **1987**, *71*, 289–298.',
 'Gimarc, B. M.; Zhao, M. Oxygen Ring Strain Energies Revisited—Effects of Terminal Atoms of the Chain Reference Structure. *J. Phys. Chem.* **1994**, *98* (6), 1596–1600.',
 'Zhao, M.; Gimarc, B. M. Strain Energies in Cyclic Oxygen, *n* = 3–8. *J. Phys. Chem.* **1993**, *97* (16), 4023–4030.',
 'Warren, D. S.; Gimarc, B. M. Strain Energies in Sulfur Rings, S_{n}, *n* = 3–8. *J. Phys. Chem.* **1993**, *97* (16), 4031–4035.',
 'Ramabhadran, R. O.; Raghavachari, K. Connectivity-Based Hierarchy for Theoretical Thermochemistry: Assessment Using Wave Function-Based Methods. *J. Phys. Chem. A* **2012**, *116* (28), 7531–7537.',
 'Ramabhadran, R. O.; Raghavachari, K. Theoretical Thermochemistry for Organic Molecules: Development of the Generalized Connectivity-Based Hierarchy. *J. Chem. Theory Comput.* **2011**, *7* (7), 2094–2103.',
 'Raghavachari, K.; Maier, S.; Collins, E. M.; Debnath, S.; Sengupta, A. Approaching Coupled Cluster Accuracy with Density Functional Theory Using the Generalized Connectivity-Based Hierarchy. *J. Chem. Theory Comput.* **2023**, *19* (13), 3763–3778.',
 'Wheeler, S. E. Homodesmotic Reactions for Thermochemistry. *WIREs Comput. Mol. Sci.* **2012**, *2* (2), 204–220.',
 'Ruscic, B.; Pinzón, R. E.; von Laszewski, G.; Kodeboyina, D.; Burcat, A.; Leahy, D.; Montoya, D.; Wagner, A. F. Active Thermochemical Tables: Thermochemistry for the 21st Century. *J. Phys. Conf. Ser.* **2005**, *16*, 561–570.',
 'Wu, H.; Payne, A. M.; Pang, H. W.; Menon, A.; Grambow, C. A.; Ranasinghe, D. S.; Dong, X.; Grinberg Dana, A.; Green, W. H. Toward Accurate Quantum Mechanical Thermochemistry: (1) Extensible Implementation and Comparison of Bond Additivity Corrections and Isodesmic Reactions. *J. Phys. Chem. A* **2024**, *128* (21), 4335–4352.',
 'Sen, S. K.; Agarwal, H.; Sen, S. Chemical Equation Balancing: An Integer Programming Approach. *Math. Comput. Model.* **2006**, *44* (7), 678–691.',
 'Ingman, V. M.; Schaefer, A. J.; Andreola, L. R.; Wheeler, S. E. QChASM: Quantum Chemistry Automation and Structure Manipulation. *WIREs Comput. Mol. Sci.* **2021**, *11* (4), e1510.',
 'Landrum, G.; et al. RDKit: Open-Source Cheminformatics, Release 2024.03.3. https://www.rdkit.org (accessed October 27, 2026).',
 'Mitchell, S.; Peschiera, F.; Roy, J. S. PuLP: A Linear Programming Modeler, Version 2.8.0. https://github.com/coin-or/pulp (accessed October 27, 2026).',
 'Forrest, J.; Vigerske, S.; Ralphs, T.; Santos, H. G.; Hafer, L.; Kristjansson, B.; et al. COIN-OR CBC: Branch-and-Cut Solver, Version 2.10.3 (bundled with PuLP 2.8.0). https://github.com/coin-or/Cbc (accessed October 27, 2026).',
 'Becke, A. D. Density-Functional Thermochemistry. III. The Role of Exact Exchange. *J. Chem. Phys.* **1993**, *98* (7), 5648–5652.',
 'Lee, C.; Yang, W.; Parr, R. G. Development of the Colle-Salvetti Correlation-Energy Formula into a Functional of the Electron Density. *Phys. Rev. B* **1988**, *37* (2), 785–789.',
 'Zhao, Y.; Truhlar, D. G. The M06 Suite of Density Functionals for Main Group Thermochemistry, Thermochemical Kinetics, Noncovalent Interactions, Excited States, and Transition Elements. *Theor. Chem. Acc.* **2008**, *120*, 215–241.',
 'Chai, J.-D.; Head-Gordon, M. Long-Range Corrected Hybrid Density Functionals with Damped Atom–Atom Dispersion Corrections. *Phys. Chem. Chem. Phys.* **2008**, *10* (44), 6615–6620.',
 'Mardirossian, N.; Head-Gordon, M. ωB97X-V: A 10-Parameter, Range-Separated Hybrid, Generalized Gradient Approximation Density Functional with Nonlocal Correlation, Designed by a Survival-of-the-Fittest Strategy. *Phys. Chem. Chem. Phys.* **2014**, *16* (21), 9904–9924.',
 'Weigend, F.; Ahlrichs, R. Balanced Basis Sets of Split Valence, Triple Zeta Valence and Quadruple Zeta Valence Quality for H to Rn: Design and Assessment of Accuracy. *Phys. Chem. Chem. Phys.* **2005**, *7* (18), 3297–3305.',
 'Frisch, M. J.; Trucks, G. W.; Schlegel, H. B.; Scuseria, G. E.; Robb, M. A.; Cheeseman, J. R.; Scalmani, G.; Barone, V.; Petersson, G. A.; Nakatsuji, H.; et al. *Gaussian 16*, Revision C.01; Gaussian, Inc.: Wallingford, CT, 2019.',
 'Smith, D. G. A.; Burns, L. A.; Simmonett, A. C.; Parrish, R. M.; Schieber, M. C.; Galvelis, R.; Kraus, P.; Kruse, H.; Di Remigio, R.; Alenaizan, A.; et al. Psi4 1.4: Open-Source Software for High-Throughput Quantum Chemistry. *J. Chem. Phys.* **2020**, *152* (18), 184108.',
  'Wodrich, M. D.; Corminboeuf, C.; Schleyer, P. v. R. Systematic Errors in Computed Alkane Energies Using B3LYP and Other Popular DFT Functionals. *Org. Lett.* **2006**, *8* (17), 3631–3634.',
  'Ruscic, B.; Pinzon, R. E.; Morton, M. L.; von Laszewski, G.; Bittner, S. J.; Nijsure, S. G.; Amin, K. A.; Minkoff, M.; Wagner, A. F. Introduction to Active Thermochemical Tables: Several "Key" Enthalpies of Formation Revisited. *J. Phys. Chem. A* **2004**, *108* (45), 9979–9997.',
  'Ruscic, B.; Bross, D. H. Active Thermochemical Tables (ATcT) Thermochemical Values, ver. 1.156; Argonne National Laboratory: Lemont, IL. https://atct.anl.gov (OSTI ID 2440253). [[AUTHOR: add the DOI of ver. 1.156.]]',
  'Riniker, S.; Landrum, G. A. Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation. *J. Chem. Inf. Model.* **2015**, *55* (12), 2562–2574.',
  'Halgren, T. A. MMFF VI. MMFF94s Option for Energy Minimization Studies. *J. Comput. Chem.* **1999**, *20* (7), 720–729.',
  'Hehre, W. J.; Ditchfield, R.; Pople, J. A. Self-Consistent Molecular Orbital Methods. XII. Further Extensions of Gaussian-Type Basis Sets for Use in Molecular Orbital Studies of Organic Molecules. *J. Chem. Phys.* **1972**, *56* (5), 2257–2261.',
  'Hariharan, P. C.; Pople, J. A. The Influence of Polarization Functions on Molecular Orbital Hydrogenation Energies. *Theor. Chim. Acta* **1973**, *28* (3), 213–222.',
  'Clark, T.; Chandrasekhar, J.; Spitznagel, G. W.; Schleyer, P. v. R. Efficient Diffuse Function-Augmented Basis Sets for Anion Calculations. III. The 3-21+G Basis Set for First-Row Elements, Li–F. *J. Comput. Chem.* **1983**, *4* (3), 294–301.',
];
const refPara = (t) => new Paragraph({ numbering: { reference: 'refs', level: 0 }, spacing: { after: 60, line: 252 }, alignment: AlignmentType.LEFT, children: runs(t, { size: 20 }) });

// ---- CONTENT
const c = [];
c.push(new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 160 }, children: [new TextRun({ text: 'Automated Computational Thermochemistry Using Homodesmotic Reactions', font: FONT, size: 34, bold: true })] }));
c.push(new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 60 }, children: [...runs('Aaron D. Seymour and Steven E. Wheeler^{*}', { size: 24 })] }));
c.push(new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 60 }, children: runs('Department of Chemistry, University of Georgia, Athens, GA 30602, USA', { size: 22, italics: true }) }));
c.push(new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 200 }, children: runs('E-mail: swheele2@uga.edu', { size: 22 }) }));

c.push(H1('Abstract'));
c.push(P('Accurate enthalpies of formation are needed throughout chemistry and engineering, but experimental values are limited and computational approaches are often expensive or not accurate enough. We present an open-source framework that formulates isodesmic, hypohomodesmotic and homodesmotic bond separation reactions as integer linear programs, generates them automatically for organic molecules, and combines them with quantum-chemical energies and Active Thermochemical Tables (ATcT) reference values. Homodesmotic reactions preserve more of the local bonding environment of the target than the lower rungs of the hierarchy, and the framework removes the manual effort of constructing them. On 77 ATcT species containing C, H, N, O, F and Cl (55 acyclic, 22 with three- to five-membered rings), homodesmotic reactions give an overall root-mean-square error of 1.10 kcal mol^{−1} after single-point refinement (0.35 for the acyclic species), compared with 1.21 for hypohomodesmotic and 1.83 for isodesmic reactions. The improvement over isodesmic reactions holds molecule by molecule; that over hypohomodesmotic reactions is comparable to the uncertainty of the reference data. All schemes degrade for strained rings, most of all for polycyclic cages, whose ring strain no scheme conserves. The number of reference species depends only on local bonding motifs and does not grow with molecular size.'));
c.push(P('**Keywords:** computational chemistry · density functional calculations · homodesmotic reactions · integer linear programming · thermochemistry', { align: AlignmentType.LEFT }));

// 1 Introduction
c.push(H1('1. Introduction'));
c.push(P('Accurate computational enthalpies of formation are important for fields ranging from combustion modeling[#1,2] and atmospheric chemistry modeling[#3,4] to materials design.[#5,6] Small errors in reference enthalpies can propagate into large uncertainties in predictions, underscoring the need for reliable thermochemical data. Scores of composite approaches based on rigorous *ab initio* methods have been developed, many of which can deliver highly accurate enthalpies of formation.[#7–15] These methods rely extensively on high-level coupled-cluster theory and computations with large basis sets. The need for these high levels of theory is due to the use of atomization energies as reference reactions, which exhibit slow convergence with respect to both electron correlation effects and basis set. The result is that the more rigorous of these methods are only applicable to small molecules.'));
c.push(P('An alternative approach to computational thermochemistry is to use error-balanced reactions.[#16] This concept dates back to 1970, when Hehre, Ditchfield, Radom and Pople introduced isodesmic reactions in which the number of each bond type (C–C, C=C, C–O, etc.) is equal in reactants and products.[#17] This balancing of bond types often leads to a systematic cancellation of errors in computed reaction enthalpies, enabling relatively accurate predictions even with the modest levels of theory that were practical at the time. The following decades witnessed a Cambrian explosion of reaction types (homodesmotic, quasihomodesmotic, hyperhomodesmotic, homomolecular homodesmotic, isogeitonic, homoplesiotic, etc.).[#18–27] These reaction schemes are based on the balancing of various factors including bond types, hybridization and connectivity. However, two non-equivalent definitions of "homodesmotic reaction" were extant in the literature, leading to considerable confusion. In 2009, Wheeler and co-workers[#16] refined some of these definitions to provide a systematic hierarchy of homodesmotic reactions that offer increasingly balanced errors in thermochemical computations: isogyric < isodesmic < hypohomodesmotic < homodesmotic < hyperhomodesmotic. Wheeler *et al.* also extended the concept of isodesmic bond separation reactions to the other rungs on the homodesmotic hierarchy. In this way, for most target molecules there is a unique decomposition reaction corresponding to each rung on the hierarchy.'));
c.push(P('Moreover, Wheeler *et al.* noted that if one excludes homodesmotic reactions, these reaction types can be classified as *n*-homodesmotic reactions, in which fragments of *n* atoms are preserved in reactants and products (0-homodesmotic is isogyric; 1-homodesmotic is isodesmic; 2-homodesmotic is hypohomodesmotic; and 3-homodesmotic is hyperhomodesmotic).[#16] Raghavachari and co-workers[#28–30] have employed a similar strategy with their connectivity-based hierarchy, where the CBH*n* of Raghavachari *et al.* corresponds to *n*-homodesmotic of Wheeler *et al.* Homodesmotic reactions, using the refined definition of Wheeler *et al.*, are effectively CBH-2 reactions (preserving three-heavy-atom fragments) plus conjugated fragments, and thus tend to fall somewhere between CBH-2 and CBH-3.'));
c.push(P('Critically, enthalpies of higher-level error-balanced reaction classes generally provide considerably faster convergence with respect to electron correlation and basis set effects.[#16] This means that, in many cases, the judicious choice of reference reaction can provide converged, accurate enthalpies of formation even at relatively low levels of theory. Similarly, error-balanced reactions often obviate the need for the sundry corrections often included in high-accuracy computational thermochemical studies (anharmonic corrections to the zero-point vibrational energy as well as corrections for relativistic effects, core–valence correlation, the Born–Oppenheimer approximation, etc.).[#16,31]'));
c.push(P('A traditional drawback of these schemes was that they require accurate enthalpies of formation for relatively large reference species, which were not widely available in 2009. However, the Active Thermochemical Tables (ATcT) of Ruscic *et al.*[#48,32] now contain accurate enthalpies of formation (including rigorous uncertainties) for a large and growing set of relatively large molecules. A second, practical shortcoming is that some of these reaction schemes, particularly homodesmotic reactions, are difficult to generate in an automated way. Indeed, the literature is replete with putative homodesmotic reactions that do not satisfy either definition of a homodesmotic reaction. Critically, it is these higher-order reaction schemes that provide the most rigorous error cancellation by accounting for subtle electronic and structural features.'));
c.push(P('While the CBH hierarchy of Raghavachari *et al.*[#28–30] can be readily automated, it does not include homodesmotic reactions, instead skipping from hypohomodesmotic (CBH-2) to hyperhomodesmotic (CBH-3) reactions. As Wheeler *et al.*[#16] have previously shown, homodesmotic reactions can provide greater error cancellation than hypohomodesmotic reactions (CBH-2), without the added complexity of the larger reference molecules required for hyperhomodesmotic reactions (CBH-3). Green and co-workers[#33] recently introduced an open-source package to automatically evaluate enthalpies of formation using isodesmic reactions and bond additivity corrections; the present work targets the next rung up, automating homodesmotic reactions for general organic molecules.'));
c.push(P('In this work we present a new framework for generating rigorous homodesmotic reactions based on ILP. Integer programming has been used in a chemistry context by Sen *et al.*,[#34] and our work expands upon that idea with an application to error-balanced reactions. By encoding bond hybridization types and hydrogen atom environments as explicit constraint sets, our approach provides a general and automated generator of homodesmotic (and isodesmic and hypohomodesmotic) reactions for molecules including but not limited to hydrocarbons. When coupled with quantum chemical calculations and thermochemical reference data from the ATcT, this framework enables accurate and efficient determination of enthalpies of formation. We show how this, combined with the automation tools provided by QChASM,[#35] enables automated high-throughput predictions of enthalpies of formation using homodesmotic reactions. We benchmark it on 77 ATcT molecules and report where the homodesmotic scheme helps and where it does not. We also use C_{60} as an extreme example to show that the reference species of a homodesmotic reaction do not grow with molecular size, so that only the stoichiometric coefficients change as the target becomes larger.'));

// 2 Methods
c.push(H1('2. Computational Methods'));
c.push(H2('2.1 Reaction generation'));
c.push(P('A homodesmotic reaction conserves bond hybridization types and hydrogen environments simultaneously. A bond type is defined by the hybridization of the two bonded atoms and the bond order; a hydrogen environment is defined by the hybridization of a heavy atom and the number of hydrogens attached to it. Each scheme in the hierarchy is defined by the set of quantities it conserves (Table 1), and each higher scheme conserves more of the local environment of the target. Figure 1 shows the three reactions generated for butan-2-one.'));
c.push(TCAP('**Table 1.** Quantities conserved by each scheme and the fragments used to build its reactions (Figures S2 and S3, Table S5). Atoms are counted by element, hybridization (sp, sp^{2}, sp^{3}) and number of attached hydrogens (nH); wildcard atoms \\* match any element.'));
{
  const cols = [2000, 3500, 3526];
  const rows = [["Isogyric", "atoms by element", "single atoms; H_{2}"], ["Isodesmic", "atoms by element (including H); bonds by element pair and bond order", "X–Y, X=Y and X≡Y fragments; element hydrides; H_{2}"], ["Hypohomodesmotic", "heavy atoms by element and hybridization; heavy atoms by element and nH", "three- to five-atom fragments (7 patterns); two-atom bonds (3 patterns)"], ["Homodesmotic", "bonds by element pair, hybridization of both atoms and bond order; heavy atoms by element, hybridization and nH", "three- to five-atom fragments (17 patterns); smaller bonded and conjugated fragments (6 patterns)"]];
  const hdr = new TableRow({ tableHeader: true, children: ['Scheme', 'Conserved quantities', 'Fragments'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: AlignmentType.LEFT })) });
  const body = rows.map((r, ri) => new TableRow({ cantSplit: true, children: r.map((t, i) => cell(t, cols[i], { align: AlignmentType.LEFT, bottom: ri === rows.length - 1, keep: ri < rows.length - 1 })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
  c.push(P('', { para: { spacing: { after: 60 } } }));
}
c.push(P('Balancing these hybridization-specific quantities by hand is tedious and error-prone, particularly for homodesmotic reactions. We therefore formulate the problem as an integer linear program (ILP) that balances reactions by enforcing the conservation of atom counts, bond hybridization types and hydrogen environments. While similar in structure to previous ILP-based formulations of chemical reaction balancing,[#34] our implementation focuses on error-balanced homodesmotic relationships, in which structural and energetic errors are minimized by matching atom and bond environments across reactants and products.'));
c.push(fig('fig1_schemes.png', 600));
c.push(CAP('**Figure 1.** Reactions generated for butan-2-one (dashed box) at the isodesmic, hypohomodesmotic and homodesmotic levels. Each rung conserves more of the local environment of the target; in the homodesmotic reaction the carbonyl group is transferred between an acetaldehyde and a propanal unit while hybridization and hydrogen environments are preserved.'));
c.push(P('Reaction generation has three steps. First, *candidate species* are constructed. The input molecule M is kekulized (aromatic rings become alternating single and double bonds, so aromatic stabilization is not conserved; Section 3.7), hydrogens are added, and RDKit[#36] perceives the hybridization and hydrogen count of every atom. The wildcard fragments of the scheme (Figures S2 and S3, Table S5) are then matched onto M by substructure search; a wildcard atom matches any element. Each match becomes a real molecule by assigning the matched elements to the wildcards and completing the valences with hydrogen, and duplicates are removed using canonical SMILES. This gives a set of candidate products P_{1}…P_{n}, the fragments that can appear on the product side. Matching the smaller reactant fragments of the library onto these products gives a set of candidate co-reactants R_{1}…R_{m}, and M itself is added to the reactant side. The isodesmic candidates are built directly rather than mapped from a library: one product for each bond type present in M (X–Y, X=Y or X≡Y, saturated with hydrogen), the hydride of each element present, and H_{2}.'));
c.push(P('Second, the candidates are *balanced*. Let x_{i} and y_{j} be the non-negative integer stoichiometric coefficients of P_{i} and R_{j}, and let n_{d}(S) be the number of times descriptor d occurs in species S (for example, the number of C(sp^{2})=O(sp^{2}) bonds, or of sp^{3} carbons bearing three hydrogens). The reaction Σ_{j} y_{j}R_{j} → Σ_{i} x_{i}P_{i} is balanced when'));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 60, after: 120 }, children: runs('Σ_{i} x_{i} n_{d}(P_{i}) = Σ_{j} y_{j} n_{d}(R_{j})  for every descriptor d in the set D of the scheme, and y_{M} = 1.        (1)') }));
c.push(P('Each descriptor in Table 1 therefore contributes one linear equality, and the atom and hydrogen counts included in D guarantee mass balance. Fixing y_{M} = 1 makes the solution a decomposition of exactly one molecule of the target. For butan-2-one the homodesmotic reaction (butan-2-one + acetaldehyde → acetone + propanal) has, on each side, one C(sp^{3})–C(sp^{3}) single bond, three C(sp^{3})–C(sp^{2}) single bonds and two C(sp^{2})=O double bonds, and three CH_{3} groups, one CH_{2} group, one sp^{2} carbon without and one with a hydrogen, and two sp^{2} oxygens. The isodesmic reaction (butan-2-one + 3 CH_{4} → 3 C_{2}H_{6} + CH_{2}O) and the hypohomodesmotic reaction (butan-2-one + C_{2}H_{6} → propane + acetone) conserve correspondingly fewer quantities.'));
c.push(P('Third, the integer linear program is *solved* with the PuLP Python library (version 2.8.0)[#37] and its CBC solver (version 2.10.3).[#38] The identity reaction M → M satisfies Eq. 1 for every scheme and is uninformative, so its product variable is set to zero; if no other solution exists, the generator reports that no balanced reaction exists at that level, as happens for molecules that are their own smallest fragment. Among the feasible reactions the current release returns the one with the fewest species, minimizing Σx_{i} + Σy_{j}. The benchmark reactions were generated without this objective; re-solving all 77 benchmark molecules with it gave identical isodesmic and homodesmotic reactions for all 77 and identical hypohomodesmotic reactions for 73, and the four exceptions, for which more than one balanced hypohomodesmotic reaction exists (the benchmark uses the solution returned without an objective), are listed in Table S4. Without an objective the solver does not always return the same reaction when several exist: for glycine and cyclopropenone (homodesmotic) and HClO_{3} (hypohomodesmotic and homodesmotic, where the runs at different levels returned different reactions), the benchmark reactions differ from those returned by the current code, and Table S6 lists the reactions actually used. Species can also be required on either side of the reaction, or a substructure of M replaced before candidates are generated. Because Eq. 1 is linear and exact, the result is not a heuristic approximation: every returned reaction satisfies the conservation laws of its scheme, and the reactions for the molecules in Figure 6 were each generated in under 0.1 s.'));
c.push(H2('2.2 Quantum-chemical calculations'));
c.push(P('For each molecule in the resulting balanced reactions, conformer generation is carried out as a function of the number of rotatable bonds: n = max(5, 5 × n_{rot}) conformers were embedded with RDKit (version 2024.03.3; EmbedMultipleConfs, random seed 42), pruned by best-fit all-atom RMSD (0.15 Å), minimized with MMFF94s[#51] and pruned again.[#36] EmbedMultipleConfs uses the ETKDG method.[#50] All conformers were converted to input files with AaronTools/QChASM[#35] and optimized, with harmonic frequency calculations, using Gaussian 16 (Revision C.01)[#45] at B3LYP/6-31+G(d),[#39,40,52–54] M06-2X/6-31+G(d)[#41] and ωB97X-D/def2-TZVP;[#42,44] the lowest-energy conformer of each species is selected for further thermochemical analysis (the workflow also supports Psi4 and ORCA). All species were treated as closed-shell singlets, and the optimizations and frequency calculations used the Gaussian 16 default convergence criteria and integration grid (route Opt Freq, with no grid or convergence keywords). Six species first converged to structures with imaginary frequencies (2-butyne, acetone, butan-2-one and cyclopentane with small methyl-rotor or ring-pucker modes, and 1,2,3-trioxolane and ethynediol at saddle points); these were displaced along the imaginary modes and reoptimized with Opt=(CalcFC,Tight) and NoSymm, and every optimized structure used here is a minimum (Table S7). The lowest-energy conformer of each species, judged by E + ZPVE at the optimization level, is used in every reaction in which that species appears. Reference species that were computed in more than one run (in the original benchmark runs and again in the later runs for the new species) use their lowest-energy copy throughout. The total energy is the sum of the electronic energy and the unscaled zero-point vibrational energy (ZPVE).'));
c.push(P('Alternatively, a composite approach (SPR) can be employed, in which a single-point energy calculation is performed at a higher level of theory on the geometry optimized at the initial level; the single-point energy and ZPVE are then combined to yield a refined total energy. Here an ωB97X-V/def2-TZVP single-point energy[#43,44] was computed with Psi4 (version 1.9.1)[#46] on the lowest-energy optimized geometry of each species, using the Psi4 defaults for the SCF (density fitting) and the integration grid, and combined with the ZPVE from the lower-level frequency calculation. [[AUTHOR: confirm the Psi4 input contained no options beyond method and basis set.]]'));
c.push(H2('2.3 Enthalpies of formation and error statistics'));
c.push(P('The reaction enthalpy ΔH_{rxn} = Σ_{products} ν_{i}E_{i} − Σ_{reactants} ν_{j}E_{j} is converted from hartree (1 E_{h} = 2625.4996 kJ mol^{−1}) and combined with ATcT 0 K enthalpies of formation[#32] for every species other than the input molecule M:'));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 60, after: 120 }, children: runs('Δ_{f}H°_{0}(M) = [ Σ_{products} ν_{i}Δ_{f}H°_{0}(i) − Σ_{reactants, j≠M} ν_{j}Δ_{f}H°_{0}(j) − ΔH_{rxn} ] / ν_{M}') }));
c.push(P('Because the energies include ZPVE and no thermal corrections, the 0 K ATcT values (version 1.156)[#49] are the consistent reference; the 626 values used are distributed with the software (hfrpkg/data/ATcT_lib.txt) and, for the benchmark molecules, listed in Table S2. Errors are defined as computed minus ATcT values and summarized by the root-mean-square error (RMSE), mean absolute error (MAE), mean signed error (ME) and sample standard deviation (STD), pooled over the three levels of theory (*N* counts reactions, molecules × 3). After refinement the three geometries give near-replicate reactions, so inference is made at the level of molecules, using bootstrap resampling of molecules and an exact sign test on per-molecule mean absolute errors (Supporting Information, Note S1). The ATcT uncertainties (95% confidence intervals; Table S2) of the reference species were propagated in quadrature, ignoring covariance.'));

// 3 Results
c.push(H1('3. Results and Discussion'));
c.push(H2('3.1 Test set'));
c.push(P('To demonstrate the utility of this tool, we assembled a test set of 77 molecules with three to five heavy atoms, extracted from ATcT, for which accurate enthalpies of formation are known (Table S2; selection criteria in Section 3.8), spanning −138 to +156 kcal mol^{−1}. The set was subdivided into 55 unstrained (acyclic) and 22 strained species, each containing a three-, four- or five-membered ring (for example cyclopropane, cyclobutene, cyclopentadiene, aziridine and spiropentadiene, and the polycyclic tetrahedrane and [1.1.1]propellane). Each molecule was treated with the isodesmic (ID), hypohomodesmotic (HY) and homodesmotic (HO) schemes at three levels of theory (B3LYP/6-31+G(d), M06-2X/6-31+G(d) and ωB97X-D/def2-TZVP) and, after single-point refinement, on the three geometries, giving 231 reaction evaluations (165 unstrained, 66 strained). This hierarchical scheme allows direct assessment of (i) functional dependence, (ii) geometric effects and (iii) the intrinsic error cancellation achieved by each reaction type. Because the molecules are small, the three schemes often coincide (Section 3.5).'));
c.push(H2('3.2 Performance of the reaction schemes before refinement'));
c.push(P('Table 2 gives error statistics for the three levels of theory before single-point refinement. Overall, homodesmotic reactions have the lowest RMSE (2.027 kcal mol^{−1}) and MAE (0.978), followed by hypohomodesmotic (2.118; 1.098) and isodesmic (2.428; 1.853). Homodesmotic reactions have the lowest MAE at each level of theory and the lowest RMSE with M06-2X/6-31+G(d) and ωB97X-D/def2-TZVP geometries; the exception is the RMSE with B3LYP/6-31+G(d), where ID (2.69 kcal mol^{−1}) is below HO (2.91) and HY (3.04) (Table S1). The B3LYP exception is caused by three polycyclic species (Section 3.6). The ordering is driven by the acyclic species: for these, RMSE falls from 2.045 (ID) to 0.658 (HY) and 0.488 kcal mol^{−1} (HO). For the strained species all schemes perform poorly (RMSE 3.19, 3.82 and 3.71 kcal mol^{−1}), and before refinement the higher-order schemes do not improve on the isodesmic one. The corresponding error distributions are shown in Figure S1.'));
c.push(TCAP('**Table 2.** Error statistics (kcal mol^{−1}) relative to ATcT (0 K) for ID, HY and HO reactions, before single-point refinement, pooled over B3LYP/6-31+G(d), M06-2X/6-31+G(d) and ωB97X-D/def2-TZVP. *N* is the number of reactions (molecules × 3 levels).'));
c.push(table(T_COLS, T_HDR, statRows('opt')));
c.push(P('', { para: { spacing: { after: 60 } } }));
c.push(H2('3.3 Effect of single-point refinement'));
c.push(P('Replacing the lower-level electronic energy with an ωB97X-V/def2-TZVP single-point energy on the same geometry (Table 3, Figures 2 and 3) reduces the overall RMSE for every scheme, mainly by reducing the spread of the errors. The one exception is the hypohomodesmotic scheme for acyclic molecules, whose RMSE changes little (0.658 to 0.674 kcal mol^{−1}). Overall RMSE decreases from 2.428 to 1.827 (ID, 25%), from 2.118 to 1.211 (HY, 43%) and from 2.027 to 1.099 kcal mol^{−1} (HO, 46%); MAE decreases from 1.853 to 1.270, 1.098 to 0.768, and 0.978 to 0.641 kcal mol^{−1}, respectively. The hypohomodesmotic and homodesmotic schemes improve most, because refinement removes the large functional errors of their reactions for the polycyclic species (Section 3.6), and homodesmotic reactions remain the most accurate both before and after refinement. After refinement, 85.7% of homodesmotic errors lie within 1 kcal mol^{−1} of ATcT (HY: 78.4%; ID: 53.7%), and the largest homodesmotic error is 4.9 kcal mol^{−1} (HY: 4.9; ID: 8.1). The mean signed error of the HY scheme changes from −0.08 to +0.23 kcal mol^{−1}, indicating a small shift in bias.'));
c.push(TCAP('**Table 3.** As Table 2, after ωB97X-V/def2-TZVP single-point refinement (SPR) on the geometries from the three optimization levels.'));
c.push(table(T_COLS, T_HDR, statRows('sp')));
c.push(P('', { para: { spacing: { after: 60 } } }));
c.push(P('The violin and box plots show narrower error distributions and smaller ranges for HO reactions, consistent with the additional conservation constraints improving thermochemical consistency: the homodesmotic distribution is the narrowest and is centered closest to zero, especially for acyclic molecules (Figure 2a), where it lies almost entirely within ±1 kcal mol^{−1}. For strained rings the distributions are broad for all schemes (Figure 2b). Figure 3 compares RMSE before and after refinement with bootstrap intervals; the intervals for the strained set are wide because it contains only 22 molecules.'));
c.push(fig('fig2_distributions.png', 600));
c.push(CAP('**Figure 2.** Distribution of post-SPR errors relative to ATcT for (a) the 55 unstrained (acyclic) and (b) the 22 strained molecules, pooled over the three geometry levels (165 and 66 reactions). Violins are Gaussian kernel density estimates, which taper slightly beyond the extreme observations; boxes show the interquartile range and median; points are individual reactions; the shaded band marks ±1 kcal mol^{−1}.'));
c.push(fig('fig3_rmse.png', 600));
c.push(CAP('**Figure 3.** RMSE before (hatched) and after (solid) single-point refinement, pooled over the three geometry levels, for (a) unstrained, (b) strained and (c) all molecules. Whiskers are 95% bootstrap intervals obtained by resampling molecules.'));
c.push(H2('3.4 Dependence on the optimization level'));
c.push(P('Refinement also removes most of the dependence on the level used for the geometry. Across the three geometries, the spread (maximum minus minimum) of the homodesmotic error for a given molecule averages 1.47 kcal mol^{−1} before refinement (maximum 14.8, for [1.1.1]propellane) but only 0.15 kcal mol^{−1} after (maximum 1.4; Figure 4). Once a better electronic energy is used, inexpensive geometries are therefore sufficient. Per-level statistics are given in the Supporting Information (Table S1).'));
c.push(fig('fig4_geometry.png', 600));
c.push(CAP('**Figure 4.** (a, b) Homodesmotic errors of each molecule obtained with B3LYP/6-31+G(d) and ωB97X-D/def2-TZVP geometries, before and after single-point refinement (note the different axis ranges); max |Δ| is the largest difference between the two. (c) Cumulative distribution of the per-molecule spread of the error across the three geometry levels for each scheme.'));
c.push(H2('3.5 What does the homodesmotic scheme add?'));
c.push(P('Across all 231 reaction evaluations, homodesmotic reactions yield the lowest pooled RMSE and MAE after refinement. The advantage over isodesmic reactions is large and holds molecule by molecule; the advantage over hypohomodesmotic reactions is smaller and less certain. Relative to isodesmic reactions, HO lowers the overall RMSE by 17% before and 38% after refinement (MAE: 46% and 47%) and has the smaller mean absolute error for 47 of 55 acyclic molecules after refinement (*p* < 0.001, sign test over molecules), but for only 12 of 22 strained molecules.'));
c.push(P('The comparison with hypohomodesmotic reactions is more nuanced. HY and HO return the same reaction for 51 of the 77 molecules and differ for 26 (16 acyclic, 10 strained), so any advantage of HO rests on these. This is largely a consequence of molecular size: with three to five heavy atoms, the three-atom fragments that define the two schemes often coincide, and they separate only when conjugation or branching creates distinct fragments (as for C_{60}; Section 3.7). The RMSE advantage of HO is correspondingly modest, 4% before and 9% after refinement. Counted over molecules, HO is better for 22 of the 26 molecules before refinement (*p* = 0.001) but for only 16 after it (*p* = 0.33); for acyclic molecules the post-SPR RMSE falls from 0.67 to 0.35 kcal mol^{−1}, an interval for the difference that only just excludes zero, whereas the interval for the strained set includes it (Note S1). We therefore regard the advantage over hypohomodesmotic reactions as suggestive for acyclic molecules and unresolved otherwise, and a benchmark of larger molecules, where the schemes diverge more often, is needed.'));
c.push(P('These differences should be read against the uncertainty of the reference data. The ATcT uncertainties of the 77 targets have a median of 0.20 kcal mol^{−1} (maximum 0.55), and propagating those of the reference species adds a median of 0.14, together a root-mean-square uncertainty of about 0.31 kcal mol^{−1} (Table S2). The post-SPR RMSE of 0.35 kcal mol^{−1} for acyclic molecules is therefore close to the limit set by the reference values, whereas the errors of the strained rings, up to 4.9 kcal mol^{−1}, are far larger than any ATcT uncertainty.'));
c.push(H2('3.6 Strained molecules'));
c.push(P('All schemes degrade for strained rings (Figure 5; Tables 2 and 3). After refinement the pooled RMSE for the strained set is 2.84 (ID), 2.00 (HY) and 1.98 kcal mol^{−1} (HO), compared with 1.21, 0.67 and 0.35 kcal mol^{−1} for the acyclic set. For the strained species the HY and HO mean signed errors are negative (about −1.2 kcal mol^{−1} before refinement; HO −0.42 after), whereas the isodesmic errors are positive, so the schemes err in different directions. The largest post-SPR homodesmotic errors (2.7–4.9 kcal mol^{−1}) all belong to strained species, namely tetrahedrane, methylenecyclopropene, bicyclo[1.1.0]but-1(3)-ene, spiropentadiene and cyclopropenone. Since increasing molecular strain commonly increases the deviation from the typical behavior underlying the criteria of the homodesmotic family, the likelihood of accurate estimates should diminish; ring strain is not conserved by fragment-based schemes unless the reference species contain comparable strained rings. While the errors do increase, the higher-order schemes nevertheless remain centered near the reference value, with several of the rings (cyclobutane, cyclobutene, oxetene, cyclopentane, cyclopentene) within 1 kcal mol^{−1}.'));
c.push(fig('fig5_strained.png', 470));
c.push(CAP('**Figure 5.** Post-SPR errors (mean over the three geometries) for the 22 strained molecules, sorted by the homodesmotic error. The shaded band marks ±1 kcal mol^{−1}; grey lines connect the three schemes for each molecule.'));
c.push(P('Three of the strained species, tetrahedrane (C_{4}H_{4}), bicyclo[1.1.0]but-1-ene (C_{4}H_{4}) and [1.1.1]propellane (C_{5}H_{6}), are polycyclic and show this limitation most clearly. No reference species of their hypohomodesmotic or homodesmotic reaction contains a ring: the cage is opened completely into acyclic fragments (for tetrahedrane, C_{4}H_{4} + 6 ethane → 4 isobutane), so the whole strain energy appears in the reaction enthalpy and must be described by the electronic-structure method, and for tetrahedrane and [1.1.1]propellane the HY and HO reactions are identical. Before refinement, the B3LYP/6-31+G(d) HY and HO errors of these three species are large and negative, −7.7 to −13.7 kcal mol^{−1}, while their isodesmic errors lie between −1.4 and +4.7 kcal mol^{−1}; with M06-2X/6-31+G(d) and ωB97X-D/def2-TZVP the HY and HO errors are between −8.0 (tetrahedrane, M06-2X) and +2.4 kcal mol^{−1}. The higher-order reactions produce several branched species (isobutane, isobutene, neopentane), and the B3LYP errors are consistent with its known underestimation of the stability of branched alkanes,[#47] which accumulates with every branched product; the isodesmic products (ethane, ethylene) are unbranched. These three species are why HY and HO do not improve on ID for the strained set before refinement (Table 2). Single-point refinement removes most of this error: after refinement the HO errors (mean over the three geometries) are −4.6 (tetrahedrane), +2.5 (bicyclo[1.1.0]but-1-ene) and +1.5 kcal mol^{−1} ([1.1.1]propellane), whereas the isodesmic errors are +0.2, +5.1 and +7.7 kcal mol^{−1}; tetrahedrane retains the largest homodesmotic error of the benchmark. The hierarchy therefore conserves bond types and hydrogen environments but not ring strain, and for polycyclic targets reliable values would require reference species that contain comparable rings.'));
c.push(H2('3.7 Scalability of reaction generation: the C₆₀ example'));
c.push(P('An advantage of this methodology is its scalability. For a large organic molecule, the computation space expands only for new fragments, which allows fast and efficient energy calculation; for example, C_{60} draws on the same small reference species as much smaller molecules, scaled only by the coefficients, which is significantly cheaper than computing a large set of C_{60}-sized species. In principle this logic extends to all balanceable organic molecules, not just hydrocarbons. We tested it by generating homodesmotic reactions for linear polyenes (C_{6}–C_{38}), polycyclic aromatic hydrocarbons (benzene to coronene) and C_{60} (Figure 6; Table S3). For every polyene the reaction involves only ethylene and 1,3-butadiene, with coefficients that grow linearly with chain length. The aromatic hydrocarbons and C_{60} are different: because the generator kekulizes aromatic rings, their reactions compare an aromatic target with non-aromatic, localized references and are not error-balanced (for benzene, benzene + 3 ethylene → 3 butadiene, whose reaction energy contains the aromatic stabilization energy[#18]). For these systems at most four reference species (ethylene, propene, 1,3-butadiene and isoprene) are needed, and the C_{60} reaction is'));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 60, after: 120 }, children: runs('C_{60} + 60 propene + 30 ethylene → 60 isoprene', { italics: false }) }));
c.push(P('which conserves all bond hybridization types and hydrogen environments and uses three small reference species that all have ATcT enthalpies of formation. Each reaction was generated in under 0.1 s; this does not test solver scaling but shows that the list of reference species is short and does not grow with the target (for C_{60}, the isodesmic and hypohomodesmotic reactions use methane/ethane/ethylene and ethane/ethylene/isobutene). Scaling the target changes only the stoichiometric coefficients, not the set of species that must be computed. Three caveats apply: the target itself must still be computed, accuracy for large systems is untested because the benchmark has at most five heavy atoms, and aromatic stabilization, curvature and ring strain are not represented by the reference species. The benzene, polycyclic aromatic and C_{60} reactions are therefore demonstrations of reaction generation, not validated predictions (Note S2).'));
c.push(fig('fig6_scalability.png', 600));
c.push(CAP('**Figure 6.** Scalability of homodesmotic reaction generation. (a) Number of distinct reference species and (b) total number of reference molecules (sum of stoichiometric coefficients) versus the number of heavy atoms in the target for linear polyenes, polycyclic aromatic hydrocarbons and C_{60}; all 15 reactions were each generated in under 0.1 s. Only the polyene reactions are error-balanced; the aromatic and C_{60} reactions are not (see text). (c) The homodesmotic reaction for C_{60}.'));

// 4 Conclusions
c.push(H2('3.8 Scope of the benchmark'));
c.push(P('The 77 benchmark molecules are small and chemically narrow: 75 contain carbon, 32 oxygen and 17 fluorine or chlorine, but only four contain nitrogen (glycine, propylamine, *N*-methylethylamine and aziridine), and none contains a nitrile, amide or nitro group or is a benzenoid aromatic. They were selected from ATcT using two criteria: each can be decomposed by all three schemes (a non-trivial reaction exists at each level), and every reference species in those reactions has an ATcT enthalpy of formation, so that the computed value can be verified. Of the 229 neutral, closed-shell ATcT species with three to five heavy atoms that contain only C, H, N, O, F and Cl, 83 satisfied both criteria. Seventy-seven form the benchmark, and six were excluded: two enols whose *E*/*Z* isomer is not specified in ATcT, the stereo-unspecified entries of 1,2-dichloro- and 1,2-difluoroethylene, which duplicate the (*Z*) isomers in the benchmark, the methanol dimer, and pyramidane (C_{5}H_{4}), for which the optimizations of two of its three starting conformers did not complete at ωB97X-D/def2-TZVP (the status of every species is listed in data/selection_check.csv of the repository). The benchmark is limited by the availability of ATcT reference values and the cost of the calculations, not by the method: reactions are built from the fragments of the target, so in principle the approach applies to any covalent molecule whose reference species have reliable enthalpies of formation. Homodesmotic reactions cannot be written for the smallest molecules, which are their own smallest fragment, and among the 229 small species surveyed a homodesmotic reaction exists for 103, compared with 228 for the isodesmic scheme. Larger molecules contain more distinct fragments, so we expect the gap to narrow with size, but this was not tested here, and unusual bonding, stereochemical or SMILES cases may require the wildcard fragment library to be extended. The conclusions nevertheless apply to small compounds of C, H, O, F and Cl, and the 22 strained molecules are a small and varied set.'));
c.push(P('The comparison made here is deliberately between rungs of the hierarchy generated by one framework: isodesmic reactions, on which automated workflows such as the package of Wu *et al.*[#33] are built, and the homodesmotic reactions that the present framework makes attainable. It shows the difference between the schemes; it is not a comparison between software packages, CBH-*n*,[#28–30] or composite methods. Nor is a cost comparison possible at this size: across the 77 benchmark reactions the isodesmic scheme needs 13 distinct reference species of at most two heavy atoms, whereas the hypohomodesmotic and homodesmotic schemes need 44 and 52 species of up to five heavy atoms, as large as the targets; any saving from small references appears only for larger targets.'));
c.push(H1('4. Conclusion'));
c.push(P('This framework demonstrates that automatically generated homodesmotic reactions give the lowest errors on a 77-molecule ATcT benchmark: an overall RMSE of 1.10 kcal mol^{−1} after single-point refinement (0.35 for acyclic molecules), an improvement over isodesmic reactions that holds for individual molecules (47 of 55 acyclic molecules) and over hypohomodesmotic reactions that is supported before refinement but comparable to the uncertainty of the reference data afterwards. For small strained rings all schemes err by about 2–3 kcal mol^{−1} and cannot be statistically separated; for polycyclic cages, whose ring strain no scheme conserves, the higher-order schemes depend strongly on the functional before refinement and still err by up to 4.9 kcal mol^{−1} after it. Single-point refinement reduces errors for all schemes and makes the results nearly independent of the optimization level. The reference species depend on local bonding motifs, not molecular size, as shown for polyenes; for C_{60} this demonstrates reaction generation, not prediction, because accuracy for large systems is untested. These data show that homodesmotic reactions are more accurate than isodesmic reactions at modest levels of theory, especially for systems that are not severely strained. They suggest a smaller advantage over hypohomodesmotic reactions, which a benchmark of larger molecules is needed to confirm. By redefining error-balanced reactions as a structured set of linear constraints, the method enables the automated generation of rigorously balanced reactions that capture detailed bonding and hybridization environments. This formulation simplifies a historically tedious and error-prone process and provides a foundation for systematically extending error-balanced methodologies to higher-order schemes.'));

c.push(H1('Supporting Information'));
c.push(P('The Supporting Information contains per-level error statistics (Table S1), the 77 test molecules with their ATcT values and uncertainties (Table S2), the scalability reactions (Table S3), the hypohomodesmotic reactions that depend on the ILP objective (Table S4), the wildcard fragment libraries (Figures S2 and S3, Table S5), every benchmark reaction (Table S6), the pre-refinement error distributions (Figure S1), two notes (S1, statistical comparison of the schemes; S2, limitations of the C₆₀ example), and the absolute energies (Table S7) and Cartesian coordinates of all optimized structures.'));
c.push(H1('Acknowledgements'));
c.push(P('This work was conducted using high-performance computing resources provided by the Georgia Advanced Computing Resource Center (http://gacrc.uga.edu). [[AUTHOR: add funding sources and grant numbers.]]'));
c.push(H1('Conflict of Interest'));
c.push(P('The authors declare no conflict of interest.'));
c.push(H1('Data Availability Statement'));
c.push(P('The software is openly available at https://github.com/aaronseymour7/HFR_test [[AUTHOR: license, release tag and archive DOI (e.g., Zenodo); cite the archive in the reference list]]. The paper/ directory of the repository contains the errors for all 462 benchmark reactions (data/benchmark_errors.csv), the scalability data (data/scaling.csv), the wildcard libraries (data/wildcard_libraries.csv) the selection check (data/selection_check.csv) and the coordinates and absolute energies of the optimized structures (data/si/), from which recompute_from_structures.py reproduces every benchmark error, together with the scripts that generate every figure (make_figures.py), the error statistics quoted in the text (paper_statistics.py), the scalability test (scaling_demo.py), the selection check (selection_check.py) and the objective check (objective_check.py). The conformer generation and selection code reproduces the procedure used for the benchmark and is tested against the original scripts. The benchmark reactions were generated with the original HFR_PACKAGE (https://github.com/aaronseymour7/HFR_PACKAGE). [[AUTHOR: deposit the output files (energies, frequencies and coordinates) in ioChem-BD or Zenodo and give the link and citation here.]]'));
c.push(H1('References'));
const REF_AT = c.length; // reference list is inserted here once every citation has been numbered

c.push(new Paragraph({ children: [new PageBreak()] }));
c.push(H1('Graphical Abstract'));
c.push(P('Homodesmotic reactions, written automatically as integer linear programs, turn modest DFT energies into enthalpies of formation within 0.4 kcal mol^{−1} of Active Thermochemical Tables values for acyclic molecules. The reference species depend only on local bonding motifs, so the same small molecules serve targets from butanone to C_{60}.'));
c.push(P('[[AUTHOR: graphical abstract image, 5.5 × 5.0 cm or 11.5 × 2.5 cm; for example the homodesmotic reaction of butan-2-one from Figure 1 next to the C₆₀ reaction of Figure 6c.]]'));

// ---- SI
c.push(new Paragraph({ children: [new PageBreak()] }));
c.push(H1('Supporting Information'));
c.push(P('**Contents.** Note S1. Statistical comparison of the reaction schemes and reference uncertainty · Note S2. Limitations of the C₆₀ and aromatic examples · Table S1. Error statistics by level of theory · Table S2. The 77 test molecules · Table S3. Scalability reactions · Table S4. Objective-dependent hypohomodesmotic reactions · Figure S1. Pre-refinement error distributions · Figures S2 and S3, Table S5. Wildcard fragment libraries · Table S6. All benchmark reactions · Table S7. Absolute energies · Cartesian coordinates', { align: AlignmentType.LEFT }));
c.push(H2('Note S1. Statistical comparison of the reaction schemes and reference uncertainty'));
c.push(P('The pooled RMSE difference (ID − HO) after refinement is 0.73 kcal mol^{−1} (95% bootstrap interval 0.18–1.28); for acyclic molecules it is 0.86 (0.62–1.10). The HY and HO schemes return the same reaction (identical errors) for 51 of the 77 molecules; they differ for 26 (16 acyclic and 10 strained), so any advantage of HO rests on these. This is largely a consequence of molecular size: with only three to five heavy atoms, the three-atom fragments that define hypohomodesmotic and homodesmotic reactions often coincide, and the schemes separate only when conjugation or branching creates distinct fragments (for C_{60} they give different reference species; Section 3.7). The overall RMSE advantage of HO is therefore modest, 4% before and 9% after refinement (MAE: 11% and 17%). Counting molecules rather than reactions, HO has the smaller mean absolute error for 14 of the 16 acyclic molecules on which the schemes differ before refinement (*p* = 0.004, two-sided sign test) but for only 12 of 16 after refinement (*p* = 0.08). For the strained molecules the corresponding counts are 8 of 10 and 4 of 10 (*p* = 0.11 and 0.75), and over all 26 molecules they are 22 and 16 (*p* = 0.001 and 0.33). After refinement the RMSE for acyclic molecules falls from 0.67 to 0.35 kcal mol^{−1} (HY − HO = 0.33 kcal mol^{−1}); the bootstrap interval, which resamples molecules, is 0.07–0.60 kcal mol^{−1} and only just excludes zero. For the strained species the interval is −0.35 to 0.44, and over the full set after refinement it is −0.10 to 0.36 kcal mol^{−1} (0.11 kcal mol^{−1}); both include zero. We therefore regard the advantage of homodesmotic over hypohomodesmotic reactions as suggestive for acyclic molecules and unresolved otherwise. Against isodesmic reactions the evidence is stronger: after refinement HO has the smaller mean absolute error for 47 of 55 acyclic molecules (*p* < 0.001) and for 59 of 77 overall, but for only 12 of 22 strained molecules. A benchmark of larger molecules, where the schemes diverge more often, is needed to establish the homodesmotic advantage more broadly. These comparisons should be read against the uncertainty of the reference data. The ATcT uncertainties of the 77 target molecules have a median of 0.20 kcal mol^{−1} (maximum 0.55, for spiropentadiene; cyclobuta-1,2-diene (0.53) and 1,2,4-trioxolane (0.50) are the only others at 0.5 or above), and propagating the uncertainties of the reference species of the homodesmotic reactions adds a median of 0.14 kcal mol^{−1} (maximum 0.39); together they correspond to a root-mean-square uncertainty of about 0.31 kcal mol^{−1} (Table S2). The RMSE of 0.35 kcal mol^{−1} for acyclic molecules after refinement is therefore close to the limit set by the reference values, and differences of a few tenths of a kcal mol^{−1} between schemes (such as the 0.33 kcal mol^{−1} difference between HY and HO for acyclic molecules) are of the same size as that uncertainty. The errors of the strained rings, up to 4.9 kcal mol^{−1}, are far larger than any ATcT uncertainty, so the conclusions for those species are not affected.'));
c.push(H2('Note S2. Limitations of the C₆₀ and aromatic examples'));
c.push(P('Several limitations should be stated clearly. First, the target itself must still be computed: what is saved is the cost of the *reference* species, which are small and can be shared across a whole dataset (the unique-molecule workflow of the software does this), not the cost of the target\'s own optimization and frequency calculation. Second, accuracy is not established for large systems. The benchmark contains molecules with at most five heavy atoms, and for a target like C_{60} the coefficients of 30 to 60 multiply the residual error of each reference species, so absolute errors in the enthalpy of formation may grow with size even when the error per fragment is small. Third, aromatic stabilization, curvature and ring strain (C_{60} has twelve pentagons) are not represented by the non-aromatic, unstrained reference species, and our strained-ring results (errors of 2–3 kcal mol^{−1}, and up to 4.9 for polycyclic cages after refinement) suggest that this cannot be neglected. We therefore present the benzene, polycyclic aromatic and C_{60} reactions as demonstrations that the reference set stays short, not as validated predictions, and C_{60} as an extreme case of reaction generation, not of accuracy.'));
c.push(TCAP('**Table S1.** Error statistics (kcal mol^{−1}) by level of theory. "Opt" = lower-level energies; "SPR" = ωB97X-V/def2-TZVP single points on the geometry of that level. *N* = 55 (unstrained), 22 (strained), 77 (combined) per level.'));
const lv = { B3LYP631GD: 'B3LYP/6-31+G(d)', M062X631GD: 'M06-2X/6-31+G(d)', wB97DEF2tzvp: 'ωB97X-D/def2-TZVP' };
const lvSP = { SPB3L: 'B3LYP/6-31+G(d)', SPM06: 'M06-2X/6-31+G(d)', SPwB97: 'ωB97X-D/def2-TZVP' };
const s1rows = [];
for (const [lab, map] of [['opt', lv], ['sp', lvSP]]) {
  for (const [key, name] of Object.entries(map)) {
    for (const [ds, dn] of [['unstrained', 'Unstr.'], ['strained', 'Str.'], ['combined', 'All']]) {
      for (const [k, i] of [['ID', 0], ['HY', 1], ['HO', 2]]) {
        const x = stats[`${lab}|${ds}|${key}|${k}`];
        s1rows.push([i === 0 && ds === 'unstrained' ? `${lab === 'opt' ? 'Opt' : 'SPR'}: ${name}` : '', i === 0 ? dn : '', k, f3(x.rmse), f3(x.mae), f3(x.me), f3(x.std)]);
      }
    }
  }
}
// Custom 7-col table for S1
{
  const cols = [2600, 900, 800, 1181, 1181, 1182, 1182];
  const hdr = new TableRow({ tableHeader: true, children: ['Level', 'Set', 'Scheme', 'RMSE', 'MAE', 'ME', 'STD'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: i < 2 ? AlignmentType.LEFT : AlignmentType.CENTER })) });
  const body = s1rows.map((r, ri) => new TableRow({ cantSplit: true, children: r.map((t, i) => cell(t, cols[i], { align: i < 2 ? AlignmentType.LEFT : AlignmentType.CENTER, bottom: ri === s1rows.length - 1, keep: false })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
}
c.push(new Paragraph({ children: [new PageBreak()] }));
c.push(TCAP('**Table S2.** The 77 test molecules: classification (U = unstrained/acyclic, S = strained, ring sizes given), SMILES, formula, ATcT 0 K enthalpy of formation and its listed uncertainty (±), and the uncertainty propagated from the ATcT values of the reference species of the homodesmotic reaction (ref. ±; quadrature, covariance ignored). All energies in kcal mol^{−1}.'));
{
  const cols = [500, 650, 2576, 1400, 800, 1100, 1000, 1000];
  const hdr = new TableRow({ tableHeader: true, children: ['#', 'Class', 'SMILES', 'Formula', 'Rings', 'ATcT', '±', 'Ref. ±'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: i === 2 ? AlignmentType.LEFT : AlignmentType.CENTER })) });
  const body = mols.map((m, ri) => new TableRow({ cantSplit: true, children: [String(m.i), m.cls === 'strained' ? 'S' : 'U', m.smiles, m.formula, m.rings.length ? m.rings.join(',') : '–', m.atct.toFixed(2), m.sigma.toFixed(2), m.sigma_ref.toFixed(2)].map((t, i) => cell(t, cols[i], { align: i === 2 ? AlignmentType.LEFT : AlignmentType.CENTER, bottom: ri === mols.length - 1, keep: false })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
}

c.push(new Paragraph({ children: [new PageBreak()] }));
c.push(fig('figS1_distributions_preSPR.png', 600));
c.push(CAP('**Figure S1.** As Figure 2 (violins are kernel density estimates, so their tails extend slightly past the data), before single-point refinement (lower-level energies).'));
c.push(TCAP('**Table S3.** Homodesmotic reactions generated for the scalability test (Figure 6). Species are given as SMILES. For every reaction, conservation of bond hybridization types and hydrogen environments was verified independently of the ILP.'));
{
  const rowsS3 = fs.readFileSync('data/scaling.csv', 'utf8').trim().split(/\r?\n/).slice(1).map(l => {
    const f = l.split(',');
    return [f[1], f[3], f[4], f.slice(9).join(',')];
  });
  const cols = [1350, 800, 800, 6076];
  const hdr = new TableRow({ tableHeader: true, children: ['Target', 'Heavy atoms', 'Solve (s)', 'Reaction'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: i === 3 ? AlignmentType.LEFT : AlignmentType.CENTER })) });
  const body = rowsS3.map((r, ri) => new TableRow({ cantSplit: true, children: r.map((t, i) => cell(i === 3 ? t.replace(/->/g, '→') : t, cols[i], { align: i === 3 ? AlignmentType.LEFT : AlignmentType.CENTER, bottom: ri === rowsS3.length - 1, keep: false })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
}

c.push(TCAP('**Table S4.** Hypohomodesmotic reactions for the four molecules whose reaction depends on the ILP objective. "No objective" is the reaction used in the benchmark; "minimal" is the reaction with the fewest species. All other 73 hypohomodesmotic reactions, and all isodesmic and homodesmotic reactions, are identical in both formulations.'));
{
  const rec = JSON.parse(fs.readFileSync('data/hy_ambiguous.json', 'utf8'));
  const cols = [1400, 3813, 3813];
  const hdr = new TableRow({ tableHeader: true, children: ['Molecule (SMILES)', 'No objective (benchmark)', 'Minimal number of species'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER })) });
  const body = rec.map((r, ri) => new TableRow({ cantSplit: true, children: [r[0], r[1].replace(/->/g, '→'), r[2].replace(/->/g, '→')].map((t, i) => cell(t, cols[i], { align: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER, bottom: ri === rec.length - 1, keep: false })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
}

c.push(new Paragraph({ children: [new PageBreak()] }));
c.push(fig('figS2_wildcard_libraries.png', 600));
c.push(CAP('**Figure S2.** Wildcard fragment libraries of the isogyric (level 1), isodesmic (level 2) and hypohomodesmotic (level 3) schemes drawn as Kekulé structures, with the SMILES/SMARTS of each fragment beneath it. X stands for any atom (a [\\*] in the pattern); hydrogens are not drawn and are added to complete the valence. Co-reactants are the candidate left-hand-side species and products the candidate right-hand-side species of the ILP. The isodesmic code builds its candidates directly from the bond types of the target instead of mapping this library, and adds H_{2} explicitly (see Table S5).'));
c.push(fig('figS3_wildcard_library_homodesmotic.png', 600));
c.push(CAP('**Figure S3.** Wildcard fragment library of the homodesmotic scheme (level 4), drawn as in Figure S2. Y stands for any non-carbon atom ([!#6] in the SMARTS patterns).'));
c.push(TCAP('**Table S5.** Wildcard fragment libraries of each level as SMILES/SMARTS (drawn in Figures S2 and S3; hfrpkg/core.py, *get_fragments*). Each [\\*] stands for any atom: a pattern is mapped onto the kekulized target, each [\\*] is replaced by the element of the matched atom, hydrogens are added, and the distinct species obtained are the candidate right-hand-side ("products") or left-hand-side ("reactants") species of the ILP. Patterns are SMILES except where marked SMARTS (which allow, for example, a non-carbon atom [!#6]). Level 2 lists the isodesmic candidates in generic form: the isodesmic code builds them directly from the bond types of the target (the same patterns are defined, but not mapped, in *get_fragments*) and adds the H_{2} co-reactant explicitly.'));
{
  const rows = JSON.parse(fs.readFileSync('data/wildcards.json', 'utf8'));
  const cols = [800, 2200, 1500, 3526, 1000];
  const hdr = new TableRow({ tableHeader: true, children: ['Level', 'Scheme', 'Role', 'Pattern', 'Syntax'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: i === 3 ? AlignmentType.LEFT : AlignmentType.CENTER })) });
  const body = rows.map((r, ri) => new TableRow({ cantSplit: true, children: [String(r[0]), r[1], r[2], r[3].replace(/\*/g, '\\*'), r[4]].map((t, i) => cell(t, cols[i], { align: i === 3 ? AlignmentType.LEFT : AlignmentType.CENTER, bottom: ri === rows.length - 1, keep: false })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
}

c.push(TCAP('**Table S6.** The isodesmic (ID), hypohomodesmotic (HY) and homodesmotic (HO) reaction of every benchmark molecule, as SMILES (*n* X = *n* molecules of X), as run in the benchmark. data/objective_check.csv lists the reactions the current code returns, which are the same except for the homodesmotic reactions of glycine and cyclopropenone and the hypohomodesmotic and homodesmotic reactions of HClO_{3}; for HClO_{3}, the runs at different levels returned different (equally balanced) reactions, as indicated. Recomputing all 1386 errors from these reactions, the energies of Table S7 and the ATcT values reproduces data/benchmark_errors.csv exactly (paper/recompute_from_structures.py).'));
{
  const lines = fs.readFileSync('data/objective_check.csv', 'utf8').trim().split(/\r?\n/).slice(1).map(l => l.split(','));
  const byMol = new Map();
  for (const [smi, scheme, , noImine, current] of lines) {
    if (!byMol.has(smi)) byMol.set(smi, {});
    byMol.get(smi)[scheme] = (scheme === 'HY' ? noImine : current).replace(/ -> /g, ' → ');
  }
  // reactions actually run in the benchmark where they differ from the regenerated ones (index.txt of the run folders)
  const A = 'OCl + OOOCl → 2 OOCl', B = 'OO + OOOCl → OOCl + OOO';
  byMol.get('NCC(=O)O').HO = 'CC + CO + NCC(=O)O → CN + CC(=O)O + CCO';
  byMol.get('O=c1cc1').HO = 'CC=O + C=O + C=CC + O=c1cc1 → 2 C=CC=O + CC(C)=O';
  byMol.get('OOOCl').HY = `${A} (B3LYP, M06-2X); ${B} (ωB97X-D)`;
  byMol.get('OOOCl').HO = `${A} (B3LYP, ωB97X-D); ${B} (M06-2X)`;
  const order = mols.map(m => m.smiles);
  const rows = order.map(smi => [smi, byMol.get(smi).ID, byMol.get(smi).HY, byMol.get(smi).HO]);
  const cols = [1300, 2242, 2742, 2742];
  const hdr = new TableRow({ tableHeader: true, children: ['Molecule', 'ID', 'HY', 'HO'].map((h, i) => cell(h, cols[i], { bold: true, top: true, bottom: true, align: AlignmentType.LEFT })) });
  const body = rows.map((r, ri) => new TableRow({ cantSplit: true, children: r.map((t, i) => cell(t, cols[i], { align: AlignmentType.LEFT, bottom: ri === rows.length - 1, keep: false })) }));
  c.push(new Table({ width: { size: W, type: WidthType.DXA }, columnWidths: cols, rows: [hdr, ...body] }));
}

// ---- coordinates and absolute energies, from the export_geometries.py output in data/si/ (species of the benchmark reactions only)
const PYRAMIDANE = 'InChI=1S/C5H4/c1-2-4-3(1)5(1,2)4/h1-4H'; // excluded from the benchmark (Section 3.8)
const LEVELS = [[/b3lyp/i, 'B3LYP/6-31+G(d)'], [/m06/i, 'M06-2X/6-31+G(d)'], [/b97/i, 'ωB97X-D/def2-TZVP']];
const structures = [];
{
  // one structure per species and level: the lowest E + ZPVE over all runs (on a tie the later run), as in
  // recompute_from_structures.py. Reruns (reoptimized saddle points, and species whose copy from the later
  // new-species runs is used in place of the original one) are marked with an asterisk.
  const REOPT = /^(b3lyp|m062x|wb97xd)$/;
  const best = new Map(), inFirst = new Set();
  for (const f of fs.readdirSync('data/si').filter(f => f.endsWith('.xyz')).sort()) {
    const L = fs.readFileSync('data/si/' + f, 'utf8').split(/\r?\n/);
    for (let i = 0; i < L.length && L[i].trim(); ) {
      const n = parseInt(L[i], 10), comment = L[i + 1], atoms = L.slice(i + 2, i + 2 + n);
      i += 2 + n;
      const [runDir, conformer, smiles, inchi] = comment.split(' ');
      const kv = Object.fromEntries([...comment.matchAll(/(\S+)=(\S+)/g)].filter(m => !m[1].startsWith('InChI')).map(m => [m[1], m[2]]));
      const level = LEVELS.find(([re]) => re.test(runDir))[1];
      if (inchi === PYRAMIDANE) continue;
      const key = level + ' ' + inchi, later = f.includes('new_species');
      if (!later) inFirst.add(key);
      const prev = best.get(key), e = Number(kv['E+ZPVE']);
      if (!prev || e <= Number(prev['E+ZPVE'])) best.set(key, { smiles, inchi, level, atoms, conformer, runDir, later, ...kv });
    }
  }
  for (const [key, s] of best) {
    s.rerun = REOPT.test(s.runDir) || (s.later && inFirst.has(key));
    s.label = s.smiles + (s.rerun ? '*' : '');
    structures.push(s);
  }
  const lv = (s) => LEVELS.findIndex(([, name]) => name === s.level);
  structures.sort((a, b) => a.smiles.localeCompare(b.smiles) || lv(a) - lv(b));
}
const nImag = structures.filter(s => s.n_imag !== '0').length;

c.push(TCAP(`**Table S7.** Absolute energies (hartree) of the lowest-energy conformer of every species: electronic energy (*E*) and unscaled zero-point vibrational energy (ZPVE) at the optimization level, and the ωB97X-V/def2-TZVP single-point energy (*E*_{SP}) on that geometry; *n*_{imag} is the number of imaginary frequencies. Where a species was computed in more than one run, the lowest-energy copy is listed and used in every reaction; an asterisk marks a rerun, either a structure reoptimized after first converging to a saddle point or a reference species whose copy from the later runs for the new species is at least as low in energy as the original. ${structures.length} structures (${structures.filter(s => s.rerun).length} reruns); Cartesian coordinates follow.`));
{
  const cols = [1900, 2026, 1500, 1100, 1500, 1000];
  c.push(table(cols, ['Species', 'Level', '*E*', 'ZPVE', '*E*_{SP}', '*n*_{imag}'],
    structures.map(s => [s.label, s.level, s.E, s.ZPVE, s.E_sp || '–', s.n_imag])));
}

c.push(H2('Cartesian coordinates (Å)'));
{
  const mono = (t, o = {}) => new Paragraph({ keepNext: o.keepNext, keepLines: true, spacing: { after: o.after || 0, line: 240 }, children: [new TextRun({ text: t, font: 'Courier New', size: 16, bold: !!o.bold })] });
  for (const s of structures) {
    c.push(mono(`${s.label}  ${s.level}  E = ${s.E}  E_SP = ${s.E_sp || '–'}  n_imag = ${s.n_imag}`, { bold: true, keepNext: true }));
    s.atoms.forEach((a, k) => {
      const [el, x, y, z] = a.trim().split(/\s+/);
      const fmt = (v) => Number(v).toFixed(6).padStart(12);
      c.push(mono(`${el.padEnd(2)} ${fmt(x)} ${fmt(y)} ${fmt(z)}`, { keepNext: k < s.atoms.length - 1, after: k === s.atoms.length - 1 ? 160 : 0 }));
    });
  }
}

{
  const uncited = REFS.map((_, i) => i + 1).filter(k => !citeOrder.includes(k));
  if (uncited.length) throw new Error('uncited references: ' + uncited.join(', '));
  c.splice(REF_AT, 0, ...citeOrder.map(k => refPara(REFS[k - 1])));
}

const doc = new Document({
  creator: 'Aaron D. Seymour', title: 'Automated Computational Thermochemistry Using Homodesmotic Reactions',
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 26, bold: true, font: FONT }, paragraph: { spacing: { before: 280, after: 120 }, outlineLevel: 0 } },
      { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 22, bold: true, italics: true, font: FONT }, paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 1 } },
    ],
  },
  numbering: { config: [{ reference: 'refs', levels: [{ level: 0, format: LevelFormat.DECIMAL, text: '[%1]', alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 540 } } } }] }] },
  sections: [{
    properties: { page: { margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] }) },
    children: c,
  }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync('Seymour_Wheeler_Homodesmotic_manuscript.docx', b); console.log('written', b.length); });
