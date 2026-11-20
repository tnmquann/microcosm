/* microcosm command builder: a decision tree over nextflow_schema.json that writes a command, params.yaml or a
   Nextflow config. Plain JavaScript, no build step. Parameter types and defaults come from schema.json, written
   by docs/src/conf.py at every docs build. */
(function () {
  'use strict';

  var SCRIPT = document.currentScript;
  var mount = document.getElementById('mcb');
  if (!mount || !SCRIPT) return;

  var SCHEMA_URL = new URL('schema.json', SCRIPT.src).href;
  var PIPELINE = 'tnmquann/microcosm';
  var STORE_KEY = 'microcosm-builder-v1';
  var DOCS = SCRIPT.src.replace(/_static\/builder\/builder\.js.*$/, '');
  var REDUCED = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  /* ------------------------------------------------------------------ schema */

  var P = {};          // name -> {type, default, enum, description, format, group}
  var GROUPS = {};     // group id -> [names]
  var VERSION = '2.0.0';
  var HIDDEN = ['pipeline_root', 'help', 'help_full', 'show_hidden', 'config_profile_name',
    'config_profile_description', 'cleanup_patterns'];

  function def(name) { return P[name] ? P[name].default : undefined; }
  function same(a, b) {
    if ((a === undefined || a === null || a === '') && (b === undefined || b === null || b === '')) return true;
    return String(a) === String(b);
  }

  /* ------------------------------------------------------------------ state */

  function freshState() {
    return {
      route: null, step: 'route', seen: {},
      p: {},
      ui: { db: {}, lineage: 'filled', hostRef: 'auto', asmStart: 'reads', memUnit: 'GB', timeUnit: 'h' },
      nf: { profile: 'singularity', revision: VERSION, workDir: '', resume: true, bg: false, executor: 'local',
        queue: '', account: '', queueSize: '', cacheDir: '' },
      fmt: 'cmd'
    };
  }
  var S = freshState();

  function get(name) { return Object.prototype.hasOwnProperty.call(S.p, name) ? S.p[name] : def(name); }
  function on(name) { var v = get(name); return v === true || v === 'true'; }
  function setP(name, value) {
    if (value === '' || value === undefined || value === null || same(value, def(name))) delete S.p[name];
    else S.p[name] = value;
  }
  function db(key) { return S.ui.db[key] || 'auto'; }

  /* route helpers */
  function reads() { return S.route === 'reads'; }
  function asm() { return S.route === 'assembly'; }
  function genomes() { return S.route === 'genomes'; }
  function fromAssemblies() { return asm() && S.ui.asmStart === 'assemblies'; }
  function taxotools() {
    return String(get('readbased_taxotools') || 'microcosm').split(',').map(function (t) { return t.trim(); })
      .filter(Boolean);
  }
  function useMicrocosm() { return reads() && taxotools().indexOf('microcosm') >= 0; }
  function useSylph() { return reads() && taxotools().indexOf('sylph') >= 0; }
  function yachtOn() { return useMicrocosm() && !on('skip_yacht'); }
  // The pipeline also skips YACHT when several sourmash databases are given.
  function yachtRuns() { return yachtOn() && count(get('sourmash_database')) <= 1; }
  function binningOn() { return asm() && !on('skip_binning'); }
  function downstreamOn() { return genomes() || (asm() && on('enable_bacterial_downstream')); }
  function classifyOn() {
    if (genomes()) return S.ui.lineage === 'blank';
    return binningOn() && !on('skip_bin_classification');
  }
  function binqcTools() {
    return String(get('binqc_tool') || '').split(',').map(function (t) { return t.trim(); }).filter(Boolean);
  }
  function needCheckm2() {
    return binningOn() && ((!on('skip_binqc') && binqcTools().indexOf('checkm2') >= 0) ||
      (!on('skip_binning_refinement') && get('refine_tool') === 'binette'));
  }
  function tool(name) { return downstreamOn() && !on('skip_' + name); }
  function count(value) {
    return String(value || '').split(',').map(function (t) { return t.trim(); }).filter(Boolean).length;
  }

  /* ------------------------------------------------------------------ DOM helpers */

  function h(tag, attrs) {
    var el = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (k) {
        var v = attrs[k];
        if (v === null || v === undefined || v === false) return;
        if (k === 'class') el.className = v;
        else if (k === 'text') el.textContent = v;
        else if (k === 'html') el.innerHTML = v;
        else if (k.slice(0, 2) === 'on') el.addEventListener(k.slice(2), v);
        else el.setAttribute(k, v === true ? '' : v);
      });
    }
    for (var i = 2; i < arguments.length; i++) append(el, arguments[i]);
    return el;
  }
  function append(el, child) {
    if (child === null || child === undefined || child === false) return;
    if (Array.isArray(child)) { child.forEach(function (c) { append(el, c); }); return; }
    el.appendChild(typeof child === 'string' ? document.createTextNode(child) : child);
  }
  function code(text) { return h('code', null, text); }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }
  var uid = 0;
  function nextId(base) { uid += 1; return 'mcb-' + base + '-' + uid; }

  /* ------------------------------------------------------------------ controls */

  // Structural changes re-render the step; typing only refreshes the output.
  function commit(structural) { if (structural) renderStep(); refresh(); save(); }

  function field(opts, control) {
    var id = control.id || (control.id = nextId('f'));
    return h('div', { class: 'mcb-field' + (opts.wide ? ' is-wide' : '') },
      h('label', { class: 'mcb-label', for: id },
        opts.label, opts.required ? h('span', { class: 'mcb-req', 'aria-hidden': 'true' }, ' required') : null,
        opts.param ? h('span', { class: 'mcb-flag' }, '--' + opts.param) : null),
      control,
      opts.help ? h('p', { class: 'mcb-help' }, opts.help) : null);
  }

  // Free text: paths, names and pasted tool arguments (passed through as typed).
  function textInput(name, opts) {
    opts = opts || {};
    var input = h('input', {
      type: 'text', id: 'in-' + name, class: 'mcb-input' + (opts.mono !== false ? ' is-mono' : ''),
      value: Object.prototype.hasOwnProperty.call(S.p, name) ? S.p[name] : '',
      placeholder: opts.placeholder || (def(name) !== undefined && def(name) !== null ? String(def(name)) : ''),
      spellcheck: 'false', autocomplete: 'off', 'data-req': opts.required ? name : null,
      oninput: function (e) { setP(name, e.target.value.trim()); commit(false); markRequired(); }
    });
    return field({ label: opts.label, help: opts.help, param: name, required: opts.required, wide: opts.wide }, input);
  }

  function numberInput(name, opts) {
    opts = opts || {};
    var cur = Object.prototype.hasOwnProperty.call(S.p, name) ? S.p[name] : '';
    var input = h('input', {
      type: 'number', id: 'in-' + name, class: 'mcb-input is-num', value: cur, step: opts.step || 'any',
      min: opts.min, placeholder: def(name) !== undefined && def(name) !== null ? String(def(name)) : '',
      oninput: function (e) {
        var v = e.target.value === '' ? '' : Number(e.target.value);
        if (P[name] && P[name].type === 'integer' && v !== '') v = Math.round(v);
        setP(name, v); commit(false);
      }
    });
    return field({ label: opts.label, help: opts.help, param: name }, input);
  }

  function toggle(label, desc, checked, onChange, opts) {
    opts = opts || {};
    var id = nextId('t');
    return h('div', { class: 'mcb-toggle' + (opts.disabled ? ' is-disabled' : '') },
      h('button', {
        type: 'button', role: 'switch', id: id, class: 'mcb-switch', 'aria-checked': checked ? 'true' : 'false',
        disabled: opts.disabled || null, onclick: function () { onChange(!checked); commit(true); }
      }, h('span', { class: 'mcb-knob' })),
      h('label', { for: id, class: 'mcb-toggle-text' },
        h('span', { class: 'mcb-toggle-label' }, label),
        opts.param ? h('span', { class: 'mcb-flag' }, opts.param) : null,
        desc ? h('span', { class: 'mcb-toggle-desc' }, desc) : null));
  }
  // A boolean parameter, optionally shown inverted (skip_x as "Run x").
  function paramToggle(name, label, desc, invert) {
    var value = on(name);
    return toggle(label, desc, invert ? !value : value, function (v) { setP(name, invert ? !v : v); },
      { param: '--' + name });
  }

  function segmented(label, options, current, onPick, opts) {
    opts = opts || {};
    var group = h('div', { class: 'mcb-seg', role: 'radiogroup', 'aria-label': label });
    options.forEach(function (o) {
      var active = o.value === current;
      group.appendChild(h('button', {
        type: 'button', role: 'radio', 'aria-checked': active ? 'true' : 'false',
        class: 'mcb-seg-btn' + (active ? ' is-on' : ''), disabled: o.disabled || null, title: o.why || null,
        onclick: function () { if (!active) { onPick(o.value); commit(true); } }
      }, o.label));
    });
    var blocked = options.filter(function (o) { return o.disabled && o.why; });
    return h('div', { class: 'mcb-field' + (opts.wide ? ' is-wide' : '') },
      h('span', { class: 'mcb-label' }, label, opts.param ? h('span', { class: 'mcb-flag' }, '--' + opts.param) : null),
      group,
      blocked.map(function (o) { return h('p', { class: 'mcb-help mcb-why' }, o.label + ' is not available: ' + o.why + '.'); }),
      opts.help ? h('p', { class: 'mcb-help' }, opts.help) : null);
  }
  function paramSeg(name, label, options, help) {
    return segmented(label, options, String(get(name)), function (v) { setP(name, v); }, { param: name, help: help });
  }

  // Multi-select tiles: each item has on(), flip(v), label, note.
  function chips(label, items, help) {
    var wrap = h('div', { class: 'mcb-chips', role: 'group', 'aria-label': label });
    items.forEach(function (it) {
      var active = it.on();
      wrap.appendChild(h('button', {
        type: 'button', class: 'mcb-chip' + (active ? ' is-on' : ''), 'aria-pressed': active ? 'true' : 'false',
        onclick: function () { it.flip(!active); commit(true); }
      }, h('span', { class: 'mcb-chip-mark', 'aria-hidden': 'true' }),
      h('span', { class: 'mcb-chip-text' }, h('span', { class: 'mcb-chip-label' }, it.label),
        it.note ? h('span', { class: 'mcb-chip-note' }, it.note) : null)));
    });
    return h('div', { class: 'mcb-field is-wide' }, h('span', { class: 'mcb-label' }, label), wrap,
      help ? h('p', { class: 'mcb-help' }, help) : null);
  }
  function skipChip(name, label, note) {
    return { label: label, note: note, on: function () { return !on('skip_' + name); },
      flip: function (v) { setP('skip_' + name, !v); } };
  }

  // Database source: let the pipeline download it (and keep a copy or not), or point to a local copy.
  function dbSource(cfg) {
    var mode = cfg.manual ? 'local' : db(cfg.key);
    var body = [];
    if (!cfg.manual) {
      body.push(h('div', { class: 'mcb-seg mcb-seg-sm', role: 'radiogroup', 'aria-label': cfg.label + ' source' },
        [['auto', 'Download in the run'], ['local', 'I have it locally']].map(function (o) {
          var active = mode === o[0];
          return h('button', {
            type: 'button', role: 'radio', 'aria-checked': active ? 'true' : 'false',
            class: 'mcb-seg-btn' + (active ? ' is-on' : ''),
            onclick: function () { if (!active) { S.ui.db[cfg.key] = o[0]; commit(true); } }
          }, o[1]);
        })));
    }
    if (mode === 'local') {
      cfg.paths.forEach(function (pp) {
        body.push(textInput(pp.param, { label: pp.label, placeholder: pp.placeholder, required: true, help: pp.help }));
      });
    } else {
      if (cfg.save) {
        body.push(paramToggle(cfg.save, 'Keep a copy in ' + (cfg.saveTo || '<outdir>/Databases/'),
          'Pass that folder next time to skip the download.'));
      }
      if (cfg.autoExtra) body.push(cfg.autoExtra());
    }
    return h('div', { class: 'mcb-db' + (mode === 'local' ? ' is-local' : '') },
      h('div', { class: 'mcb-db-head' }, h('span', { class: 'mcb-db-name' }, cfg.label),
        cfg.size ? h('span', { class: 'mcb-db-size' }, cfg.size) : null),
      cfg.note ? h('p', { class: 'mcb-help' }, cfg.note) : null,
      h('div', { class: 'mcb-db-body' }, body));
  }

  function note(kind, content) { return h('div', { class: 'mcb-note is-' + kind, role: kind === 'error' ? 'alert' : null }, content); }
  function section(title, content) {
    return h('div', { class: 'mcb-section' }, title ? h('h3', { class: 'mcb-section-title' }, title) : null, content);
  }
  function link(href, text) { return h('a', { href: DOCS + href }, text); }

  /* ------------------------------------------------------------------ requirements */

  // What a choice pulls in: shown in the step where the choice is made and gathered on the last step.
  // Each item is [text, ok]; ok === null means automatic (set by the builder).
  function memGB() {
    var m = String(get('max_memory') || '').match(/^([\d.]+)\s*\.?\s*([KMGT]?)B$/);
    if (!m) return null;
    return Number(m[1]) * ({ K: 1 / 1048576, M: 1 / 1024, G: 1, T: 1024, '': 1 / 1073741824 }[m[2]]);
  }
  function dbItem(key, label, params) {
    if (db(key) !== 'local') return [label + ': downloaded in the run', null];
    return [label + ': local path', params.every(filled)];
  }
  function requirements(stepId) {
    var r = [];
    function need(title, items) { r.push({ title: title, items: items }); }
    var pp = get('readbased_postprocess');
    switch (stepId) {
      case 'input':
        if (fromAssemblies()) need('Your own assemblies', [['Reads given as a samplesheet (CSV)', get('input_format') === 'csv'],
          ['Assemblies samplesheet', filled('assembly_input')], ['Reads already cleaned (cleaning is skipped)', null]]);
        if (genomes() && S.ui.lineage === 'blank') need('Genomes without a lineage', [['A sourmash or GTDB-Tk database (Taxonomy of genomes step)',
          get('bin_classifier') === 'gtdbtk' ? filled('gtdbtk_db') : filled('sourmash_database') && filled('sourmash_taxonomy_csv')]]);
        break;
      case 'preprocess':
        if (!on('skip_preprocess') && S.ui.hostRef === 'local') need('Own host index', [['hostile index folder', filled('hostile_reference')], ['Index name matching the files', null]]);
        break;
      case 'profile':
        if (useMicrocosm()) need('sourmash', [['sourmash database', filled('sourmash_database')], ['Lineages file', filled('sourmash_taxonomy_csv')]]);
        if (yachtOn()) need('YACHT filtering', [['YACHT database (*_config.json)', filled('yacht_database')],
          ['A single sourmash database (YACHT takes one)', count(get('sourmash_database')) <= 1],
          ['Same k-mer size as sourmash (set together: k=' + get('sourmash_ksize') + ')', null]]);
        if (useSylph()) need('sylph', [['sylph database', filled('sylph_database')]]);
        if (useSylph() && on('enable_sylphtax')) need('sylph-tax', [['sylph-tax metadata', filled('sylph_taxonomy')],
          ['One metadata file per sylph database (' + count(get('sylph_database')) + ')', count(get('sylph_taxonomy')) === count(get('sylph_database'))]]);
        break;
      case 'readsamr':
        if (on('enable_rgi_bwt')) need('AMR genes in reads', [dbItem('card', 'CARD database', ['rgi_preparecarddb_dir'])]);
        break;
      case 'postprocess':
        if (pp) {
          var items = [];
          if (useMicrocosm()) items.push(['YACHT on, with one sourmash database', yachtRuns()]);
          if (pp === 'rgi_bwt') items.push(['AMR genes in reads (RGI bwt)', on('enable_rgi_bwt')]);
          if (pp === 'all' && useMicrocosm()) items.push(['AMR genes in reads (RGI bwt), only for the merged AMR table', on('enable_rgi_bwt') || null]);
          if (pp === 'taxburst' || pp === 'rgi_bwt') items.push(['sourmash + YACHT among the profilers', useMicrocosm()]);
          if ((pp === 'phyloseq' || pp === 'all') && useSylph()) items.push(['sylph-tax on, for the sylph lineages', on('enable_sylphtax')]);
          need('Post-processing: ' + pp, items);
        }
        break;
      case 'contigs':
        if (!on('skip_contig_analysis') && !on('skip_skani_search')) need('skani', [['skani database', filled('skani_db')]]);
        if (!on('skip_contig_analysis') && !on('skip_genomad_e2e')) need('geNomad', [dbItem('genomad', 'geNomad database', ['genomad_db'])]);
        break;
      case 'refine':
        if (!on('skip_binning_refinement') && get('refine_tool') === 'binette') need('Binette', [dbItem('checkm2', 'CheckM2 database', ['checkm2_db'])]);
        break;
      case 'binqc':
        if (qcTool('checkm2')) need('CheckM2', [dbItem('checkm2', 'CheckM2 database (shared with Binette)', ['checkm2_db'])]);
        if (qcTool('checkm')) need('CheckM', [dbItem('checkm', 'CheckM database', ['checkm_db'])]);
        if (qcTool('busco')) need('BUSCO', [dbItem('busco', 'BUSCO lineages', ['busco_db']), ['A lineage (auto picks one per bin)', filled('busco_db_lineage')]]);
        if (!on('skip_binqc') && on('run_gunc')) need('GUNC', [dbItem('gunc', 'GUNC database', ['gunc_db'])]);
        if (on('bin_domain_classification')) need('Domain classification', [['Runs on the bins kept after refinement (' + String(get('postbinning_input')).replace(/_/g, ' ') + ')', null]]);
        break;
      case 'annotation':
        if (bakta()) need('Bakta', [dbItem('bakta', 'Bakta database', ['annotation_bakta_db'])]);
        break;
      case 'classify':
        if (classifyOn() && get('bin_classifier') === 'gtdbtk') {
          var gb = memGB();
          need('GTDB-Tk', [['GTDB-Tk reference data, extracted', filled('gtdbtk_db')],
            ['At least 64 GB for the largest task (now ' + (get('max_memory') || '?') + ')', gb === null ? false : gb >= 64]]);
        } else if (classifyOn()) {
          need('sourmash', [['sourmash database', filled('sourmash_database')], ['Lineages file', filled('sourmash_taxonomy_csv')],
            ['k-mer size of the database (k=' + get('sourmash_ksize') + ')', null]]);
        }
        break;
      case 'bacterial':
        if (asm() && on('enable_bacterial_downstream')) need('Typing the bins', [['Binning on', !on('skip_binning')],
          ['Read depths for binning (--skip_binning_bamabund off)', !on('skip_binning_bamabund')],
          ['Bin classification on (the lineage picks the tools)', !on('skip_bin_classification')]]);
        if (downstreamOn()) {
          if (tool('genotyphi') || tool('paratype')) need('Salmonella lineage typing', [['SISTR on (its serovar routes Typhi to GenoTyphi and Paratyphi A to Paratype)', !on('skip_sistr')]]);
          if (tool('resfinder') && db('resfinder') === 'local') need('ResFinder, local', [['resfinder_db, pointfinder_db and disinfinder_db, all three', filled('resfinder_db') && filled('pointfinder_db') && filled('disinfinder_db')]]);
          if (tool('argnorm')) need('ARGnorm', [['At least one AMR tool besides RGI', ['amrfinderplus', 'abricate', 'resfinder', 'deeparg'].some(tool)]]);
        }
        break;
      case 'run':
        if (S.nf.executor !== 'local') need('Cluster', [['Queue or partition (unless your cluster has a default)', !!S.nf.queue], ['Nextflow itself stays running until the end (use -bg, screen or tmux)', null]]);
        if (S.nf.profile === 'singularity' || S.nf.profile === 'apptainer') need('Containers', [['Optional: an image cache folder, so images are not pulled again for every project', !!S.nf.cacheDir || null]]);
        break;
    }
    return r;
  }
  function renderRequirements(stepId, withTitle) {
    var reqs = requirements(stepId);
    if (!reqs.length) return null;
    return h('div', { class: 'mcb-reqs' }, reqs.map(function (q) {
      return h('div', { class: 'mcb-req-group' },
        h('p', { class: 'mcb-req-title' }, (withTitle ? BYID[stepId].title + ' · ' : '') + q.title + ' also needs'),
        h('ul', null, q.items.map(function (it) {
          var state = it[1] === null ? 'auto' : it[1] ? 'ok' : 'missing';
          return h('li', { class: 'is-' + state },
            h('span', { class: 'mcb-req-mark', 'aria-hidden': 'true' }),
            h('span', { class: 'mcb-offscreen' }, state === 'ok' ? 'Done: ' : state === 'missing' ? 'Missing: ' : 'Automatic: '),
            it[0]);
        })));
    }));
  }

  /* ------------------------------------------------------------------ steps */

  // Each step: when() decides if it is on the route; params lists what it may write (with conditions);
  // adv lists schema groups whose other parameters appear under "More options".
  function cond(name, test) { return { name: name, test: test }; }
  function always(name) { return { name: name, test: null }; }

  var STEPS = [
    {
      id: 'route', title: 'Goal', lead: 'What should this run produce?',
      when: function () { return true; },
      summary: function () { return { reads: 'Read-based profiling', assembly: 'Genomes from metagenomes', genomes: 'Typing my own genomes' }[S.route] || ''; },
      params: [cond('enable_readbase', reads)],
      render: renderRoute
    },
    {
      id: 'input', title: 'Input and output', lead: 'Where are the samples, and where should results go?',
      when: function () { return !!S.route; },
      summary: function () {
        if (genomes()) return get('genome_input') ? 'genome CSV set' : '';
        return get('input') ? (get('input_format') === 'directory' ? 'FASTQ folder' : 'samplesheet') : '';
      },
      params: [cond('input', function () { return !genomes(); }), cond('input_format', function () { return !genomes(); }),
        cond('assembly_input', fromAssemblies), cond('genome_input', genomes), always('outdir')],
      adv: ['input_output_options'],
      render: renderInput
    },
    {
      id: 'preprocess', title: 'Read cleaning', lead: 'Trim adapters with fastp and remove host reads with hostile.',
      when: function () { return (reads() || asm()) && !fromAssemblies(); },
      summary: function () { return on('skip_preprocess') ? 'off' : (S.ui.hostRef === 'local' ? 'own host index' : 'human host'); },
      params: [always('skip_preprocess'),
        cond('hostile_reference', function () { return !on('skip_preprocess') && S.ui.hostRef === 'local'; }),
        cond('hostile_index', function () { return !on('skip_preprocess') && S.ui.hostRef === 'local'; }),
        cond('hostile_ref_name', function () { return !on('skip_preprocess') && S.ui.hostRef === 'local'; })],
      adv: ['preprocessing_options'],
      render: renderPreprocess
    },
    {
      id: 'profile', title: 'Taxonomic profiling', lead: 'Who is in each sample, straight from the reads.',
      when: reads,
      summary: function () {
        return taxotools().map(function (t) { return t === 'microcosm' ? (yachtOn() ? 'sourmash + YACHT' : 'sourmash') : 'sylph'; }).join(' + ');
      },
      params: [always('readbased_taxotools'),
        cond('sourmash_database', useMicrocosm), cond('sourmash_taxonomy_csv', useMicrocosm),
        cond('sourmash_ksize', useMicrocosm), cond('skip_yacht', useMicrocosm), cond('yacht_database', yachtOn),
        cond('yacht_ksize', yachtOn), cond('enable_singlesketch', useMicrocosm),
        cond('sylph_database', useSylph), cond('enable_sylphtax', useSylph),
        cond('sylph_taxonomy', function () { return useSylph() && on('enable_sylphtax'); }),
        cond('sylph_datatype', useSylph)],
      advNames: ['sourmash_output_format', 'sourmash_thresholdbp', 'sourmash_profile', 'enable_copysketch',
        'enable_copysketch_singles', 'sylph_args'],
      render: renderProfile
    },
    {
      id: 'readsamr', title: 'AMR genes in reads', lead: 'Map reads to CARD with RGI bwt.',
      when: reads,
      summary: function () { return on('enable_rgi_bwt') ? 'RGI bwt' : 'off'; },
      params: [always('enable_rgi_bwt'),
        cond('rgi_preparecarddb_dir', function () { return on('enable_rgi_bwt') && db('card') === 'local'; }),
        cond('save_rgi_db', function () { return on('enable_rgi_bwt') && db('card') !== 'local'; }),
        cond('rgi_aligner', function () { return on('enable_rgi_bwt'); })],
      adv: ['rgi_options'],
      advTest: function () { return on('enable_rgi_bwt'); },
      render: renderReadsAmr
    },
    {
      id: 'postprocess', title: 'Tables for analysis', lead: 'Tables and charts that bring all samples of the run together.',
      when: reads,
      summary: function () { return get('readbased_postprocess') || 'none'; },
      params: [always('readbased_postprocess'),
        cond('readbased_postprocess_options', function () { return !!get('readbased_postprocess') && useMicrocosm(); })],
      render: renderPostprocess
    },
    {
      id: 'assembly', title: 'Assembly', lead: 'Assemble each sample into contigs.',
      when: function () { return asm() && !fromAssemblies(); },
      summary: function () {
        return [!on('skip_megahit') ? 'MEGAHIT' : null, !on('skip_spades') ? 'metaSPAdes' : null].filter(Boolean).join(' + ');
      },
      params: [always('skip_megahit'), always('skip_spades'), always('skip_quast')],
      adv: ['assembly_options'],
      render: renderAssembly
    },
    {
      id: 'contigs', title: 'Contig analysis', lead: 'Genes, viruses and plasmids, domains and close relatives of the contigs.',
      when: asm,
      summary: function () { return on('skip_contig_analysis') ? 'off' : 'on'; },
      params: [always('skip_contig_analysis'),
        cond('skip_contig_annotation', function () { return !on('skip_contig_analysis'); }),
        cond('annotation_tool', function () { return !on('skip_contig_analysis') && !on('skip_contig_annotation'); }),
        cond('skip_genomad_e2e', function () { return !on('skip_contig_analysis'); }),
        cond('genomad_db', function () { return !on('skip_contig_analysis') && !on('skip_genomad_e2e') && db('genomad') === 'local'; }),
        cond('save_genomad_db', function () { return !on('skip_contig_analysis') && !on('skip_genomad_e2e') && db('genomad') !== 'local'; }),
        cond('skip_tiara_contigqc', function () { return !on('skip_contig_analysis'); }),
        cond('skip_skani_search', function () { return !on('skip_contig_analysis'); }),
        cond('skani_db', function () { return !on('skip_contig_analysis') && !on('skip_skani_search'); })],
      adv: ['contig_analysis_options', 'pyrodigal_options', 'prokka_options'],
      advTest: function () { return !on('skip_contig_analysis'); },
      render: renderContigs
    },
    {
      id: 'binning', title: 'Binning', lead: 'Group contigs into genomes (bins).',
      when: asm,
      summary: function () {
        if (on('skip_binning')) return 'off';
        return BINNERS.filter(function (b) { return !on('skip_' + b[0]); }).length + ' binners';
      },
      params: [always('skip_binning'),
        cond('skip_metabat2', binningOn), cond('skip_maxbin2', binningOn), cond('skip_comebin', binningOn),
        cond('skip_semibin', binningOn), cond('skip_concoct', binningOn), cond('skip_vamb', binningOn),
        cond('semibin_environment', function () { return binningOn() && !on('skip_semibin'); }),
        cond('binprepare_mappertool', asm)],
      adv: ['binning_options'],
      advTest: binningOn,
      render: renderBinning
    },
    {
      id: 'refine', title: 'Bin refinement', lead: 'Combine the binners into one better set of bins.',
      when: binningOn,
      summary: function () { return on('skip_binning_refinement') ? 'off' : get('refine_tool'); },
      params: [always('skip_binning_refinement'),
        cond('refine_tool', function () { return !on('skip_binning_refinement'); }),
        cond('postbinning_input', function () { return !on('skip_binning_refinement'); }),
        cond('checkm2_db', function () { return needCheckm2() && db('checkm2') === 'local'; }),
        cond('save_checkm2_data', function () { return needCheckm2() && db('checkm2') !== 'local'; })],
      adv: ['binning_refinement_options'],
      advTest: function () { return !on('skip_binning_refinement'); },
      render: renderRefine
    },
    {
      id: 'binqc', title: 'Bin quality', lead: 'Completeness, contamination and chimeras of every bin.',
      when: binningOn,
      summary: function () { return on('skip_binqc') ? 'off' : binqcTools().join(', '); },
      params: [always('skip_binqc'),
        cond('binqc_tool', function () { return !on('skip_binqc'); }),
        cond('checkm_db', function () { return qcTool('checkm') && db('checkm') === 'local'; }),
        cond('save_checkm_data', function () { return qcTool('checkm') && db('checkm') !== 'local'; }),
        cond('busco_db', function () { return qcTool('busco') && db('busco') === 'local'; }),
        cond('save_busco_db', function () { return qcTool('busco') && db('busco') !== 'local'; }),
        cond('busco_db_lineage', function () { return qcTool('busco'); }),
        cond('run_gunc', function () { return !on('skip_binqc'); }),
        cond('gunc_db', function () { return !on('skip_binqc') && on('run_gunc') && db('gunc') === 'local'; }),
        cond('gunc_save_db', function () { return !on('skip_binqc') && on('run_gunc') && db('gunc') !== 'local'; }),
        cond('gunc_database_type', function () { return !on('skip_binqc') && on('run_gunc') && db('gunc') !== 'local'; }),
        always('bin_domain_classification')],
      adv: ['bin_qc_options', 'bin_domain_classification_options'],
      render: renderBinqc
    },
    {
      id: 'annotation', title: 'Bin annotation', lead: 'Annotate genes on every bin.',
      when: binningOn,
      summary: function () { return on('skip_bin_annotation') ? 'off' : get('bin_annotation_tool'); },
      params: [always('skip_bin_annotation'),
        cond('bin_annotation_tool', function () { return !on('skip_bin_annotation'); }),
        cond('annotation_bakta_db', function () { return bakta() && db('bakta') === 'local'; }),
        cond('save_bakta_db', function () { return bakta() && db('bakta') !== 'local'; }),
        cond('annotation_bakta_db_downloadtype', function () { return bakta() && db('bakta') !== 'local'; })],
      adv: ['bin_annotation_options'],
      advTest: function () { return !on('skip_bin_annotation'); },
      render: renderAnnotation
    },
    {
      id: 'classify', title: 'Taxonomy of genomes', lead: 'Give each genome a GTDB lineage.',
      when: function () { return binningOn() || genomes(); },
      summary: function () {
        if (genomes() && S.ui.lineage !== 'blank') return 'lineages given';
        if (asm() && on('skip_bin_classification')) return 'off';
        return get('bin_classifier') === 'gtdbtk' ? 'GTDB-Tk' : 'sourmash';
      },
      params: [cond('skip_bin_classification', binningOn),
        cond('bin_classifier', classifyOn),
        cond('sourmash_database', function () { return classifyOn() && get('bin_classifier') !== 'gtdbtk'; }),
        cond('sourmash_taxonomy_csv', function () { return classifyOn() && get('bin_classifier') !== 'gtdbtk'; }),
        cond('sourmash_ksize', function () { return classifyOn() && get('bin_classifier') !== 'gtdbtk'; }),
        cond('gtdbtk_db', function () { return classifyOn() && get('bin_classifier') === 'gtdbtk'; }),
        cond('gtdbtk_pplacer_cpus', function () { return classifyOn() && get('bin_classifier') === 'gtdbtk'; })],
      advNames: ['gtdbtk_args', 'gtdbtk_pplacer_scratch', 'sourmash_output_format', 'sourmash_thresholdbp',
        'sourmash_profile', 'enable_copysketch', 'enable_copysketch_singles'],
      advTest: classifyOn,
      render: renderClassify
    },
    {
      id: 'bacterial', title: 'Bacterial typing and AMR', lead: 'MLST, virulence, plasmids, species typing and AMR genes for every bacterial genome.',
      when: function () { return binningOn() || genomes(); },
      summary: function () {
        if (!downstreamOn()) return 'off';
        var all = BAC_TOOLS.length;
        return (all - BAC_TOOLS.filter(function (t) { return on('skip_' + t[0]); }).length) + ' of ' + all + ' tools';
      },
      params: [cond('enable_bacterial_downstream', asm),
        cond('skip_mlst', downstreamOn), cond('skip_virulencefinder', downstreamOn), cond('skip_plasmidfinder', downstreamOn),
        cond('skip_kleborate', downstreamOn), cond('skip_sistr', downstreamOn), cond('skip_seqsero2', downstreamOn),
        cond('skip_genotyphi', downstreamOn), cond('skip_paratype', downstreamOn), cond('skip_amrfinderplus', downstreamOn),
        cond('skip_abricate', downstreamOn), cond('skip_rgi', downstreamOn), cond('skip_resfinder', downstreamOn),
        cond('skip_deeparg', downstreamOn), cond('skip_argnorm', downstreamOn),
        cond('amrfinderplus_db', function () { return tool('amrfinderplus') && db('amrfinderplus') === 'local'; }),
        cond('save_amrfinderplus_db', function () { return tool('amrfinderplus') && db('amrfinderplus') !== 'local'; }),
        cond('rgi_preparecarddb_dir', function () { return tool('rgi') && db('card') === 'local'; }),
        cond('save_rgi_db', function () { return tool('rgi') && db('card') !== 'local'; }),
        cond('resfinder_db', function () { return tool('resfinder') && db('resfinder') === 'local'; }),
        cond('pointfinder_db', function () { return tool('resfinder') && db('resfinder') === 'local'; }),
        cond('disinfinder_db', function () { return tool('resfinder') && db('resfinder') === 'local'; }),
        cond('save_resfinder_db', function () { return tool('resfinder') && db('resfinder') !== 'local'; }),
        cond('deeparg_db', function () { return tool('deeparg') && db('deeparg') === 'local'; }),
        cond('save_deeparg_db', function () { return tool('deeparg') && db('deeparg') !== 'local'; }),
        cond('virulencefinder_db', function () { return tool('virulencefinder') && db('virulencefinder') === 'local'; }),
        cond('save_virulencefinder_db', function () { return tool('virulencefinder') && db('virulencefinder') !== 'local'; }),
        cond('paratype_bundle', function () { return paratype() && db('paratype') === 'local'; }),
        cond('save_paratype_db', function () { return paratype() && db('paratype') !== 'local'; }),
        cond('bacterial_decision_table', downstreamOn)],
      adv: ['database_options'],
      advExclude: ['sourmash_database', 'sourmash_taxonomy_csv', 'yacht_database'],
      advTest: downstreamOn,
      render: renderBacterial
    },
    {
      id: 'run', title: 'Where it runs', lead: 'Software, compute and the Nextflow options of this run.',
      when: function () { return !!S.route; },
      summary: function () { return S.nf.profile + (S.nf.executor !== 'local' ? ' · ' + S.nf.executor : ''); },
      params: [always('max_cpus'), always('max_memory'), always('max_time'), always('skip_multiqc'), always('email')],
      adv: ['misc_options'],
      render: renderRun
    },
    {
      id: 'done', title: 'Run it', lead: 'Copy the command, or save a file and run it.',
      when: function () { return !!S.route; },
      summary: function () { var n = problems().filter(function (x) { return x.level === 'error'; }).length; return n ? n + ' to fix' : 'ready'; },
      params: [],
      render: renderDone
    }
  ];
  var BYID = {};
  STEPS.forEach(function (s) { BYID[s.id] = s; });

  var BINNERS = [['metabat2', 'MetaBAT2', 'robust default'], ['maxbin2', 'MaxBin2', 'marker genes'], ['comebin', 'COMEBin', 'deep learning; faster with a GPU'],
    ['semibin', 'SemiBin2', 'trained per environment'], ['concoct', 'CONCOCT', 'slow on large data'], ['vamb', 'VAMB', 'best with many samples']];
  var BAC_TOOLS = [
    ['mlst', 'MLST', 'sequence type (ST), to compare strains'], ['virulencefinder', 'VirulenceFinder', 'virulence genes'], ['plasmidfinder', 'PlasmidFinder', 'plasmid types'],
    ['kleborate', 'Kleborate', 'Klebsiella and E. coli: ST, virulence, AMR'], ['sistr', 'SISTR', 'Salmonella serovar'],
    ['seqsero2', 'SeqSero2', 'Salmonella serotype, second opinion'], ['genotyphi', 'GenoTyphi', 'S. Typhi genotype'], ['paratype', 'Paratype', 'S. Paratyphi A genotype'],
    ['amrfinderplus', 'AMRFinderPlus', 'NCBI reference; genes and point mutations'], ['abricate', 'ABRicate', 'quick screen'], ['rgi', 'RGI', 'CARD resistance models'],
    ['resfinder', 'ResFinder', 'acquired genes, mutations, disinfectants'], ['deeparg', 'DeepARG', 'finds distant, novel genes'], ['argnorm', 'ARGnorm', 'puts all results on one drug-class scale']];
  function qcTool(t) { return binningOn() && !on('skip_binqc') && binqcTools().indexOf(t) >= 0; }
  function bakta() { return binningOn() && !on('skip_bin_annotation') && get('bin_annotation_tool') === 'bakta'; }
  function paratype() { return tool('paratype') && !on('skip_sistr'); }

  /* advanced parameter lists, from the schema minus what the steps already ask */
  function buildAdvanced() {
    var curated = {};
    STEPS.forEach(function (s) { s.params.forEach(function (e) { curated[e.name] = true; }); });
    var placed = {};
    STEPS.forEach(function (s) {
      var names = [];
      (s.adv || []).forEach(function (g) { names = names.concat(GROUPS[g] || []); });
      names = names.concat(s.advNames || []);
      s.advList = names.filter(function (n, i) {
        return P[n] && names.indexOf(n) === i && !curated[n] && HIDDEN.indexOf(n) < 0 &&
          (s.advExclude || []).indexOf(n) < 0 && P[n].type !== 'array';
      });
      s.advList.forEach(function (n) { placed[n] = true; });
    });
    var missing = Object.keys(P).filter(function (n) { return !curated[n] && !placed[n] && HIDDEN.indexOf(n) < 0; });
    if (missing.length && window.console) console.info('[microcosm builder] parameters not offered:', missing.join(', '));
  }

  /* ------------------------------------------------------------------ step views */

  function renderRoute() {
    var cards = [
      { id: 'reads', title: 'Profile my reads', text: 'Who is in each sample and how much, straight from the reads. sourmash + YACHT or sylph, and AMR genes with RGI.', meta: 'Fastest · reads in' },
      { id: 'assembly', title: 'Recover genomes from metagenomes', text: 'Assemble, bin, refine and check genomes (MAGs), then type the bacterial ones if you want.', meta: 'Most complete · reads in' },
      { id: 'genomes', title: 'Type my own genomes', text: 'Isolate assemblies or MAGs made elsewhere: MLST, AMR, Kleborate, Salmonella typing.', meta: 'FASTA in' }
    ];
    return h('div', { class: 'mcb-routes', role: 'radiogroup', 'aria-label': 'Goal of the run' }, cards.map(function (c, i) {
      var active = S.route === c.id;
      return h('button', {
        type: 'button', role: 'radio', 'aria-checked': active ? 'true' : 'false',
        class: 'mcb-route is-' + c.id + (active ? ' is-on' : ''), style: '--i:' + i,
        onclick: function () { pickRoute(c.id); }
      }, routeGlyph(c.id),
      h('span', { class: 'mcb-route-title' }, c.title),
      h('span', { class: 'mcb-route-text' }, c.text),
      h('span', { class: 'mcb-route-meta' }, c.meta));
    }));
  }
  var autoAdvance = null;
  function pickRoute(id) {
    var changed = S.route !== id;
    S.route = id;
    setP('enable_readbase', id === 'reads');
    if (changed) S.seen = { route: true };
    commit(true);
    if (changed) autoAdvance = setTimeout(function () { go(nextStep('route')); }, REDUCED ? 0 : 220);
  }
  // Small line drawings of each route, in the metro-map style of the pipeline diagram.
  function routeGlyph(id) {
    var paths = {
      reads: '<path d="M6 20h44"/><circle cx="12" cy="20" r="4"/><circle cx="28" cy="20" r="4"/><circle cx="44" cy="20" r="4"/><path d="M28 24v10h16" class="g2"/><circle cx="44" cy="34" r="3" class="g2"/>',
      assembly: '<path d="M6 12h14l8 8h22"/><path d="M6 28h14l8-8" /><circle cx="12" cy="12" r="3.5"/><circle cx="12" cy="28" r="3.5"/><circle cx="36" cy="20" r="4.5"/><circle cx="48" cy="20" r="3" class="g2"/>',
      genomes: '<circle cx="14" cy="20" r="8"/><circle cx="14" cy="20" r="3"/><path d="M22 20h26"/><circle cx="34" cy="20" r="3.5"/><circle cx="46" cy="20" r="3.5" class="g2"/>'
    };
    var span = h('span', { class: 'mcb-route-glyph', 'aria-hidden': 'true' });
    span.innerHTML = '<svg viewBox="0 0 56 40" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round">' + paths[id] + '</svg>';
    return span;
  }

  function templateLink(name, text, content) {
    return h('button', {
      type: 'button', class: 'mcb-linkbtn', onclick: function () { download(name, content, 'text/csv'); }
    }, text);
  }

  function renderInput() {
    var out = [];
    if (genomes()) {
      out.push(textInput('genome_input', {
        label: 'Genome samplesheet', placeholder: 'genomes.csv', required: true, wide: true,
        help: 'CSV with the columns sample_id,fasta[,lineage]; FASTA paths may be relative to the CSV.'
      }));
      out.push(h('p', { class: 'mcb-inline' }, templateLink('genomes.csv',
        'Download a genomes.csv template', 'sample_id,fasta,lineage\ngenome_1,genome_1.fasta.gz,\ngenome_2,genome_2.fna,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n')));
      out.push(segmented('Is the lineage column filled?', [
        { value: 'filled', label: 'Yes, for every genome' }, { value: 'blank', label: 'No, for some or all' }
      ], S.ui.lineage, function (v) { S.ui.lineage = v; }, {
        wide: true, help: 'Genomes with a GTDB lineage are not reclassified. Blank rows need a sourmash or GTDB-Tk database (next steps).'
      }));
    } else {
      if (asm()) {
        out.push(segmented('Start from', [{ value: 'reads', label: 'Raw reads' }, { value: 'assemblies', label: 'Reads + my own assemblies' }],
          S.ui.asmStart, function (v) {
            S.ui.asmStart = v;
            if (v === 'assemblies') setP('input_format', 'csv');
          }, { wide: true, help: 'With your own assemblies, read cleaning and assembly are skipped; the reads are only mapped back for binning, so give clean reads.' }));
      }
      out.push(segmented('Reads are given as', [
        { value: 'csv', label: 'A samplesheet (CSV)' },
        { value: 'directory', label: 'A folder of FASTQ files', disabled: fromAssemblies(), why: 'Own assemblies need a samplesheet' }
      ], get('input_format'), function (v) { setP('input_format', v); }, { param: 'input_format', wide: true }));
      if (get('input_format') === 'directory') {
        out.push(textInput('input', { label: 'FASTQ folder', placeholder: '/data/run_01/fastq', required: true, wide: true,
          help: 'Paired files named like sample_1.fq.gz / sample_2.fq.gz or sample_R1.fastq.gz / sample_R2.fastq.gz.' }));
      } else {
        out.push(textInput('input', { label: 'Samplesheet', placeholder: 'samples.csv', required: true, wide: true,
          help: 'One row per sample and run; runs of one sample are merged after cleaning.' }));
        out.push(h('p', { class: 'mcb-inline' }, templateLink('samples.csv', 'Download a samples.csv template',
          'sample_id,run_id,group,short_reads_1,short_reads_2,long_reads\nsample1,1,0,data/sample1_1.fq.gz,data/sample1_2.fq.gz,\nsample2,1,0,data/sample2_1.fq.gz,data/sample2_2.fq.gz,\n')));
      }
      if (fromAssemblies()) {
        out.push(textInput('assembly_input', { label: 'Assemblies samplesheet', placeholder: 'assemblies.csv', required: true, wide: true,
          help: 'CSV with sample_id, assembler and fasta columns, one row per assembly.' }));
      }
    }
    out.push(textInput('outdir', { label: 'Results folder', placeholder: './results', wide: true,
      help: 'Created if missing. Relative paths start from the folder you launch Nextflow in.' }));
    return h('div', { class: 'mcb-grid' }, out);
  }

  function renderPreprocess() {
    var out = [paramToggle('skip_preprocess', 'Remove adapters, poor-quality reads and host DNA?',
      'Recommended for raw reads: leftover adapters and human reads end up in the profiles and assemblies, and human reads should not be shared. Turn off only for reads you already cleaned. Uses fastp, then hostile.', true)];
    if (!on('skip_preprocess')) {
      out.push(segmented('Whose DNA should be removed?', [
        { value: 'auto', label: 'Human (downloaded automatically)' }, { value: 'local', label: 'Another host: I have an index' }
      ], S.ui.hostRef, function (v) { S.ui.hostRef = v; }, { wide: true, help: 'For animal, plant or other hosts, build a Bowtie2 index of the host genome first (see Databases).' }));
      if (S.ui.hostRef === 'local') {
        out.push(h('div', { class: 'mcb-grid' },
          textInput('hostile_reference', { label: 'Index folder', placeholder: '/db/hostile', required: true, help: 'Folder holding the Bowtie2 index files.' }),
          textInput('hostile_index', { label: 'Index name', placeholder: 'human-t2t-hla', help: 'Prefix of the index files in that folder.' })));
      }
    }
    return out;
  }

  function renderProfile() {
    var tools = taxotools();
    function flip(t, v) {
      var next = ['microcosm', 'sylph'].filter(function (x) { return x === t ? v : tools.indexOf(x) >= 0; });
      setP('readbased_taxotools', next.join(','));
    }
    var out = [chips('How should the organisms be identified?', [
      { label: 'sourmash + YACHT', note: 'default; checks each species has enough evidence', on: function () { return tools.indexOf('microcosm') >= 0; }, flip: function (v) { flip('microcosm', v); } },
      { label: 'sylph', note: 'very fast, good for large batches; its own database', on: function () { return tools.indexOf('sylph') >= 0; }, flip: function (v) { flip('sylph', v); } }
    ], 'Pick one, or both to compare them; each writes its own results. Both report relative abundance per sample.')];
    if (!tools.length) out.push(note('error', 'Choose at least one profiler.'));
    if (useMicrocosm()) {
      var k = String(get('sourmash_ksize'));
      out.push(section('sourmash', [
        h('div', { class: 'mcb-grid' },
          textInput('sourmash_database', { label: 'sourmash database', placeholder: '/db/gtdb-rs226-reps.k31.rocksdb', required: true, wide: true,
            help: 'A .rocksdb or .zip; several zip databases separated by commas (this turns YACHT off).' }),
          textInput('sourmash_taxonomy_csv', { label: 'Lineages file', placeholder: '/db/gtdb-rs226.lineages.csv', required: true, wide: true }),
          segmented('How fine should the identification be?', [{ value: '31', label: 'Species (k=31)' }, { value: '51', label: 'Strain (k=51)' }, { value: '21', label: 'Genus-level, sensitive (k=21)' }],
            k, function (v) { setP('sourmash_ksize', Number(v)); setP('yacht_ksize', Number(v)); },
            { param: 'sourmash_ksize', help: 'Must match the k-mer size the databases were built with; the database file name usually says (k31, k51).' })),
        paramToggle('skip_yacht', 'Drop organisms with too little evidence?',
          'Recommended: without it, organisms that share only a few k-mers with your reads can appear as present. Uses YACHT, which needs a database trained on the same sourmash database.', true),
        yachtOn() ? textInput('yacht_database', { label: 'YACHT database (JSON)', placeholder: '/db/gtdb-rs226-reps.k31_0.995_config.json', required: true, wide: true,
          help: 'The *_config.json of a trained YACHT database; edit its paths after download.' }) : null,
        paramToggle('enable_singlesketch', 'Process every sample on its own?', 'Off: all samples are compared to the database in one batch (fast). On: one search per sample; slower, but easier on memory with a few very large samples.')
      ]));
    }
    if (useSylph()) {
      out.push(section('sylph', [
        textInput('sylph_database', { label: 'sylph database', placeholder: '/db/gtdb-r226-c200-dbv2.syl2db', required: true, wide: true,
          help: 'A .syldb or .syl2db; several separated by commas.' }),
        paramToggle('enable_sylphtax', 'Add the full lineage (phylum to species) to sylph results?', 'sylph alone reports reference genomes; sylph-tax turns them into taxa. Needed for phyloseq tables from sylph.'),
        on('enable_sylphtax') ? textInput('sylph_taxonomy', { label: 'sylph-tax metadata', placeholder: '/db/gtdb_r226_metadata.tsv.gz', required: true, wide: true,
          help: 'One file per sylph database, in the same order, separated by commas.' }) : null,
        paramSeg('sylph_datatype', 'Which abundance should the tables show?', [{ value: 'taxonomic_abundance', label: 'Share of cells (taxonomic)' }, { value: 'sequence_abundance', label: 'Share of reads (sequence)' }])
      ]));
    }
    out.push(h('p', { class: 'mcb-help' }, 'No database yet? See ', link('getting-started/databases/', 'Databases'), ' for download links.'));
    return out;
  }

  function cardDb() {
    return dbSource({ key: 'card', label: 'CARD (RGI)', size: '~12 GB', save: 'save_rgi_db', saveTo: '<outdir>/Databases/RGI',
      paths: [{ param: 'rgi_preparecarddb_dir', label: 'Prepared CARD folder', placeholder: '/db/rgi_db', help: 'A folder prepared by RGI (rgi load), e.g. Databases/RGI from an earlier run.' }] });
  }

  function renderReadsAmr() {
    var out = [paramToggle('enable_rgi_bwt', 'Look for antimicrobial resistance genes in the reads?',
      'Gives the resistance genes of each sample and how abundant they are, without assembling, so it also catches genes of rare organisms. It does not tell which organism carries a gene. Uses RGI bwt and the CARD database.')];
    if (on('enable_rgi_bwt')) {
      out.push(cardDb());
      out.push(paramSeg('rgi_aligner', 'Read aligner', [{ value: 'kma', label: 'KMA (recommended)' }, { value: 'bowtie2', label: 'Bowtie2' }, { value: 'bwa', label: 'BWA' }]));
    }
    return out;
  }

  // What each post-processing mode gives with the profilers chosen, or why it cannot run (rules of workflow/microcosm.nf).
  function postprocessInfo(mode) {
    var mic = useMicrocosm(), syl = useSylph(), rgi = on('enable_rgi_bwt'), tax = on('enable_sylphtax');
    var toProfile = { step: 'profile', label: 'Go to Taxonomic profiling' };
    if (!mode) return { gives: 'Only the per-sample results.' };
    if (mic && !yachtRuns()) return { why: 'Post-processing reads the sourmash + YACHT result, so YACHT must be on with a single sourmash database.', fix: toProfile };
    if (mode === 'phyloseq') {
      if (syl && !tax) return { why: 'sylph results need lineages: turn on sylph-tax.', fix: toProfile };
      return { gives: 'Abundance, taxonomy and sample tables to load in R (phyloseq), from ' + [mic ? 'sourmash + YACHT' : null, syl ? 'sylph-tax' : null].filter(Boolean).join(' and ') + '.' };
    }
    if (mode === 'taxburst') {
      if (!mic) return { why: 'taxburst charts are drawn from the sourmash + YACHT result; with sylph only there is nothing to draw. Add sourmash + YACHT as a profiler.', fix: toProfile };
      return { gives: 'Interactive sunburst charts to explore in a browser, from sourmash + YACHT' + (syl ? ' (sylph results are not drawn).' : '.') };
    }
    if (mode === 'rgi_bwt') {
      if (!rgi) return { why: 'Needs the AMR genes in reads step (RGI bwt) to be on.', fix: { step: 'readsamr', label: 'Go to AMR genes in reads' } };
      if (!mic) return { why: 'RGI bwt still runs for each sample, but the merged AMR table is built with the sourmash + YACHT result. Add sourmash + YACHT as a profiler.', fix: toProfile };
      return { gives: 'One table of RGI bwt results across all samples.' };
    }
    if (syl && !tax) return { why: 'sylph results need lineages for the phyloseq tables: turn on sylph-tax.', fix: toProfile };
    var parts = ['phyloseq tables (' + [mic ? 'sourmash + YACHT' : null, syl ? 'sylph-tax' : null].filter(Boolean).join(', ') + ')'];
    if (mic) parts.push('taxburst charts');
    if (mic && rgi) parts.push('the merged AMR table');
    var missing = [];
    if (!mic) missing.push('taxburst and the merged AMR table need sourmash + YACHT');
    else if (!rgi) missing.push('no AMR table, since RGI bwt is off');
    return { gives: 'Everything that applies: ' + parts.join(', ') + '.' + (missing.length ? ' Not made: ' + missing.join('; ') + '.' : ''), same: !mic ? 'phyloseq' : null };
  }

  function renderPostprocess() {
    var current = get('readbased_postprocess') || '';
    var modes = [['', 'No'], ['phyloseq', 'Tables for R (phyloseq)'], ['taxburst', 'Interactive charts (taxburst)'],
      ['rgi_bwt', 'Merged AMR gene table'], ['all', 'All that apply']];
    var list = h('div', { class: 'mcb-opts', role: 'radiogroup', 'aria-label': 'Post-processing' });
    modes.forEach(function (m) {
      var info = postprocessInfo(m[0]), active = current === m[0], blocked = !!info.why;
      var id = nextId('opt');
      list.appendChild(h('div', { class: 'mcb-opt' + (active ? ' is-on' : '') + (blocked ? ' is-blocked' : '') },
        h('button', {
          type: 'button', role: 'radio', id: id, class: 'mcb-opt-btn', 'aria-checked': active ? 'true' : 'false',
          'aria-disabled': blocked && !active ? 'true' : null, 'aria-describedby': id + '-d',
          onclick: function () { if (!active && !blocked) { setP('readbased_postprocess', m[0]); commit(true); } }
        }, h('span', { class: 'mcb-opt-radio', 'aria-hidden': 'true' }),
        h('span', { class: 'mcb-opt-text' },
          h('span', { class: 'mcb-opt-title' }, m[1], m[0] ? h('span', { class: 'mcb-flag' }, m[0]) : null),
          h('span', { class: 'mcb-opt-desc', id: id + '-d' },
            blocked ? h('span', { class: 'mcb-opt-why' }, (active ? 'Will fail: ' : 'Not available: ') + info.why) : info.gives))),
        blocked && info.fix ? h('button', { type: 'button', class: 'mcb-linkbtn mcb-opt-fix', onclick: function () { go(info.fix.step); } }, info.fix.label) : null));
    });
    var out = [h('div', { class: 'mcb-field is-wide' },
      h('span', { class: 'mcb-label' }, 'Do you want tables and charts that merge all samples?', h('span', { class: 'mcb-flag' }, '--readbased_postprocess')),
      list)];
    var cur = postprocessInfo(current);
    if (current === 'all' && cur.same) out.push(note('info', 'With sylph only, "All that apply" gives the same files as phyloseq.'));
    if (current && useMicrocosm()) {
      out.push(textInput('readbased_postprocess_options', { label: 'Post-processing options', wide: true,
        help: 'Passed to the post-processing script as typed, e.g. --min_coverage 0.1 --remove_unclassified.' +
          (useSylph() ? ' They apply to the sourmash + YACHT tables; the sylph tables have no options.' : '') }));
    }
    return out;
  }

  function renderAssembly() {
    var out = [chips('Which assemblers should rebuild the genomes from the reads?', [
      { label: 'MEGAHIT', note: 'fast, low memory; a safe default', on: function () { return !on('skip_megahit'); }, flip: function (v) { setP('skip_megahit', !v); } },
      { label: 'metaSPAdes', note: 'often longer contigs; needs much more memory', on: function () { return !on('skip_spades'); }, flip: function (v) { setP('skip_spades', !v); } }
    ], 'Each assembler gives its own contigs and genomes, so you can compare them. On a small machine keep MEGAHIT only.')];
    if (on('skip_megahit') && on('skip_spades')) out.push(note('error', 'Keep at least one assembler.'));
    out.push(paramToggle('skip_quast', 'Report how good each assembly is?', 'Contig lengths, N50 and totals, to compare samples and assemblers. Uses QUAST.', true));
    return out;
  }

  function renderContigs() {
    var out = [paramToggle('skip_contig_analysis', 'Look at the contigs themselves, before they are grouped into genomes?',
      'Useful for what binning misses: viruses, plasmids and small genomes rarely end up in a bin.', true)];
    if (on('skip_contig_analysis')) return out;
    out.push(paramToggle('skip_contig_annotation', 'Find the genes on every contig?', 'Gene counts and protein sequences for all contigs.', true));
    if (!on('skip_contig_annotation')) {
      out.push(paramSeg('annotation_tool', 'Gene finder', [{ value: 'pyrodigal', label: 'Pyrodigal (fast, genes only)' }, { value: 'prokka', label: 'Prokka (genes and their functions)' }]));
    }
    out.push(paramToggle('skip_genomad_e2e', 'Find viruses and plasmids among the contigs?', 'Phages and plasmids often carry resistance and virulence genes. Uses geNomad.', true));
    if (!on('skip_genomad_e2e')) {
      out.push(dbSource({ key: 'genomad', label: 'geNomad database', size: '~2 GB', save: 'save_genomad_db',
        paths: [{ param: 'genomad_db', label: 'geNomad database folder', placeholder: '/db/genomad_db' }] }));
    }
    out.push(paramToggle('skip_tiara_contigqc', 'Tell bacterial, archaeal and eukaryotic contigs apart?', 'Shows how much of the sample is not bacterial (fungi, protists, host leftovers). Uses Tiara.', true));
    out.push(paramToggle('skip_skani_search', 'Find the closest known genome for each contig?', 'Average nucleotide identity to reference genomes. Needs a skani database you build or download (about 50 GB). Uses skani.', true));
    if (!on('skip_skani_search')) {
      out.push(textInput('skani_db', { label: 'skani database', placeholder: '/db/skani/gtdb_r226', required: true, wide: true }));
    }
    return out;
  }

  function renderBinning() {
    var out = [];
    var locked = asm() && on('enable_bacterial_downstream');
    out.push(toggle('Group the contigs into genomes (MAGs)?', locked ? 'Needed for bacterial typing (turned on in a later step).' : 'Needed for every genome-level result: quality, taxonomy, annotation and typing.', !on('skip_binning'),
      function (v) { setP('skip_binning', !v); }, { param: '--skip_binning', disabled: locked }));
    if (on('skip_binning')) return out;
    out.push(chips('Which binning tools?', BINNERS.map(function (b) { return skipChip(b[0], b[1], b[2]); }),
      'Each tool groups contigs differently. Several tools plus refinement (next step) usually recover more and cleaner genomes than one tool, at the cost of time.'));
    if (BINNERS.every(function (b) { return on('skip_' + b[0]); })) out.push(note('error', 'Keep at least one binner.'));
    if (!on('skip_semibin')) {
      out.push(textInput('semibin_environment', { label: 'What kind of sample is it? (SemiBin2 model)', help: 'A SemiBin2 pre-trained model: human_gut, dog_gut, ocean, soil, cat_gut, human_oral, mouse_gut, pig_gut, built_environment, wastewater, chicken_caecum or global.' }));
    }
    out.push(paramSeg('binprepare_mappertool', 'Read mapper (measures how deeply each contig is covered)', [{ value: 'bowtie2', label: 'Bowtie2' }, { value: 'bwamem2', label: 'BWA-MEM2' }],
      'Coverage is what lets binners group contigs of the same organism. Reads are mapped once and shared by all steps.'));
    return out;
  }

  function checkm2Db() {
    return dbSource({ key: 'checkm2', label: 'CheckM2 database', size: '~3 GB', save: 'save_checkm2_data',
      paths: [{ param: 'checkm2_db', label: 'CheckM2 database file', placeholder: '/db/checkm2/uniref100.KO.1.dmnd' }] });
  }

  function renderRefine() {
    var out = [paramToggle('skip_binning_refinement', 'Merge the results of the binning tools into one best set of genomes?',
      'Recommended with more than one binner: keeps, for each genome, the most complete and least contaminated version.', true)];
    if (on('skip_binning_refinement')) return out;
    out.push(paramSeg('refine_tool', 'Refinement tool', [{ value: 'binette', label: 'Binette (recommended)' }, { value: 'dastool', label: 'DAS Tool' }],
      'Binette judges genomes with CheckM2; DAS Tool with single-copy genes and can give no result on small data.'));
    out.push(paramSeg('postbinning_input', 'Which genomes go on to quality, annotation and typing?', [
      { value: 'both', label: 'Raw and refined' }, { value: 'refined_bins_only', label: 'Refined only' }, { value: 'raw_bins_only', label: 'Raw only' }]));
    if (get('refine_tool') === 'binette') out.push(checkm2Db());
    return out;
  }

  function renderBinqc() {
    var out = [paramToggle('skip_binqc', 'Check how complete and clean each genome is?',
      'Tells you which genomes are good enough to report; MAG standards ask for completeness and contamination (e.g. over 90% and under 5% for high quality).', true)];
    if (!on('skip_binqc')) {
      var tools = binqcTools();
      out.push(chips('How should quality be measured?', [['checkm2', 'CheckM2', 'fast; recommended for bacteria and archaea'], ['checkm', 'CheckM', 'marker genes; the classic reference'], ['busco', 'BUSCO', 'also for eukaryotes; large database']].map(function (t) {
        return { label: t[1], note: t[2], on: function () { return tools.indexOf(t[0]) >= 0; },
          flip: function (v) {
            var next = ['busco', 'checkm', 'checkm2'].filter(function (x) { return x === t[0] ? v : tools.indexOf(x) >= 0; });
            setP('binqc_tool', next.join(','));
          } };
      })));
      if (!tools.length) out.push(note('error', 'Choose at least one quality tool, or turn bin quality off.'));
      if (tools.indexOf('checkm2') >= 0 && !(get('refine_tool') === 'binette' && !on('skip_binning_refinement'))) out.push(checkm2Db());
      else if (tools.indexOf('checkm2') >= 0) out.push(h('p', { class: 'mcb-help' }, 'CheckM2 uses the database set in the refinement step.'));
      if (tools.indexOf('checkm') >= 0) {
        out.push(dbSource({ key: 'checkm', label: 'CheckM database', size: '~1.5 GB', save: 'save_checkm_data',
          paths: [{ param: 'checkm_db', label: 'CheckM data folder', placeholder: '/db/checkm' }] }));
      }
      if (tools.indexOf('busco') >= 0) {
        out.push(dbSource({ key: 'busco', label: 'BUSCO lineages', size: 'up to ~60 GB', save: 'save_busco_db',
          paths: [{ param: 'busco_db', label: 'BUSCO downloads folder', placeholder: '/db/busco_downloads' }] }));
        out.push(textInput('busco_db_lineage', { label: 'BUSCO lineage', help: 'auto picks the lineage per bin; or a dataset name such as bacteria_odb10.' }));
      }
      out.push(paramToggle('run_gunc', 'Check whether genomes are chimeras (pieces of different organisms)?',
        'Completeness tools can miss contamination from distant organisms. Worth it before you publish or deposit MAGs. Uses GUNC (about 13 GB database).'));
      if (on('run_gunc')) {
        out.push(dbSource({ key: 'gunc', label: 'GUNC database', size: '~13 GB', save: 'gunc_save_db',
          paths: [{ param: 'gunc_db', label: 'GUNC database file', placeholder: '/db/gunc/gunc_db_progenomes2.1.dmnd' }],
          autoExtra: function () {
            return paramSeg('gunc_database_type', 'Which reference should GUNC compare against?',
              [{ value: 'progenomes', label: 'proGenomes (default)' }, { value: 'gtdb', label: 'GTDB' }],
              'GTDB matches the taxonomy used everywhere else in the pipeline.');
          } }));
      }
    }
    out.push(paramToggle('bin_domain_classification', 'Sort genomes into bacteria, archaea and eukaryotes?', 'Useful for environmental or food samples that hold fungi or protists: their bins are labelled before annotation. Uses Tiara.'));
    return out;
  }

  function renderAnnotation() {
    var out = [paramToggle('skip_bin_annotation', 'Name the genes of each genome and their functions?', 'Needed for functional questions (pathways, enzymes) and for submitting genomes.', true)];
    if (on('skip_bin_annotation')) return out;
    out.push(paramSeg('bin_annotation_tool', 'Annotation tool', [{ value: 'prokka', label: 'Prokka (quick, no download)' }, { value: 'bakta', label: 'Bakta (richer, large database)' }],
      'Bakta gives standardised gene names and database links; its database is 4 GB (light) to 85 GB (full).'));
    if (get('bin_annotation_tool') === 'bakta') {
      out.push(dbSource({ key: 'bakta', label: 'Bakta database', size: '4 GB (light) to 85 GB (full)', save: 'save_bakta_db',
        paths: [{ param: 'annotation_bakta_db', label: 'Bakta database folder', placeholder: '/db/bakta/db' }],
        autoExtra: function () { return paramSeg('annotation_bakta_db_downloadtype', 'Database to download', [{ value: 'full', label: 'Full' }, { value: 'light', label: 'Light' }]); } }));
    }
    return out;
  }

  function renderClassify() {
    var out = [];
    if (genomes() && S.ui.lineage !== 'blank') {
      return [note('info', 'Every genome has a lineage in the samplesheet, so no classification runs. Change this in Input and output.')];
    }
    if (asm()) {
      var locked = on('enable_bacterial_downstream');
      out.push(toggle('Name each genome (which species is it)?', locked ? 'Needed for bacterial typing (turned on in the next step).' : 'Gives a GTDB lineage per genome; needed for typing and for the species tables.',
        !on('skip_bin_classification'), function (v) { setP('skip_bin_classification', !v); },
        { param: '--skip_bin_classification', disabled: locked }));
      if (on('skip_bin_classification')) return out;
    }
    out.push(paramSeg('bin_classifier', 'How should genomes be named?', [{ value: 'sourmash', label: 'sourmash (fast, light)' }, { value: 'gtdbtk', label: 'GTDB-Tk (reference method, heavy)' }],
      'sourmash compares each genome with a sourmash database and is enough for most uses. GTDB-Tk places genomes in the GTDB tree, the standard for publications; it needs over 100 GB of reference data and 64-128 GB of memory.'));
    if (get('bin_classifier') === 'gtdbtk') {
      out.push(h('div', { class: 'mcb-grid' },
        textInput('gtdbtk_db', { label: 'GTDB-Tk reference data', placeholder: '/db/gtdbtk/release226', required: true, wide: true,
          help: 'The extracted folder; the pipeline does not download it.' }),
        numberInput('gtdbtk_pplacer_cpus', { label: 'pplacer CPUs', min: 1, step: 1, help: 'Memory grows with each CPU.' })));
    } else {
      out.push(h('div', { class: 'mcb-grid' },
        textInput('sourmash_database', { label: 'sourmash database', placeholder: '/db/gtdb-rs226-reps.k31.rocksdb', required: true, wide: true }),
        textInput('sourmash_taxonomy_csv', { label: 'Lineages file', placeholder: '/db/gtdb-rs226.lineages.csv', required: true, wide: true }),
        segmented('k-mer size', [{ value: '21', label: 'k=21' }, { value: '31', label: 'k=31' }, { value: '51', label: 'k=51' }],
          String(get('sourmash_ksize')), function (v) { setP('sourmash_ksize', Number(v)); },
          { param: 'sourmash_ksize', help: 'Must match the database.' })));
    }
    return out;
  }

  function renderBacterial() {
    var out = [];
    if (asm()) {
      out.push(toggle('Type the bacterial genomes and screen them for antimicrobial resistance?',
        'For surveillance and outbreak work: sequence type, virulence and plasmids for every bacterium, extra typing for Klebsiella, E. coli and Salmonella, and AMR genes. Turns binning and genome naming on.',
        on('enable_bacterial_downstream'), function (v) {
          setP('enable_bacterial_downstream', v);
          if (v) { setP('skip_binning', false); setP('skip_binning_bamabund', false); setP('skip_bin_classification', false); }
        }, { param: '--enable_bacterial_downstream' }));
      if (!on('enable_bacterial_downstream')) return out;
    }
    out.push(chips('What should be known about every bacterium?', BAC_TOOLS.slice(0, 3).map(function (t) { return skipChip(t[0], t[1], t[2]); })));
    out.push(chips('Extra typing for specific pathogens', BAC_TOOLS.slice(3, 8).map(function (t) { return skipChip(t[0], t[1], t[2]); }),
      'Each runs only on the genus it is made for, so leaving them on costs nothing for other bacteria. Without SISTR, GenoTyphi and Paratype do not run either.'));
    out.push(chips('Which resistance-gene tools?', BAC_TOOLS.slice(8).map(function (t) { return skipChip(t[0], t[1], t[2]); }),
      'Tools disagree at the edges; several tools give a more reliable picture. All results are merged into one table (hAMRonization).'));
    var dbs = [];
    if (tool('amrfinderplus')) dbs.push(dbSource({ key: 'amrfinderplus', label: 'AMRFinderPlus', size: '~250 MB', save: 'save_amrfinderplus_db',
      paths: [{ param: 'amrfinderplus_db', label: 'AMRFinderPlus database folder', placeholder: '/db/amrfinderplus/latest', help: 'A 4.x database (a dated folder or its latest link).' }] }));
    if (tool('rgi')) dbs.push(cardDb());
    if (tool('resfinder')) dbs.push(dbSource({ key: 'resfinder', label: 'ResFinder, PointFinder, DisinFinder', size: '~100 MB', save: 'save_resfinder_db',
      paths: [{ param: 'resfinder_db', label: 'resfinder_db', placeholder: '/db/resfinder/resfinder_db' },
        { param: 'pointfinder_db', label: 'pointfinder_db', placeholder: '/db/resfinder/pointfinder_db' },
        { param: 'disinfinder_db', label: 'disinfinder_db', placeholder: '/db/resfinder/disinfinder_db' }] }));
    if (tool('deeparg')) dbs.push(dbSource({ key: 'deeparg', label: 'DeepARG', size: '~5 GB', save: 'save_deeparg_db',
      paths: [{ param: 'deeparg_db', label: 'DeepARG data folder', placeholder: '/db/deeparg', help: 'Contains database/ and model/.' }] }));
    if (tool('virulencefinder')) dbs.push(dbSource({ key: 'virulencefinder', label: 'VirulenceFinder', size: '~10 MB', save: 'save_virulencefinder_db',
      paths: [{ param: 'virulencefinder_db', label: 'VirulenceFinder database folder', placeholder: '/db/virulencefinder_db' }] }));
    if (paratype()) dbs.push(dbSource({ key: 'paratype', label: 'Paratype', size: '~50 MB', save: 'save_paratype_db',
      paths: [{ param: 'paratype_bundle', label: 'Paratype checkout', placeholder: '/db/Paratype' }] }));
    if (dbs.length) out.push(section('Databases', [h('p', { class: 'mcb-help' }, 'Downloaded only when a bacterial genome reaches the tool. MLST, PlasmidFinder, ABRicate, Kleborate, SISTR, SeqSero2 and Mykrobe use the databases inside their containers.'), h('div', { class: 'mcb-dbs' }, dbs)]));
    out.push(textInput('bacterial_decision_table', { label: 'Per-lineage decision table (optional)', placeholder: 'decisions.tsv', wide: true,
      help: 'TSV choosing steps and options per GTDB lineage. See the Cookbook for the format.' }));
    return out;
  }

  function renderRun() {
    var nf = S.nf, out = [];
    out.push(segmented('How should the tools be installed?', [
      { value: 'singularity', label: 'Singularity' }, { value: 'apptainer', label: 'Apptainer' }, { value: 'docker', label: 'Docker' },
      { value: 'conda', label: 'Conda' }, { value: 'mamba', label: 'Mamba' }
    ], nf.profile, function (v) { nf.profile = v; }, { wide: true, help: 'Containers (Singularity, Apptainer, Docker) give the same tool versions everywhere; on a shared cluster it is usually Singularity or Apptainer. Use Conda only when containers are not possible.' }));
    if (nf.profile !== 'docker') {
      out.push(nfText('cacheDir', nf.profile === 'conda' || nf.profile === 'mamba' ? 'Conda environment cache (optional)' : 'Image cache folder (optional)',
        '/shared/containers', 'Reuse downloaded images or environments across runs and projects.'));
    }
    out.push(segmented('Where should the tasks run?', [
      { value: 'local', label: 'This machine' }, { value: 'slurm', label: 'SLURM' }, { value: 'pbspro', label: 'PBS Pro' },
      { value: 'sge', label: 'SGE' }, { value: 'lsf', label: 'LSF' }
    ], nf.executor, function (v) { nf.executor = v; }, { wide: true }));
    if (nf.executor !== 'local') {
      out.push(h('div', { class: 'mcb-grid' },
        nfText('queue', 'Queue or partition', 'normal'),
        nf.executor === 'slurm' ? nfText('account', 'Account (optional)', 'my_lab') : null,
        nfText('queueSize', 'Jobs at a time (optional)', '50')));
    }
    var memNum = String(get('max_memory') || '').match(/^([\d.]+)\s*\.?\s*([KMGT]?B)$/);
    var timeNum = String(get('max_time') || '').match(/^([\d.]+)\s*\.?\s*([smhd])$/);
    out.push(section('How big may a single task be?', [h('div', { class: 'mcb-grid is-3' },
      numberInput('max_cpus', { label: 'CPUs', min: 1, step: 1 }),
      unitInput('max_memory', 'Memory', memNum ? memNum[1] : '', memNum ? memNum[2] : 'GB', ['GB', 'TB']),
      unitInput('max_time', 'Time', timeNum ? timeNum[1] : '', timeNum ? timeNum[2] : 'h', ['h', 'd'])),
    h('p', { class: 'mcb-help' }, 'Caps for any single task. Set them to what one node of your machine or cluster can give.')]));
    out.push(section('Nextflow', [
      nfText('workDir', 'Work folder (optional)', '/scratch/me/work', 'Intermediate files; put it on fast scratch space, not in backed-up storage.'),
      toggle('Resume', 'Reuse the finished tasks of an earlier run in the same folder.', nf.resume, function (v) { nf.resume = v; }, { param: '-resume' }),
      toggle('Run in the background', 'Return the prompt at once; follow progress in .nextflow.log.', nf.bg, function (v) { nf.bg = v; }, { param: '-bg' }),
      nfText('revision', 'Pipeline version', VERSION, 'A release tag, or dev for the latest code.')
    ]));
    out.push(section('Report', [
      paramToggle('skip_multiqc', 'Make one summary report for the whole run?', 'An HTML page with every step, sample and tool version. Uses MultiQC.', true),
      textInput('email', { label: 'Email when done (optional)', placeholder: 'me@example.org', mono: false, help: 'Needs a working mail setup on the machine running Nextflow.' })
    ]));
    return out;
  }
  function nfText(key, label, placeholder, help) {
    var input = h('input', {
      type: 'text', id: 'nf-' + key, class: 'mcb-input is-mono', value: S.nf[key] || '', placeholder: placeholder || '',
      spellcheck: 'false', autocomplete: 'off',
      oninput: function (e) { S.nf[key] = e.target.value.trim(); commit(false); }
    });
    return field({ label: label, help: help }, input);
  }
  function unitInput(name, label, num, unit, units) {
    var input = h('input', {
      type: 'number', id: 'in-' + name, class: 'mcb-input is-num', value: Object.prototype.hasOwnProperty.call(S.p, name) ? num : '',
      min: '1', step: 'any', placeholder: String(def(name) || '').replace(/[^\d.]/g, ''),
      oninput: function () { write(); }
    });
    var select = h('select', { class: 'mcb-select', 'aria-label': label + ' unit', onchange: function () { write(); } },
      units.map(function (u) { return h('option', { value: u, selected: u === unit ? true : null }, u); }));
    function write() { setP(name, input.value === '' ? '' : input.value + '.' + select.value); commit(false); }
    return field({ label: label, param: name }, h('div', { class: 'mcb-unit' }, input, select));
  }

  /* "More options": every other schema parameter of the step, typed from the schema */
  function renderAdvanced(step) {
    var list = step.advList || [];
    if (!list.length || (step.advTest && !step.advTest())) return null;
    var changed = list.filter(function (n) { return Object.prototype.hasOwnProperty.call(S.p, n); }).length;
    var body = h('div', { class: 'mcb-adv-body' });
    list.forEach(function (n) {
      var p = P[n], el;
      if (p.type === 'boolean') {
        el = toggle(n.replace(/_/g, ' '), p.description, on(n), (function (name) { return function (v) { setP(name, v); }; })(n), { param: '--' + n });
      } else if (p.enum) {
        var select = h('select', { class: 'mcb-select', id: 'in-' + n, onchange: (function (name) { return function (e) { setP(name, e.target.value); commit(false); }; })(n) },
          p.enum.map(function (v) { return h('option', { value: v, selected: String(get(n)) === String(v) ? true : null }, String(v) + (same(v, p.default) ? ' (default)' : '')); }));
        el = field({ label: p.description || n, param: n }, select);
      } else if (p.type === 'integer' || p.type === 'number') {
        el = numberInput(n, { label: p.description || n, step: p.type === 'integer' ? 1 : 'any' });
      } else {
        var argsLike = /_(args|options)$/.test(n);
        el = textInput(n, { label: p.description || n, wide: argsLike, help: argsLike ? 'Pasted as is and passed to the tool.' : null });
      }
      body.appendChild(el);
    });
    var details = h('details', { class: 'mcb-adv', open: S.ui['adv_' + step.id] ? true : null,
      ontoggle: function (e) { S.ui['adv_' + step.id] = e.target.open; } },
    h('summary', null, h('span', null, 'More options'),
      h('span', { class: 'mcb-adv-count' }, list.length + ' settings' + (changed ? ' · ' + changed + ' changed' : ''))),
    body);
    return details;
  }

  function renderDone() {
    var errs = problems().filter(function (p) { return p.level === 'error'; });
    var out = [];
    if (errs.length) {
      out.push(note('error', [h('strong', null, errs.length === 1 ? 'One thing to fix before running:' : errs.length + ' things to fix before running:'),
        h('ul', null, errs.map(function (e) { return h('li', null, h('button', { type: 'button', class: 'mcb-linkbtn', onclick: function () { go(e.step); } }, BYID[e.step].title), ': ' + e.msg); }))]));
    } else {
      out.push(note('ok', 'Everything needed is set. Pick how you want to run it on the right, or below on a phone.'));
    }
    var allReqs = activeSteps().map(function (s) { return renderRequirements(s.id, true); }).filter(Boolean);
    if (allReqs.length) out.push(section('What your choices need', allReqs));
    out.push(section('Three ways to run the same thing', h('dl', { class: 'mcb-ways' },
      h('dt', null, 'Command'), h('dd', null, 'Everything on one command line. Good for a single run.'),
      h('dt', null, 'params.yaml'), h('dd', null, 'The settings in a file, given with -params-file. Easy to keep with the project and edit.'),
      h('dt', null, '+ config'), h('dd', null, 'params.yaml for the pipeline settings, plus microcosm.config (given with -c) for the cluster and the image cache. Reuse the config for every project on the same machine. Pipeline parameters do not go in the -c file: some are read before it is loaded.'))));
    out.push(section('Before you start', h('ul', { class: 'mcb-checklist' },
      h('li', null, 'Nextflow 26.04 or later and Java 17 or later: ', code('nextflow -version'), '.'),
      h('li', null, 'Check every path exists on the machine that runs Nextflow.'),
      h('li', null, 'Databases you must download yourself: ', link('getting-started/databases/', 'Databases'), '.'),
      h('li', null, 'Results: open ', code((get('outdir') || './results') + '/MultiQC/'), ' first; every folder is described in ', link('reference/outputs/', 'Outputs'), '.'))));
    return out;
  }

  /* ------------------------------------------------------------------ validation */

  function filled(name) { var v = get(name); return v !== undefined && v !== null && String(v).trim() !== ''; }
  function problems() {
    var list = [];
    function err(step, msg) { list.push({ level: 'error', step: step, msg: msg }); }
    function warn(step, msg) { list.push({ level: 'warn', step: step, msg: msg }); }
    if (!S.route) { err('route', 'choose what the run should do'); return list; }
    if (genomes()) {
      if (!filled('genome_input')) err('input', 'give the genome samplesheet');
    } else {
      if (!filled('input')) err('input', get('input_format') === 'directory' ? 'give the FASTQ folder' : 'give the samplesheet');
      else if (get('input_format') !== 'directory' && !/\.csv$/i.test(get('input'))) warn('input', 'the samplesheet should be a .csv file');
      if (fromAssemblies() && !filled('assembly_input')) err('input', 'give the assemblies samplesheet');
    }
    if ((reads() || asm()) && !fromAssemblies() && !on('skip_preprocess') && S.ui.hostRef === 'local' && !filled('hostile_reference')) err('preprocess', 'give the hostile index folder');
    if (reads()) {
      if (!taxotools().length) err('profile', 'choose a profiler');
      if (useMicrocosm()) {
        if (!filled('sourmash_database')) err('profile', 'give the sourmash database');
        if (!filled('sourmash_taxonomy_csv')) err('profile', 'give the lineages file');
        if (yachtOn() && !filled('yacht_database')) err('profile', 'give the YACHT database, or turn YACHT off');
        if (yachtOn() && count(get('sourmash_database')) > 1) warn('profile', 'several sourmash databases: YACHT will be skipped');
        if (yachtOn() && filled('yacht_database') && !/\.json$/i.test(get('yacht_database'))) warn('profile', 'the YACHT database is usually its *_config.json');
      }
      if (useSylph()) {
        if (!filled('sylph_database')) err('profile', 'give the sylph database');
        if (on('enable_sylphtax') && !filled('sylph_taxonomy')) err('profile', 'give the sylph-tax metadata');
        if (on('enable_sylphtax') && filled('sylph_taxonomy') && count(get('sylph_taxonomy')) !== count(get('sylph_database'))) err('profile', 'give one sylph-tax file per sylph database');
      }
      if (on('enable_rgi_bwt') && db('card') === 'local' && !filled('rgi_preparecarddb_dir')) err('readsamr', 'give the CARD folder');
      var pp = get('readbased_postprocess');
      if (pp) {
        if (useMicrocosm() && !yachtRuns()) err('postprocess', 'post-processing needs YACHT on and one sourmash database');
        if (pp === 'rgi_bwt' && !on('enable_rgi_bwt')) err('postprocess', 'RGI tables need AMR genes in reads');
        if ((pp === 'taxburst' || pp === 'rgi_bwt') && !useMicrocosm()) err('postprocess', pp + ' uses the sourmash + YACHT result');
        if ((pp === 'all' || pp === 'phyloseq') && useSylph() && !on('enable_sylphtax')) err('postprocess', 'sylph phyloseq tables need sylph-tax');
        if (pp === 'all' && useMicrocosm() && !on('enable_rgi_bwt')) warn('postprocess', 'without AMR genes in reads, "all" writes no AMR table');
      }
    }
    if (asm()) {
      if (!fromAssemblies() && on('skip_megahit') && on('skip_spades')) err('assembly', 'keep at least one assembler');
      if (!on('skip_contig_analysis')) {
        if (!on('skip_skani_search') && !filled('skani_db')) err('contigs', 'give the skani database, or turn skani off');
        if (!on('skip_genomad_e2e') && db('genomad') === 'local' && !filled('genomad_db')) err('contigs', 'give the geNomad database');
      }
      if (binningOn()) {
        if (BINNERS.every(function (b) { return on('skip_' + b[0]); })) err('binning', 'keep at least one binner');
        if (needCheckm2() && db('checkm2') === 'local' && !filled('checkm2_db')) err(!on('skip_binning_refinement') && get('refine_tool') === 'binette' ? 'refine' : 'binqc', 'give the CheckM2 database');
        if (!on('skip_binqc') && !binqcTools().length) err('binqc', 'choose a quality tool');
        if (qcTool('checkm') && db('checkm') === 'local' && !filled('checkm_db')) err('binqc', 'give the CheckM database');
        if (qcTool('busco') && db('busco') === 'local' && !filled('busco_db')) err('binqc', 'give the BUSCO folder');
        if (!on('skip_binqc') && on('run_gunc') && db('gunc') === 'local' && !filled('gunc_db')) err('binqc', 'give the GUNC database');
        if (bakta() && db('bakta') === 'local' && !filled('annotation_bakta_db')) err('annotation', 'give the Bakta database');
        if (bakta() && db('bakta') !== 'local' && get('annotation_bakta_db_downloadtype') === 'full') warn('annotation', 'the full Bakta database is about 85 GB');
      } else if (on('enable_bacterial_downstream')) err('binning', 'bacterial typing needs binning');
    }
    if (classifyOn()) {
      if (get('bin_classifier') === 'gtdbtk') {
        if (!filled('gtdbtk_db')) err('classify', 'give the GTDB-Tk reference data');
        if (memGB() !== null && memGB() < 64) warn('run', 'GTDB-Tk needs 64-128 GB; raise the memory of the largest task');
      }
      else {
        if (!filled('sourmash_database')) err('classify', 'give the sourmash database');
        if (!filled('sourmash_taxonomy_csv')) err('classify', 'give the lineages file');
      }
    }
    if (downstreamOn()) {
      [['amrfinderplus', 'amrfinderplus', 'amrfinderplus_db'], ['rgi', 'card', 'rgi_preparecarddb_dir'], ['resfinder', 'resfinder', 'resfinder_db'],
        ['resfinder', 'resfinder', 'pointfinder_db'], ['resfinder', 'resfinder', 'disinfinder_db'], ['deeparg', 'deeparg', 'deeparg_db'],
        ['virulencefinder', 'virulencefinder', 'virulencefinder_db']].forEach(function (d) {
        if (tool(d[0]) && db(d[1]) === 'local' && !filled(d[2])) err('bacterial', 'give ' + d[2].replace(/_/g, ' '));
      });
      if (paratype() && db('paratype') === 'local' && !filled('paratype_bundle')) err('bacterial', 'give the Paratype checkout');
      if (tool('argnorm') && !['amrfinderplus', 'abricate', 'resfinder', 'deeparg'].some(tool)) warn('bacterial', 'ARGnorm has no AMR table to normalise');
      if ((tool('genotyphi') || tool('paratype')) && on('skip_sistr')) warn('bacterial', 'GenoTyphi and Paratype only run after SISTR');
      if (filled('bacterial_decision_table') && !/\.tsv$/i.test(get('bacterial_decision_table'))) err('bacterial', 'the decision table must be a .tsv file');
    }
    if (S.nf.profile === 'conda' || S.nf.profile === 'mamba') warn('run', 'Conda builds every environment on first use; containers are faster and more reproducible');
    return list;
  }

  /* ------------------------------------------------------------------ output */

  function activeSteps() { return STEPS.filter(function (s) { return s.when(); }); }
  function collect() {
    var out = [], seen = {};
    activeSteps().forEach(function (s) {
      var names = s.params.filter(function (e) { return !e.test || e.test(); }).map(function (e) { return e.name; });
      if (!s.advTest || s.advTest()) names = names.concat(s.advList || []);
      names.forEach(function (n) {
        if (seen[n] || !Object.prototype.hasOwnProperty.call(S.p, n)) return;
        seen[n] = true;
        out.push([n, S.p[n]]);
      });
    });
    return out;
  }
  function typed(name, v) {
    var t = P[name] && P[name].type;
    if (t === 'boolean') return v === true || v === 'true';
    if ((t === 'integer' || t === 'number') && v !== '' && !isNaN(Number(v))) return Number(v);
    return String(v);
  }

  function shq(s) { return /^[A-Za-z0-9_\/.,:@%+=~-]+$/.test(s) ? s : "'" + s.replace(/'/g, "'\\''") + "'"; }
  // Nextflow reads "--a --plus" as a flag, so a value starting with "-" is joined with "=".
  function cliParam(name, v) {
    v = typed(name, v);
    if (v === true) return '--' + name;
    if (v === false) return '--' + name + ' false';
    var s = String(v);
    return s.charAt(0) === '-' ? '--' + name + '=' + shq(s) : '--' + name + ' ' + shq(s);
  }
  function yamlValue(name, v) {
    v = typed(name, v);
    if (typeof v === 'boolean' || typeof v === 'number') return String(v);
    return "'" + String(v).replace(/'/g, "''") + "'";
  }

  function cacheEnv() {
    var nf = S.nf;
    if (!nf.cacheDir || nf.profile === 'docker') return null;
    var v = { singularity: 'NXF_SINGULARITY_CACHEDIR', apptainer: 'NXF_APPTAINER_CACHEDIR', conda: 'NXF_CONDA_CACHEDIR', mamba: 'NXF_CONDA_CACHEDIR' }[nf.profile];
    return 'export ' + v + '=' + shq(nf.cacheDir);
  }
  function baseLines(withExecutor) {
    var nf = S.nf, lines = ['nextflow run ' + PIPELINE + ' -r ' + shq(nf.revision || VERSION), '-profile ' + nf.profile];
    if (nf.workDir) lines.push('-work-dir ' + shq(nf.workDir));
    if (withExecutor && nf.executor !== 'local') {
      lines.push('-process.executor ' + nf.executor);
      if (nf.queue) lines.push('-process.queue ' + shq(nf.queue));
      if (nf.executor === 'slurm' && nf.account) lines.push("-process.clusterOptions=" + shq('--account=' + nf.account));
      if (nf.queueSize) lines.push('-executor.queueSize ' + shq(nf.queueSize));
    }
    if (nf.resume) lines.push('-resume');
    if (nf.bg) lines.push('-bg');
    return lines;
  }
  function joinCmd(lines) { return lines.join(' \\\n    '); }

  function outputs() {
    var params = collect();
    var env = cacheEnv();
    var cmd = joinCmd(baseLines(true).concat(params.map(function (p) { return cliParam(p[0], p[1]); })));
    var yamlFile = ['# microcosm ' + VERSION + ' parameters (made with the command builder)']
      .concat(params.map(function (p) { return p[0] + ': ' + yamlValue(p[0], p[1]); })).join('\n') + '\n';
    var yamlCmd = joinCmd(baseLines(true).concat(['-params-file params.yaml']));
    // Parameters stay in params.yaml: the pipeline config reads some of them (outdir, save_*, max_*) while it is
    // parsed, before a -c file is merged, so a params block in microcosm.config would be partly ignored.
    var cfg = ['// microcosm ' + VERSION + ' compute settings (made with the command builder)',
      '// nextflow run ' + PIPELINE + ' -r ' + (S.nf.revision || VERSION) + ' -profile ' + S.nf.profile + ' -params-file params.yaml -c microcosm.config'];
    var nf = S.nf;
    var infra = nf.executor !== 'local' || (nf.cacheDir && nf.profile !== 'docker');
    if (!infra) cfg.push('', '// Nothing to set yet: choose a cluster or an image cache folder in "Where it runs".', '// Keep this file for your own process settings (e.g. withName selectors).');
    if (nf.executor !== 'local') {
      cfg.push('', 'process {', "    executor = '" + nf.executor + "'");
      if (nf.queue) cfg.push('    queue = ' + gq(nf.queue));
      if (nf.executor === 'slurm' && nf.account) cfg.push('    clusterOptions = ' + gq('--account=' + nf.account));
      cfg.push('}');
      if (nf.queueSize && !isNaN(Number(nf.queueSize))) cfg.push('', 'executor {', '    queueSize = ' + Number(nf.queueSize), '}');
    }
    if (nf.cacheDir && nf.profile !== 'docker') {
      var scope = nf.profile === 'mamba' ? 'conda' : nf.profile;
      cfg.push('', scope + ' {', '    cacheDir = ' + gq(nf.cacheDir), '}');
    }
    var cfgFile = cfg.join('\n') + '\n';
    var cfgCmd = joinCmd(baseLines(false).concat(['-params-file params.yaml', '-c microcosm.config']));
    return {
      count: params.length,
      cmd: { run: (env ? env + '\n' : '') + cmd },
      yaml: { file: yamlFile, run: (env ? env + '\n' : '') + yamlCmd, name: 'params.yaml' },
      cfg: { file: cfgFile, run: cfgCmd, name: 'microcosm.config', yaml: yamlFile }
    };
  }
  function gq(s) { return "'" + String(s).replace(/\\/g, '\\\\').replace(/'/g, "\\'") + "'"; }

  /* tiny highlighters: flags, values, keys, comments */
  function hlShell(text) {
    return text.split('\n').map(function (line) {
      return esc(line)
        .replace(/^(export)(\s+)([A-Z_]+)/, '<span class="t-kw">$1</span>$2<span class="t-key">$3</span>')
        .replace(/^(nextflow run)(\s+)(\S+)/, '<span class="t-kw">$1</span>$2<span class="t-str">$3</span>')
        .replace(/(^|\s)(--?[A-Za-z][\w.]*)(?=[\s=]|$)/g, '$1<span class="t-flag">$2</span>')
        .replace(/\\$/, '<span class="t-dim">\\</span>');
    });
  }
  function hlYaml(text) {
    return text.replace(/\n$/, '').split('\n').map(function (line) {
      if (/^#/.test(line)) return '<span class="t-dim">' + esc(line) + '</span>';
      var m = line.match(/^([\w]+)(:\s)(.*)$/);
      return m ? '<span class="t-key">' + m[1] + '</span>' + m[2] + valueSpan(m[3]) : esc(line);
    });
  }
  function hlGroovy(text) {
    return text.replace(/\n$/, '').split('\n').map(function (line) {
      if (/^\s*\/\//.test(line)) return '<span class="t-dim">' + esc(line) + '</span>';
      var m = line.match(/^(\s*)([\w.]+)(\s*=\s*)(.*)$/);
      if (m) return m[1] + '<span class="t-key">' + m[2] + '</span>' + m[3] + valueSpan(m[4]);
      return esc(line).replace(/^(\w+)(\s*\{)/, '<span class="t-kw">$1</span>$2');
    });
  }
  function valueSpan(v) { return /^'/.test(v) ? '<span class="t-str">' + esc(v) + '</span>' : '<span class="t-num">' + esc(v) + '</span>'; }

  /* ------------------------------------------------------------------ layout and rendering */

  var els = {};
  var lastLines = {};

  function build() {
    mount.innerHTML = '';
    mount.removeAttribute('data-state');
    els.rail = h('nav', { class: 'mcb-rail', 'aria-label': 'Builder steps' });
    els.stage = h('section', { class: 'mcb-stage', 'aria-live': 'polite' });
    els.out = h('aside', { class: 'mcb-out', 'aria-label': 'Generated command' });
    els.progress = h('div', { class: 'mcb-progress', 'aria-hidden': 'true' }, h('span', { class: 'mcb-progress-bar' }));
    mount.appendChild(els.progress);
    mount.appendChild(els.rail);
    mount.appendChild(els.stage);
    mount.appendChild(els.out);
    buildOutput();
    renderStep(true);
    refresh();
  }

  function stepIndex(id) { var a = activeSteps(); for (var i = 0; i < a.length; i++) if (a[i].id === id) return i; return -1; }
  function nextStep(id) { var a = activeSteps(), i = stepIndex(id); return a[Math.min(i + 1, a.length - 1)].id; }
  function prevStep(id) { var a = activeSteps(), i = stepIndex(id); return a[Math.max(i - 1, 0)].id; }
  function go(id) {
    clearTimeout(autoAdvance);
    if (!BYID[id] || !BYID[id].when()) return;
    S.seen[S.step] = true;
    S.step = id;
    S.seen[id] = true;
    renderStep(true);
    refresh();
    save();
    var top = mount.getBoundingClientRect().top + window.scrollY - 80;
    if (window.scrollY > top) window.scrollTo({ top: top, behavior: REDUCED ? 'auto' : 'smooth' });
  }

  function renderRail() {
    var errsBy = {};
    problems().forEach(function (p) { if (p.level === 'error') errsBy[p.step] = (errsBy[p.step] || 0) + 1; });
    var act = activeSteps();
    var cur = stepIndex(S.step);
    var list = h('ol', { class: 'mcb-stations' + (S.route ? ' route-' + S.route : '') });
    act.forEach(function (s, i) {
      var state = i === cur ? 'is-current' : (S.seen[s.id] ? 'is-done' : 'is-todo');
      var sum = S.seen[s.id] || i === cur ? s.summary() : '';
      list.appendChild(h('li', { class: 'mcb-station ' + state + (errsBy[s.id] && S.seen[s.id] && i !== cur ? ' has-error' : ''), style: '--i:' + i },
        h('button', { type: 'button', class: 'mcb-station-btn', 'aria-current': i === cur ? 'step' : null,
          disabled: !S.route && s.id !== 'route' ? true : null, onclick: function () { go(s.id); } },
        h('span', { class: 'mcb-dot', 'aria-hidden': 'true' }),
        h('span', { class: 'mcb-station-text' }, h('span', { class: 'mcb-station-title' }, s.title),
          sum ? h('span', { class: 'mcb-station-sum' }, sum) : null))));
    });
    els.rail.innerHTML = '';
    els.rail.appendChild(list);
    var pct = act.length > 1 ? Math.round(cur / (act.length - 1) * 100) : 0;
    els.progress.firstChild.style.width = pct + '%';
  }

  function renderStep(animate) {
    var focusId = document.activeElement && mount.contains(document.activeElement) ? document.activeElement.id : null;
    var step = BYID[S.step];
    if (!step || !step.when()) { S.step = 'route'; step = BYID.route; }
    var act = activeSteps(), i = stepIndex(step.id);
    var card = h('div', { class: 'mcb-card' + (animate && !REDUCED ? ' is-entering' : '') },
      h('header', { class: 'mcb-card-head' },
        h('p', { class: 'mcb-kicker' }, S.route ? 'Step ' + (i + 1) + ' of ' + act.length : 'Start here'),
        h('h2', { class: 'mcb-card-title', tabindex: '-1', id: 'mcb-title' }, step.lead),
        null),
      h('div', { class: 'mcb-card-body' }, step.render(), h('div', { class: 'mcb-reqs-slot' }), renderAdvanced(step)),
      h('footer', { class: 'mcb-card-foot' },
        i > 0 ? h('button', { type: 'button', class: 'mcb-btn is-ghost', onclick: function () { go(prevStep(step.id)); } }, '← Back') : h('span'),
        step.id === 'done' ? h('button', { type: 'button', class: 'mcb-btn is-ghost', onclick: resetAll }, 'Start over')
          : h('button', { type: 'button', class: 'mcb-btn is-primary', disabled: !S.route ? true : null, onclick: function () { go(nextStep(step.id)); } },
            S.route ? 'Next: ' + BYID[nextStep(step.id)].title + ' →' : 'Choose a goal')));
    els.stage.innerHTML = '';
    els.stage.appendChild(card);
    markRequired();
    if (focusId && document.getElementById(focusId)) document.getElementById(focusId).focus();
    else if (animate && S.route && mount.dataset.ready) document.getElementById('mcb-title').focus({ preventScroll: true });
    mount.dataset.ready = '1';
  }

  function markRequired() {
    var seen = S.seen[S.step];
    Array.prototype.forEach.call(mount.querySelectorAll('[data-req]'), function (input) {
      var empty = !input.value.trim();
      input.setAttribute('aria-invalid', empty && seen ? 'true' : 'false');
      input.closest('.mcb-field').classList.toggle('is-missing', empty);
    });
  }

  function buildOutput() {
    var tabs = h('div', { class: 'mcb-tabs', role: 'tablist', 'aria-label': 'Output format' });
    [['cmd', 'Command'], ['yaml', 'params.yaml'], ['cfg', '+ config']].forEach(function (t) {
      tabs.appendChild(h('button', { type: 'button', role: 'tab', id: 'mcb-tab-' + t[0], 'aria-controls': 'mcb-panel',
        class: 'mcb-tab', onclick: function () { S.fmt = t[0]; refresh(); save(); },
        onkeydown: function (e) {
          var order = ['cmd', 'yaml', 'cfg'], i = order.indexOf(S.fmt);
          if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') {
            S.fmt = order[(i + (e.key === 'ArrowRight' ? 1 : 2)) % 3]; refresh(); save();
            document.getElementById('mcb-tab-' + S.fmt).focus();
          }
        } }, t[1]));
    });
    els.tabs = tabs;
    els.status = h('div', { class: 'mcb-status' });
    els.panel = h('div', { class: 'mcb-panel', role: 'tabpanel', id: 'mcb-panel' });
    els.share = h('button', { type: 'button', class: 'mcb-linkbtn', onclick: shareLink }, 'Copy a link to these settings');
    els.out.appendChild(h('div', { class: 'mcb-out-inner' },
      h('div', { class: 'mcb-out-head' }, h('h2', { class: 'mcb-out-title' }, 'Your run'), els.status),
      tabs, els.panel,
      h('div', { class: 'mcb-out-foot' }, els.share, h('span', { class: 'mcb-sep', 'aria-hidden': 'true' }, '·'),
        h('button', { type: 'button', class: 'mcb-linkbtn', onclick: resetAll }, 'Start over'))));
  }

  function codeBlock(key, title, text, lang, actions) {
    var lines = lang === 'yaml' ? hlYaml(text) : lang === 'groovy' ? hlGroovy(text) : hlShell(text);
    var raw = text.replace(/\n$/, '').split('\n');
    var before = lastLines[key] || null;
    var pre = h('pre', { class: 'mcb-code', tabindex: '0', 'aria-label': title });
    pre.innerHTML = lines.map(function (l, i) {
      var fresh = before && before.indexOf(raw[i]) < 0 && !REDUCED;
      return '<span class="mcb-line' + (fresh ? ' is-new' : '') + '">' + l + '</span>';
    }).join('');
    lastLines[key] = raw;
    return h('div', { class: 'mcb-block' },
      h('div', { class: 'mcb-block-head' }, h('span', { class: 'mcb-block-title' }, title), h('span', { class: 'mcb-actions' }, actions)),
      pre);
  }
  function copyBtn(text, label) {
    return h('button', { type: 'button', class: 'mcb-btn is-small is-primary', onclick: function (e) { copy(text, e.currentTarget); } }, label || 'Copy');
  }
  function dlBtn(name, text) {
    return h('button', { type: 'button', class: 'mcb-btn is-small', onclick: function () { download(name, text, 'text/plain'); } }, 'Download');
  }

  function refresh() {
    renderRail();
    var slot = els.stage.querySelector('.mcb-reqs-slot');
    if (slot) { slot.innerHTML = ''; append(slot, S.step === 'done' ? null : renderRequirements(S.step)); }
    var o = outputs();
    var probs = problems();
    var errs = probs.filter(function (p) { return p.level === 'error'; });
    var warns = probs.filter(function (p) { return p.level === 'warn'; });
    els.status.className = 'mcb-status ' + (!S.route ? 'is-idle' : errs.length ? 'is-error' : 'is-ok');
    els.status.textContent = !S.route ? 'Choose a goal to start' : errs.length ? errs.length + (errs.length === 1 ? ' thing' : ' things') + ' to fill in' : 'Ready to run · ' + o.count + ' settings changed';
    Array.prototype.forEach.call(els.tabs.children, function (b) {
      var active = b.id === 'mcb-tab-' + S.fmt;
      b.setAttribute('aria-selected', active ? 'true' : 'false');
      b.setAttribute('tabindex', active ? '0' : '-1');
      b.classList.toggle('is-on', active);
    });
    els.panel.setAttribute('aria-labelledby', 'mcb-tab-' + S.fmt);
    var blocks = [];
    if (S.fmt === 'cmd') {
      blocks.push(codeBlock('cmd', 'Run', o.cmd.run, 'sh', [copyBtn(o.cmd.run)]));
    } else if (S.fmt === 'yaml') {
      blocks.push(codeBlock('yaml-file', '1 · Save as params.yaml', o.yaml.file, 'yaml', [dlBtn('params.yaml', o.yaml.file), copyBtn(o.yaml.file)]));
      blocks.push(codeBlock('yaml-run', '2 · Run, next to the file', o.yaml.run, 'sh', [copyBtn(o.yaml.run)]));
    } else {
      blocks.push(codeBlock('cfg-yaml', '1 · Save as params.yaml', o.cfg.yaml, 'yaml', [dlBtn('params.yaml', o.cfg.yaml), copyBtn(o.cfg.yaml)]));
      blocks.push(codeBlock('cfg-file', '2 · Save as microcosm.config', o.cfg.file, 'groovy', [dlBtn('microcosm.config', o.cfg.file), copyBtn(o.cfg.file)]));
      blocks.push(codeBlock('cfg-run', '3 · Run, next to both files', o.cfg.run, 'sh', [copyBtn(o.cfg.run)]));
    }
    if (S.route && (errs.length || warns.length)) {
      blocks.push(h('ul', { class: 'mcb-issues' }, errs.concat(warns).slice(0, 8).map(function (p) {
        return h('li', { class: 'is-' + p.level },
          h('button', { type: 'button', class: 'mcb-linkbtn', onclick: function () { go(p.step); } }, BYID[p.step].title),
          ' ' + p.msg);
      })));
    }
    els.panel.innerHTML = '';
    append(els.panel, blocks);
  }

  /* ------------------------------------------------------------------ copy, download, share, storage */

  function copy(text, btn) {
    function done(ok) {
      var old = btn.textContent;
      btn.textContent = ok ? 'Copied' : 'Press Ctrl+C';
      btn.classList.add('is-flash');
      setTimeout(function () { btn.textContent = old; btn.classList.remove('is-flash'); }, 1400);
    }
    if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(text).then(function () { done(true); }, function () { done(fallback()); });
    else done(fallback());
    function fallback() {
      var ta = h('textarea', { class: 'mcb-offscreen' });
      ta.value = text; document.body.appendChild(ta); ta.select();
      var ok = false; try { ok = document.execCommand('copy'); } catch (e) { ok = false; }
      document.body.removeChild(ta); return ok;
    }
  }
  function download(name, text, type) {
    var url = URL.createObjectURL(new Blob([text], { type: type + ';charset=utf-8' }));
    var a = h('a', { href: url, download: name });
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }
  function snapshot() { return { v: 1, route: S.route, step: S.step, seen: S.seen, p: S.p, ui: S.ui, nf: S.nf, fmt: S.fmt }; }
  function encode(obj) { return btoa(unescape(encodeURIComponent(JSON.stringify(obj)))).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, ''); }
  function decode(s) { s = s.replace(/-/g, '+').replace(/_/g, '/'); return JSON.parse(decodeURIComponent(escape(atob(s)))); }
  function restore(obj) {
    if (!obj || obj.v !== 1) return false;
    var base = freshState();
    S = { route: obj.route || null, step: obj.step || 'route', seen: obj.seen || {}, p: obj.p || {},
      ui: Object.assign(base.ui, obj.ui || {}), nf: Object.assign(base.nf, obj.nf || {}), fmt: obj.fmt || 'cmd' };
    return true;
  }
  var saveTimer = null;
  function save() {
    clearTimeout(saveTimer);
    saveTimer = setTimeout(function () {
      try { localStorage.setItem(STORE_KEY, JSON.stringify(snapshot())); } catch (e) { /* storage off */ }
      if (location.hash.indexOf('#cfg=') === 0) history.replaceState(null, '', location.pathname + location.search);
    }, 250);
  }
  function shareLink(e) {
    var url = location.origin + location.pathname + '#cfg=' + encode(snapshot());
    copy(url, e.currentTarget);
  }
  function resetAll() {
    if (S.route && !window.confirm('Clear every choice and start over?')) return;
    S = freshState();
    lastLines = {};
    try { localStorage.removeItem(STORE_KEY); } catch (e) { /* storage off */ }
    renderStep(true); refresh();
  }

  /* ------------------------------------------------------------------ start */

  fetch(SCHEMA_URL).then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); }).then(function (schema) {
    VERSION = schema.version || VERSION;
    schema.groups.forEach(function (g) {
      GROUPS[g.id] = Object.keys(g.params);
      Object.keys(g.params).forEach(function (n) { P[n] = Object.assign({ group: g.id }, g.params[n]); });
    });
    S = freshState();
    buildAdvanced();
    var loaded = false;
    if (location.hash.indexOf('#cfg=') === 0) { try { loaded = restore(decode(location.hash.slice(5))); } catch (e) { loaded = false; } }
    if (!loaded) { try { restore(JSON.parse(localStorage.getItem(STORE_KEY) || 'null')); } catch (e) { /* storage off */ } }
    build();
  }).catch(function () {
    mount.innerHTML = '<p class="mcb-loading">The parameter list could not be loaded. Use the <a href="../cookbook/">Cookbook</a> and the <a href="../../reference/parameters/">Parameters</a> list instead.</p>';
  });
})();
