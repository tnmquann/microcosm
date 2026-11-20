# Sphinx configuration for the microcosm documentation (Read the Docs).
import json
import os

# -- Project information
project = "microcosm"
copyright = "2026, Minh-Quan Ton-Ngoc, Si-Nguyen Mai"
version = "2.0.0"
release = "2.0.0"
html_context = {"doc_author": "Minh-Quan Ton-Ngoc, Si-Nguyen Mai"}

# -- General configuration
extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_togglebutton",
    "sphinx_tabs.tabs",
    "sphinxcontrib.mermaid",
]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
myst_enable_extensions = ["colon_fence", "deflist", "attrs_inline", "attrs_block", "html_image", "dollarmath", "amsmath"]
myst_heading_anchors = 3
mermaid_version = "11.4.0"
copybutton_prompt_text = r"\$ "
copybutton_prompt_is_regexp = True
suppress_warnings = ["myst.header"]
smartquotes = False  # keep "--param" as typed
mermaid_d3_zoom = False
mermaid_init_js = (
    "mermaid.initialize({startOnLoad: true, theme: 'base', themeVariables: {"
    "primaryColor: '#e7f6ee', primaryBorderColor: '#23af62', primaryTextColor: '#173a27', "
    "lineColor: '#17703f', fontFamily: 'inherit'}});"
)

templates_path = ["_templates"]
# docs/images holds the logo and diagrams shared with the README
html_static_path = ["_static", "../images"]
html_css_files = ["css/custom.css"]

# -- HTML output
# Site root (GitHub Pages serves the docs under /microcosm/); used for the logo and the top navigation links
BASE = os.environ.get("DOCS_BASEURL", "https://tnmquann.github.io/microcosm/")
html_baseurl = BASE
html_title = "microcosm"
html_theme = "shibuya"
html_copy_source = False
html_show_sourcelink = False
html_theme_options = {
    "accent_color": "grass",
    "color_mode": "auto",
    "logo_target": BASE,
    "light_logo": "_static/microcosm_logo.png",
    "dark_logo": "_static/microcosm_logo_dark.png",
    "github_url": "https://github.com/tnmquann/microcosm",
    "globaltoc_expand_depth": 1,
    "toctree_collapse": False,
    "nav_links": [
        {"title": "Quick start", "url": BASE + "getting-started/quickstart/"},
        {"title": "Pipeline", "url": BASE + "pipeline/overview/"},
        {"title": "Command builder", "url": BASE + "usage/builder/"},
        {"title": "Cookbook", "url": BASE + "usage/cookbook/"},
        {"title": "Parameters", "url": BASE + "reference/parameters/"},
        {"title": "Releases", "url": "https://github.com/tnmquann/microcosm/releases"},
    ],
}

# -- Parameters page, generated from the pipeline schema at every build
SCHEMA = os.path.join(os.path.dirname(__file__), "..", "..", "nextflow_schema.json")
HIDDEN = {"pipeline_root"}


def _cell(text):
    return " ".join(str(text).split()).replace("|", "\\|")


def _default(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return f"`{str(value).lower()}`"
    if isinstance(value, list):
        return "`" + ", ".join(map(str, value)) + "`"
    return f"`{value}`"


def write_parameters(app):
    with open(SCHEMA) as handle:
        schema = json.load(handle)
    lines = [
        "# Parameters",
        "",
        "Every parameter of microcosm, generated from `nextflow_schema.json` (the same file Nextflow validates your",
        "run against). Give them on the command line (`--name value`) or in a `-params-file`.",
        "Run `nextflow run tnmquann/microcosm -r 2.0.0 --help` for the same list in the terminal.",
        "",
        ":::{tip}",
        "Use the search box (top right) or your browser search to find a parameter. Booleans can be given as",
        "`--flag` (true) or `--flag false`.",
        ":::",
        "",
    ]
    for group in schema.get("$defs", schema.get("definitions", {})).values():
        props = {k: v for k, v in group.get("properties", {}).items() if k not in HIDDEN}
        if not props:
            continue
        lines += [f"## {group.get('title', 'Options')}", ""]
        if group.get("description"):
            lines += [_cell(group["description"]), ""]
        lines += ["| Parameter | Type | Default | Description |", "| :--- | :--- | :--- | :--- |"]
        details = []
        for name, prop in props.items():
            kind = prop.get("type", "")
            desc = _cell(prop.get("description", ""))
            if prop.get("enum"):
                desc = desc.rstrip(".") + ". Choices: " + ", ".join(f"`{v}`" for v in prop["enum"]) + "."
            if prop.get("help_text"):
                desc += f" {{ref}}`more <help-{name.replace('_', '-')}>`"
                details.append((name, prop["help_text"]))
            lines.append(f"| `--{name}` | {kind} | {_default(prop.get('default'))} | {desc} |")
        lines.append("")
        for name, help_text in details:
            lines += [f"(help-{name.replace('_', '-')})=", f":::{{dropdown}} --{name}", help_text.strip(), ":::", ""]
    out = os.path.join(app.srcdir, "reference", "parameters.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as handle:
        handle.write("\n".join(lines) + "\n")


BUILDER_PAGE = "usage/builder"


def write_builder_schema(app):
    """Parameter types, defaults and help for the command builder, from the same schema."""
    with open(SCHEMA) as handle:
        schema = json.load(handle)
    keep = ("type", "default", "enum", "description", "format", "pattern")
    groups = []
    for key, group in schema.get("$defs", {}).items():
        props = {name: {k: v for k, v in prop.items() if k in keep} for name, prop in group.get("properties", {}).items()}
        groups.append({"id": key, "title": group.get("title", key), "params": props})
    out = os.path.join(app.srcdir, "_static", "builder", "schema.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as handle:
        json.dump({"version": release, "groups": groups}, handle, separators=(",", ":"))


def builder_assets(app, pagename, templatename, context, doctree):
    if pagename == BUILDER_PAGE:
        app.add_css_file("builder/builder.css")
        app.add_js_file("builder/builder.js", defer="defer")


def setup(app):
    app.connect("builder-inited", write_parameters)
    app.connect("builder-inited", write_builder_schema)
    app.connect("html-page-context", builder_assets)
