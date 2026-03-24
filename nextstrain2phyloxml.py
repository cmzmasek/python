#!/usr/bin/env python3
"""Convert Nextclade JSON (Auspice v2 format) to PhyloXML."""

import json
import sys
import argparse
import os
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, ElementTree, indent


VERSION = "1.0.0"

PHYLOXML_NS = "http://www.phyloxml.org"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
SCHEMA_LOC = "http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd"

EPILOG = """
INPUT FORMAT
  Accepts Nextclade's Auspice JSON v2 tree file, typically named
  nextclade.auspice.json and produced by:

      nextclade run --output-tree nextclade.auspice.json ...

  Two JSON layouts are supported:

  1. Full Auspice v2 envelope (recommended):
       { "version": "v2", "meta": { "title": "..." }, "tree": { ... } }

  2. Bare tree object (the "tree" value only):
       { "name": "root", "node_attrs": {...}, "children": [...] }

DIRECTORY MODE
  When INPUT is a directory, all *.json files inside it are converted.
  Each output file is written to the output directory with the same base
  name and a .xml extension (e.g. sample.auspice.json →
  sample.auspice.xml).  The output directory is created if it does
  not exist.  Use --recursive / -r to also process subdirectories
  (subdirectory structure is mirrored under the output directory).
  Use --suffix to change the output file extension (e.g. --suffix .phyloxml).
  Files that fail to parse or convert are reported and skipped; the
  process exits with code 1 if any file failed.

FIELD MAPPING
  Nextclade field                 PhyloXML element
  ─────────────────────────────── ──────────────────────────────────────────
  node.name                       <name>
  node_attrs.div (delta)          <branch_length>   (child div − parent div)
  node_attrs.clade_membership     <taxonomy><id provider="nextclade_clade">
  node_attrs.num_date             <property ref="nextclade:num_date" datatype="xsd:decimal">
  node_attrs.country              <property ref="nextclade:country">
  node_attrs.division             <property ref="nextclade:division">
  node_attrs.location             <property ref="nextclade:location">
  node_attrs.accession            <sequence><accession source="genbank">
  branch_attrs.mutations.nuc      <property ref="nextclade:nuc_mutations">
  branch_attrs.mutations.<gene>   <property ref="nextclade:aa_mutations:<gene>">
  branch_attrs.labels.*           <property ref="nextclade:label:*">
  other node_attrs                <property ref="nextclade:<key>">
  meta.title / meta.name          <phylogeny><name>
  meta.description                <phylogeny><description>

EXAMPLES
  Convert a single file and write to stdout:
      python3 nextclade2phyloxml.py nextclade.auspice.json

  Convert a single file and save the result:
      python3 nextclade2phyloxml.py nextclade.auspice.json -o tree.phyloxml

  Pipe from another command:
      cat nextclade.auspice.json | python3 nextclade2phyloxml.py > tree.phyloxml

  Convert all JSON files in a directory:
      python3 nextclade2phyloxml.py input_dir/ -o output_dir/

  Convert all JSON files recursively:
      python3 nextclade2phyloxml.py input_dir/ -o output_dir/ --recursive

  Use with Nextclade CLI in one step:
      nextclade run --input-dataset sars-cov-2 \\
          --output-tree /dev/stdout sequences.fasta \\
        | python3 nextclade2phyloxml.py > tree.phyloxml

NOTES
  - Branch lengths are computed as the difference in cumulative divergence
    (node_attrs.div) between a child and its parent. If div is absent from a
    node, that node will have no <branch_length> element.
  - Node attributes stored as { "value": X, "confidence": ... } are
    unwrapped automatically; only the scalar value is written to PhyloXML.
  - Requires Python 3.9+ (uses xml.etree.ElementTree.indent).
"""


class ConversionError(Exception):
    """Raised for recoverable, user-facing errors during conversion."""


def _die(msg, hint=None):
    """Print an error message to stderr and exit with code 1."""
    print(f"ERROR: {msg}", file=sys.stderr)
    if hint:
        print(f"HINT:  {hint}", file=sys.stderr)
    sys.exit(1)


def _scalar(value):
    """Extract scalar from an Auspice node_attr value.

    Auspice stores many attributes as either a plain scalar or an object of
    the form {"value": X, "confidence": {...}}.  This helper normalises both.
    """
    if isinstance(value, dict):
        return value.get("value")
    return value


def _build_clade(node, parent_div=0.0, path="root"):
    """Recursively build a <clade> Element from an Auspice tree node.

    Args:
        node:       dict representing one tree node.
        parent_div: cumulative divergence of the parent node (float).
        path:       dot-separated node path used in error messages.

    Returns:
        xml.etree.ElementTree.Element
    """
    if not isinstance(node, dict):
        raise ConversionError(
            f"Tree node at path '{path}' is not a JSON object "
            f"(got {type(node).__name__}).\n"
            "Each node must be a dict with at least a 'name' key."
        )

    clade = Element("clade")

    name = node.get("name", "")
    node_label = name or f"<unnamed @ {path}>"

    if name:
        SubElement(clade, "name").text = str(name)

    node_attrs = node.get("node_attrs", {})
    branch_attrs = node.get("branch_attrs", {})

    if not isinstance(node_attrs, dict):
        raise ConversionError(
            f"'node_attrs' of node '{node_label}' must be a JSON object, "
            f"got {type(node_attrs).__name__}."
        )
    if not isinstance(branch_attrs, dict):
        raise ConversionError(
            f"'branch_attrs' of node '{node_label}' must be a JSON object, "
            f"got {type(branch_attrs).__name__}."
        )

    # Branch length: difference in cumulative divergence from parent
    div_raw = node_attrs.get("div")
    div = _scalar(div_raw)
    if div is not None:
        try:
            div = float(div)
        except (TypeError, ValueError):
            raise ConversionError(
                f"'div' of node '{node_label}' cannot be converted to a number "
                f"(value: {div_raw!r}).\n"
                "Expected a numeric divergence value such as 0.0012."
            )
        branch_length = div - float(parent_div)
        SubElement(clade, "branch_length").text = str(branch_length)
        current_div = div
    else:
        current_div = parent_div

    # Clade membership → taxonomy/id
    clade_val = _scalar(node_attrs.get("clade_membership"))
    if clade_val:
        taxonomy = SubElement(clade, "taxonomy")
        SubElement(taxonomy, "id", provider="nextclade_clade").text = str(clade_val)

    # Collection date
    num_date = _scalar(node_attrs.get("num_date"))
    if num_date is not None:
        prop = SubElement(
            clade, "property",
            ref="nextclade:num_date",
            applies_to="clade",
            datatype="xsd:decimal",
        )
        prop.text = str(num_date)

    # Geographic location
    for attr in ("country", "division", "location"):
        val = _scalar(node_attrs.get(attr))
        if val:
            prop = SubElement(
                clade, "property",
                ref=f"nextclade:{attr}",
                applies_to="clade",
                datatype="xsd:string",
            )
            prop.text = str(val)

    # Accession → sequence element
    accession = _scalar(node_attrs.get("accession"))
    if accession:
        seq_el = SubElement(clade, "sequence")
        SubElement(seq_el, "accession", source="genbank").text = str(accession)

    # Nucleotide mutations
    mutations = branch_attrs.get("mutations", {})
    if not isinstance(mutations, dict):
        raise ConversionError(
            f"'branch_attrs.mutations' of node '{node_label}' must be a JSON "
            f"object, got {type(mutations).__name__}."
        )
    nuc_muts = mutations.get("nuc", [])
    if not isinstance(nuc_muts, list):
        raise ConversionError(
            f"'branch_attrs.mutations.nuc' of node '{node_label}' must be a "
            f"JSON array, got {type(nuc_muts).__name__}."
        )
    if nuc_muts:
        prop = SubElement(
            clade, "property",
            ref="nextclade:nuc_mutations",
            applies_to="clade",
            datatype="xsd:string",
        )
        prop.text = ",".join(str(m) for m in nuc_muts)

    # Amino-acid mutations (one property element per gene)
    for gene, aa_muts in mutations.items():
        if gene == "nuc" or not aa_muts:
            continue
        if not isinstance(aa_muts, list):
            raise ConversionError(
                f"'branch_attrs.mutations.{gene}' of node '{node_label}' must "
                f"be a JSON array, got {type(aa_muts).__name__}."
            )
        prop = SubElement(
            clade, "property",
            ref=f"nextclade:aa_mutations:{gene}",
            applies_to="clade",
            datatype="xsd:string",
        )
        prop.text = ",".join(str(m) for m in aa_muts)

    # Branch labels (e.g., clade, emerging_lineage)
    labels = branch_attrs.get("labels", {})
    if not isinstance(labels, dict):
        raise ConversionError(
            f"'branch_attrs.labels' of node '{node_label}' must be a JSON "
            f"object, got {type(labels).__name__}."
        )
    for label_key, label_val in labels.items():
        if label_val:
            prop = SubElement(
                clade, "property",
                ref=f"nextclade:label:{label_key}",
                applies_to="clade",
                datatype="xsd:string",
            )
            prop.text = str(label_val)

    # Remaining node_attrs → generic <property> elements
    skip_keys = {"div", "num_date", "clade_membership", "country", "division",
                 "location", "accession", "author", "url", "vaccine"}
    for key, val in node_attrs.items():
        if key in skip_keys:
            continue
        scalar = _scalar(val)
        if scalar is not None:
            prop = SubElement(
                clade, "property",
                ref=f"nextclade:{key}",
                applies_to="clade",
                datatype="xsd:string",
            )
            prop.text = str(scalar)

    # Recurse into children
    children = node.get("children", [])
    if not isinstance(children, list):
        raise ConversionError(
            f"'children' of node '{node_label}' must be a JSON array, "
            f"got {type(children).__name__}."
        )
    for i, child in enumerate(children):
        child_path = f"{path}.children[{i}]"
        clade.append(_build_clade(child, parent_div=current_div, path=child_path))

    return clade


def convert(data):
    """Convert a parsed Nextclade/Auspice JSON dict to a PhyloXML ElementTree.

    Args:
        data: dict parsed from a Nextclade Auspice JSON file.

    Returns:
        xml.etree.ElementTree.ElementTree

    Raises:
        ConversionError: if the data is structurally invalid.
    """
    if not isinstance(data, dict):
        raise ConversionError(
            f"Top-level JSON value must be an object, got {type(data).__name__}.\n"
            "Expected either a full Auspice v2 envelope "
            '({"version": "v2", "meta": {...}, "tree": {...}}) '
            "or a bare tree object."
        )

    # Support both full Auspice v2 envelope and bare tree objects
    if "tree" in data:
        tree_root = data["tree"]
        if not isinstance(tree_root, dict):
            raise ConversionError(
                f"The 'tree' key must contain a JSON object, "
                f"got {type(tree_root).__name__}."
            )
        meta = data.get("meta", {})
        if not isinstance(meta, dict):
            raise ConversionError(
                f"The 'meta' key must contain a JSON object, "
                f"got {type(meta).__name__}."
            )
        version = data.get("version", "")
        if version and version != "v2":
            print(
                f"WARNING: Auspice JSON version is '{version}'; "
                "only 'v2' is fully supported. Conversion will proceed but "
                "some fields may be missing or misread.",
                file=sys.stderr,
            )
    elif "name" in data or "node_attrs" in data or "children" in data:
        # Looks like a bare tree node
        tree_root = data
        meta = {}
    else:
        raise ConversionError(
            "Could not identify the JSON format.\n"
            "Expected either:\n"
            '  • An Auspice v2 envelope with a "tree" key, or\n'
            '  • A bare tree node with "name", "node_attrs", or "children" keys.\n'
            "Run 'nextclade run --output-tree nextclade.auspice.json ...' to "
            "produce a supported file."
        )

    ET.register_namespace("", PHYLOXML_NS)
    ET.register_namespace("xsi", XSI_NS)

    phyloxml = Element(
        f"{{{PHYLOXML_NS}}}phyloxml",
        {f"{{{XSI_NS}}}schemaLocation": SCHEMA_LOC},
    )

    rooted = str(meta.get("rooted", "true")).lower()
    phylogeny = SubElement(phyloxml, f"{{{PHYLOXML_NS}}}phylogeny", rooted=rooted)

    title = meta.get("title") or meta.get("name") or "Nextclade tree"
    SubElement(phylogeny, f"{{{PHYLOXML_NS}}}name").text = title

    desc = meta.get("description")
    if desc:
        SubElement(phylogeny, f"{{{PHYLOXML_NS}}}description").text = desc

    root_clade = _build_clade(tree_root, parent_div=0.0, path="tree")
    phylogeny.append(root_clade)

    tree = ElementTree(phyloxml)
    indent(tree, space="  ")
    return tree


def _write_tree(tree, dest, overwrite=False):
    """Write PhyloXML ElementTree to a file path or '-' for stdout."""
    xml_declaration = '<?xml version="1.0" encoding="UTF-8"?>\n'
    if dest == "-":
        sys.stdout.write(xml_declaration)
        tree.write(sys.stdout, encoding="unicode", xml_declaration=False)
        sys.stdout.write("\n")
    else:
        if os.path.exists(dest) and not overwrite:
            raise FileExistsError(
                f"Output file '{dest}' already exists. "
                "Use --overwrite to replace existing files."
            )
        try:
            with open(dest, "w", encoding="utf-8") as fh:
                fh.write(xml_declaration)
                tree.write(fh, encoding="unicode", xml_declaration=False)
                fh.write("\n")
        except OSError as exc:
            raise OSError(exc.strerror) from exc
        print(f"  Written: {dest}", file=sys.stderr)


def _output_path(src_file, src_root, out_root, suffix=".xml"):
    """Compute the output path mirroring src_file's position under src_root."""
    rel = os.path.relpath(src_file, src_root)
    base, ext = os.path.splitext(rel)
    out_rel = base + suffix if ext.lower() == ".json" else rel + suffix
    return os.path.join(out_root, out_rel)


def _load_json(path):
    """Load and return parsed JSON from path, raising descriptive errors."""
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except OSError as exc:
        raise OSError(
            f"Cannot read '{path}': {exc.strerror}. "
            "Check that you have read permission."
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Failed to parse '{path}' as JSON at line {exc.lineno}, "
            f"column {exc.colno}: {exc.msg}. "
            "Make sure this is an Auspice JSON tree file, not a CSV/TSV or Newick file."
        ) from exc


def _convert_one(src, dest, overwrite=False):
    """Load, convert, and write a single file. Returns an error string or None."""
    try:
        data = _load_json(src)
    except (OSError, ValueError) as exc:
        return str(exc)

    try:
        tree = convert(data)
    except ConversionError as exc:
        return str(exc)
    except RecursionError:
        return (
            "Recursion limit exceeded. The tree is extremely deep. "
            "Raise sys.setrecursionlimit() to handle it."
        )

    parent = os.path.dirname(dest)
    if parent:
        os.makedirs(parent, exist_ok=True)

    try:
        _write_tree(tree, dest, overwrite=overwrite)
    except FileExistsError as exc:
        return str(exc)
    except OSError as exc:
        return f"Cannot write '{dest}': {exc}. Check directory permissions."

    return None


def _collect_json_files(directory, recursive):
    """Return sorted list of .json file paths under directory."""
    results = []
    if recursive:
        for dirpath, _, filenames in os.walk(directory):
            for fname in filenames:
                if fname.lower().endswith(".json"):
                    results.append(os.path.join(dirpath, fname))
    else:
        try:
            entries = os.listdir(directory)
        except OSError as exc:
            _die(
                f"Cannot list directory '{directory}': {exc.strerror}.",
                hint="Check that you have read permission for the directory.",
            )
        for fname in entries:
            if fname.lower().endswith(".json"):
                full = os.path.join(directory, fname)
                if os.path.isfile(full):
                    results.append(full)
    return sorted(results)


def _run_directory_mode(src_dir, out_dir, recursive, overwrite=False, suffix=".xml"):
    """Convert all JSON files in src_dir, writing results to out_dir."""
    json_files = _collect_json_files(src_dir, recursive)

    if not json_files:
        _die(
            f"No .json files found in '{src_dir}'"
            + (" (searched recursively)" if recursive else "") + ".",
            hint=(
                "Make sure the directory contains Nextclade Auspice JSON files.\n"
                "       Use --recursive to also search subdirectories."
            ),
        )

    try:
        os.makedirs(out_dir, exist_ok=True)
    except OSError as exc:
        _die(
            f"Cannot create output directory '{out_dir}': {exc.strerror}.",
            hint="Check that the parent directory exists and you have write permission.",
        )

    n_total = len(json_files)
    print(
        f"Converting {n_total} file(s) from '{src_dir}' → '{out_dir}'"
        + (" (recursive)" if recursive else "") + " ...",
        file=sys.stderr,
    )

    failures = []
    for src in json_files:
        dest = _output_path(src, src_dir, out_dir, suffix=suffix)
        err = _convert_one(src, dest, overwrite=overwrite)
        if err:
            print(f"  FAILED:  {src}\n           {err}", file=sys.stderr)
            failures.append(src)

    n_ok = n_total - len(failures)
    print(
        f"\nDone: {n_ok}/{n_total} file(s) converted successfully"
        + (f", {len(failures)} failed." if failures else "."),
        file=sys.stderr,
    )

    if failures:
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        prog="nextclade2phyloxml.py",
        description=(
            "Convert Nextclade phylogenetic tree(s) (Auspice JSON v2) to PhyloXML.\n\n"
            "Accepts a single JSON file (or stdin) or a directory of JSON files.\n"
            "Writes standards-compliant PhyloXML suitable for FigTree, Archaeopteryx,\n"
            "PhyloWidget, Biopython, and similar tools."
        ),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        metavar="INPUT",
        help=(
            "Path to a Nextclade Auspice JSON file or a directory of JSON files. "
            "Use '-' or omit to read a single file from stdin."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default="-",
        metavar="OUTPUT",
        help=(
            "Output file path (single-file mode) or output directory (directory mode). "
            "Use '-' or omit to write to stdout. (default: stdout)"
        ),
    )
    parser.add_argument(
        "-r", "--recursive",
        action="store_true",
        help=(
            "In directory mode, recurse into subdirectories. "
            "The subdirectory structure is mirrored under the output directory."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files. By default, existing files cause an error.",
    )
    parser.add_argument(
        "--suffix",
        default=".xml",
        metavar="EXT",
        help=(
            "File extension for output files in directory mode or when the output "
            "is an existing directory. Must start with '.'. (default: .xml)"
        ),
    )
    args = parser.parse_args()

    # --- Directory mode ---
    if args.input != "-" and os.path.isdir(args.input):
        if args.output == "-":
            _die(
                "An output directory must be specified when the input is a directory.",
                hint=(
                    "Use -o / --output to provide a destination directory:\n"
                    f"       python3 nextclade2phyloxml.py {args.input} -o output_dir/"
                ),
            )
        if os.path.isfile(args.output):
            _die(
                f"Output path '{args.output}' is an existing file, but a directory is required.",
                hint="Provide a directory path (existing or new) as the output.",
            )
        _run_directory_mode(args.input, args.output, args.recursive, overwrite=args.overwrite, suffix=args.suffix)
        return

    # --- Single-file mode ---
    if args.recursive:
        _die(
            "--recursive only applies when the input is a directory.",
            hint=f"Remove --recursive, or pass a directory instead of '{args.input}'.",
        )

    if args.input == "-":
        if sys.stdin.isatty():
            _die(
                "No input file specified and stdin is a terminal.",
                hint=(
                    "Provide a file path:  python3 nextclade2phyloxml.py nextclade.auspice.json\n"
                    "       or a directory:  python3 nextclade2phyloxml.py input_dir/ -o output_dir/\n"
                    "       or pipe input:   cat nextclade.auspice.json | python3 nextclade2phyloxml.py\n"
                    "       or get help:     python3 nextclade2phyloxml.py --help"
                ),
            )
        try:
            data = json.load(sys.stdin)
        except json.JSONDecodeError as exc:
            _die(
                f"Failed to parse JSON from stdin: {exc}",
                hint="Make sure the input is valid JSON (not Newick, NEXUS, or TSV).",
            )
    else:
        if not os.path.exists(args.input):
            _die(
                f"Input not found: '{args.input}'",
                hint=(
                    "Check the file path.  Nextclade produces the tree file as\n"
                    "       'nextclade.auspice.json' when you pass --output-tree."
                ),
            )
        if not os.path.isfile(args.input):
            _die(
                f"'{args.input}' is not a regular file.",
                hint="Provide the path to a .json file, not a directory.",
            )
        try:
            data = _load_json(args.input)
        except (OSError, ValueError) as exc:
            _die(str(exc))

    try:
        tree = convert(data)
    except ConversionError as exc:
        _die(str(exc))
    except RecursionError:
        _die(
            "Recursion limit exceeded while traversing the tree.",
            hint=(
                "The tree is extremely deep.  You can raise Python's recursion limit by\n"
                "       adding 'import sys; sys.setrecursionlimit(10000)' at the top of\n"
                "       the script, or by processing a smaller subtree."
            ),
        )

    if args.output != "-" and os.path.isdir(args.output):
        # Convenience: if user passes an existing directory as output for a
        # single file, write the result into that directory.
        base = os.path.splitext(os.path.basename(args.input))[0] + args.suffix
        dest = os.path.join(args.output, base)
    else:
        dest = args.output

    try:
        _write_tree(tree, dest, overwrite=args.overwrite)
    except FileExistsError as exc:
        _die(str(exc))
    except OSError as exc:
        _die(
            f"Cannot write output file '{dest}': {exc}.",
            hint="Check that the directory exists and you have write permission.",
        )


if __name__ == "__main__":
    main()