#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pure-stdlib script to fetch marine/brackish tardigrade molecular records
(Alaska → Baja California) and write a minimal XLSX workbook.

- No Biopython, pandas, or xlsxwriter required.
- Uses urllib + xml.etree for NCBI E-utilities.
- Creates 3 sheets:
    * Marine_Tardigrades_Individuals
    * Marine_Tardigrades_eDNA
    * Metadata_Notes
"""

from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import xml.etree.ElementTree as ET
import json, time, re, io, zipfile, html, datetime

# ------------- CONFIGURATION -------------
ENTREZ_EMAIL = "mab361@humboldt.edu"  # REQUIRED
SEARCH_TERMS = ["Tardigrada[Organism] AND (COI OR COX1 OR 18S OR 28S OR ITS1 OR ITS2)"]
REGION_FILTER = [
    "Alaska",
    "British Columbia",
    "Washington",
    "Oregon",
    "California",
    "Baja California",
]
MAX_RECORDS = 5000  # hard cap across all search terms
ESARCH_PAGE = 500  # esearch page size
EFETCH_BATCH = 200  # efetch batch size (ids per call)
REQUEST_DELAY = 0.35  # seconds between E-utilities calls (≤3/sec)
OUT_FILE = "Pacific_Marine_Tardigrades_FILLED.xlsx"
TOOL_NAME = "tardigrade_fetcher_stdlib"
# ---- DIAGNOSTICS / FILTER TUNING ----
REGION_FILTER_ENABLED = False           # quick on/off switch
# West Coast box (Alaska → Baja California), used if lat/lon is present
USE_LATLON_BOX = True
LAT_MIN, LAT_MAX = 22.0, 72.0          # degrees North
LON_MIN, LON_MAX = -180.0, -100.0      # degrees West
DEBUG_SHOW_DROPS = 5                   # show up to N examples dropped per batch

DB = "nuccore"              # more precise alias for NCBI nucleotide
EFETCH_BATCH = 100          # keep modest; too large can trigger errors
REQUEST_DELAY = 0.35        # ≤ 3 requests/sec per NCBI guidance

# Diagnostics
WRITE_DEBUG_XML = True      # dump first efetch batch to a file
DEBUG_XML_LIMIT = 1
PRINT_FIRST_N_IDS = 10
REGION_FILTER_ENABLED = False   # keep OFF while testing

# ----------------------------------------

HEADERS = {
    "User-Agent": f"{TOOL_NAME}/1.0 (+https://www.ncbi.nlm.nih.gov/; email:{ENTREZ_EMAIL})",
    "Accept": "*/*",
}

FIELDS = [
    "Species_ID",
    "COI",
    "18S",
    "28S",
    "ITS1",
    "ITS2",
    "Latitude",
    "Longitude",
    "Locality",
    "Depth",
    "Collection_Date",
    "Population_Representation",
    "Data_Type",
    "Citation_DOI",
]

DOI_RE = re.compile(r"(10\.\d{4,9}/[^\s\"<>]+)", re.IGNORECASE)


from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

def http_request(base, params, method="GET", retry=3, timeout=60):
    q = params.copy()
    q["email"] = ENTREZ_EMAIL
    q["tool"] = TOOL_NAME
    data = None
    url = base
    if method.upper() == "GET":
        url = f"{base}?{urlencode(q)}"
    else:
        data = urlencode(q).encode("utf-8")

    last_err = None
    for _ in range(retry):
        try:
            req = Request(url, data=data, headers=HEADERS)
            with urlopen(req, timeout=timeout) as r:
                return r.read()
        except (HTTPError, URLError) as e:
            last_err = e
            time.sleep(1.0)
    raise last_err


def esearch_ids(term, retmax=ESARCH_PAGE, cap=MAX_RECORDS):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    ids = []
    retstart = 0
    while retstart < cap:
        payload = {
            "db": "nucleotide",
            "term": term,
            "retmode": "json",
            "retmax": str(retmax),
            "retstart": str(retstart),
        }
        raw = http_request(base, payload, method="GET")
        data = json.loads(raw.decode("utf-8"))
        lst = data.get("esearchresult", {}).get("idlist", [])
        if not lst:
            break
        ids.extend(lst)
        retstart += retmax
        if len(lst) < retmax:
            break
        time.sleep(REQUEST_DELAY)
        if len(ids) >= cap:
            ids = ids[:cap]
            break
    return ids


def parse_lat_lon(s):
    """Return (lat, lon) as strings; empty if unknown.
    Handles '12.3 N 45.6 W' or '12.3, -45.6' forms."""
    if not s:
        return "", ""
    s = s.strip()
    # decimal comma-separated
    m = re.match(r"^\s*([+-]?\d+(?:\.\d+)?)\s*[,;/]\s*([+-]?\d+(?:\.\d+)?)\s*$", s)
    if m:
        return m.group(1), m.group(2)
    # '12.3 N 45.6 W' (allow commas/spaces)
    m = re.search(
        r"([0-9]+(?:\.[0-9]+)?)\s*([NSns])[^0-9\-+]+([0-9]+(?:\.[0-9]+)?)\s*([EWew])", s
    )
    if m:
        lat = float(m.group(1))
        ns = m.group(2).upper()
        lon = float(m.group(3))
        ew = m.group(4).upper()
        if ns == "S":
            lat = -lat
        if ew == "W":
            lon = -lon
        return f"{lat}", f"{lon}"
    return "", ""


def pick_marker(definition, feature_text):
    """Return a dict of marker -> accession flag to fill."""
    text = f"{definition} {feature_text}".lower()
    flags = {"COI": "", "18S": "", "28S": "", "ITS1": "", "ITS2": ""}
    if ("coi" in text) or ("cox1" in text) or ("cytochrome oxidase subunit i" in text):
        flags["COI"] = "X"
    if ("18s" in text) or ("small subunit" in text) or ("ssu" in text):
        flags["18S"] = "X"
    if ("28s" in text) or ("large subunit" in text) or ("lsu" in text):
        flags["28S"] = "X"
    if "its1" in text:
        flags["ITS1"] = "X"
    if "its2" in text:
        flags["ITS2"] = "X"
    return flags


DOI_RE = re.compile(r"(10\.\d{4,9}/[^\s\"<>]+)", re.IGNORECASE)

def extract_doi_any(ref_node):
    # INSD path
    xp = ref_node.find("./INSDReference_xref")
    if xp is not None:
        for x in xp.findall("./INSDXref"):
            db = (x.findtext("./INSDXref_dbname") or "").lower()
            if db == "doi":
                val = x.findtext("./INSDXref_id") or ""
                if val:
                    return val.strip()

    # GBSeq path
    xp = ref_node.find("./GBReference_xref")
    if xp is not None:
        for x in xp.findall("./GBXref"):
            db = (x.findtext("./GBXref_dbname") or "").lower()
            if db == "doi":
                val = x.findtext("./GBXref_id") or ""
                if val:
                    return val.strip()

    # Fallback: scan text fields for a DOI
    parts = []
    for tag in ("INSDReference_title","INSDReference_journal","INSDReference_remark",
                "GBReference_title","GBReference_journal","GBReference_remark"):
        t = ref_node.findtext(f"./{tag}")
        if t:
            parts.append(t)
    m = DOI_RE.search(" ".join(parts))
    return (m.group(1) if m else "")


def region_match(locality, definition, feature_text, lat_raw):
    """
    Return (is_match: bool, reason: str).
    Match if any REGION_FILTER word appears in locality/definition/notes,
    or lat/lon falls in the West Coast bounding box.
    """
    text = " ".join([locality or "", definition or "", feature_text or ""]).lower()

    # Text match first
    if REGION_FILTER and any(r.lower() in text for r in REGION_FILTER):
        return True, "text"

    # Lat/lon fallback
    if USE_LATLON_BOX and lat_raw:
        lat_s, lon_s = parse_lat_lon(lat_raw)
        if lat_s and lon_s:
            try:
                lat_f = float(lat_s); lon_f = float(lon_s)
                if LAT_MIN <= lat_f <= LAT_MAX and LON_MIN <= lon_f <= LON_MAX:
                    return True, "latlon"
            except ValueError:
                pass

    return False, "none"


def efetch_records(id_batch, batch_index=1):
    """
    Fetch a batch via POST, parse INSDSeq or GBSeq XML, return row dicts.
    Writes the first batch XML to disk if WRITE_DEBUG_XML is True.
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    payload = {
        "db": DB,
        "id": ",".join(id_batch),
        "rettype": "gb",
        "retmode": "xml"
    }
    raw = http_request(base, payload, method="POST")
    time.sleep(REQUEST_DELAY)

    # Debug-dump the first batch XML
    if WRITE_DEBUG_XML and batch_index <= DEBUG_XML_LIMIT:
        with open(f"debug_efetch_batch{batch_index}.xml", "wb") as f:
            f.write(raw)
        print(f"    efetch(): wrote debug_efetch_batch{batch_index}.xml ({len(raw)} bytes)")

    # Detect obvious error payloads
    if b"<eFetchResult>" in raw and b"<ERROR>" in raw:
        try:
            rt = ET.fromstring(raw)
            err = rt.findtext("./ERROR") or ""
        except ET.ParseError:
            err = "(unparseable error XML)"
        print(f"    efetch(): NCBI ERROR -> {err.strip()}")
        return []

    try:
        root = ET.fromstring(raw)
    except ET.ParseError as e:
        snippet = raw[:240].decode("utf-8", "ignore").replace("\n", " ")
        print(f"    efetch(): XML parse error ({e}). First bytes: {snippet}")
        return []

    insd_list = root.findall(".//INSDSeq")
    gbseq_list = root.findall(".//GBSeq")

    rows = []
    if insd_list:
        print(f"    efetch(): format=INSDSeq, nodes={len(insd_list)}")
        rows = _parse_insdseq_nodes(insd_list)
    elif gbseq_list:
        print(f"    efetch(): format=GBSeq, nodes={len(gbseq_list)}")
        rows = _parse_gbseq_nodes(gbseq_list)
    else:
        snippet = raw[:240].decode("utf-8", "ignore").replace("\n"," ")
        print(f"    efetch(): no INSDSeq/GBSeq nodes. First bytes: {snippet}")

    print(f"    efetch(): kept (no region filter) = {len(rows)}")
    return rows


def _parse_insdseq_nodes(nodes):
    out = []
    for insd in nodes:
        accession = (insd.findtext("./INSDSeq_accession-version")
                     or insd.findtext("./INSDSeq_primary-accession") or "")
        organism = insd.findtext("./INSDSeq_organism") or ""
        definition = insd.findtext("./INSDSeq_definition") or ""

        locality = ""; lat_raw = ""; depth = ""; date = ""
        feature_notes = []
        for feat in insd.findall("./INSDSeq_feature-table/INSDFeature"):
            key = feat.findtext("./INSDFeature_key") or ""
            for q in feat.findall("./INSDFeature_quals/INSDQualifier"):
                name = (q.findtext("./INSDQualifier_name") or "").lower()
                val  = q.findtext("./INSDQualifier_value") or ""
                if key == "source":
                    if name == "country": locality = val
                    elif name in ("lat_lon","lat-lon","lat-long","latlong"): lat_raw = val
                    elif name == "depth": depth = val
                    elif name == "collection_date": date = val
                    elif name in {"isolation_source","note"} and val: feature_notes.append(val)
                else:
                    if name in {"gene","product"} and val: feature_notes.append(val)

        doi = ""
        refs = insd.findall("./INSDSeq_references/INSDReference")
        if refs:
            doi = extract_doi_any(refs[0])

        markers = pick_marker(definition, " ".join(feature_notes))
        lat, lon = parse_lat_lon(lat_raw)
        row = {
            "Species_ID": organism,
            "COI": accession if markers["COI"] else "",
            "18S": accession if markers["18S"] else "",
            "28S": accession if markers["28S"] else "",
            "ITS1": accession if markers["ITS1"] else "",
            "ITS2": accession if markers["ITS2"] else "",
            "Latitude": lat, "Longitude": lon,
            "Locality": locality, "Depth": depth, "Collection_Date": date,
            "Population_Representation": "", "Data_Type": "Individual",
            "Citation_DOI": doi
        }
        out.append(row)
    return out


def _parse_gbseq_nodes(nodes):
    out = []
    for gb in nodes:
        accession = (gb.findtext("./GBSeq_accession-version")
                     or gb.findtext("./GBSeq_primary-accession") or "")
        organism = gb.findtext("./GBSeq_organism") or ""
        definition = gb.findtext("./GBSeq_definition") or ""

        locality = ""; lat_raw = ""; depth = ""; date = ""
        feature_notes = []
        for feat in gb.findall("./GBSeq_feature-table/GBFeature"):
            key = feat.findtext("./GBFeature_key") or ""
            for q in feat.findall("./GBFeature_quals/GBQualifier"):
                name = (q.findtext("./GBQualifier_name") or "").lower()
                val  = q.findtext("./GBQualifier_value") or ""
                if key == "source":
                    if name == "country": locality = val
                    elif name in ("lat_lon","lat-lon","lat-long","latlong"): lat_raw = val
                    elif name == "depth": depth = val
                    elif name == "collection_date": date = val
                    elif name in {"isolation_source","note"} and val: feature_notes.append(val)
                else:
                    if name in {"gene","product"} and val: feature_notes.append(val)

        doi = ""
        refs = gb.findall("./GBSeq_references/GBReference")
        if refs:
            doi = extract_doi_any(refs[0])

        markers = pick_marker(definition, " ".join(feature_notes))
        lat, lon = parse_lat_lon(lat_raw)
        row = {
            "Species_ID": organism,
            "COI": accession if markers["COI"] else "",
            "18S": accession if markers["18S"] else "",
            "28S": accession if markers["28S"] else "",
            "ITS1": accession if markers["ITS1"] else "",
            "ITS2": accession if markers["ITS2"] else "",
            "Latitude": lat, "Longitude": lon,
            "Locality": locality, "Depth": depth, "Collection_Date": date,
            "Population_Representation": "", "Data_Type": "Individual",
            "Citation_DOI": doi
        }
        out.append(row)
    return out



# ----------------- Minimal XLSX writer (inline strings; no deps) -----------------


def _col_name(idx1):
    """1-based column index -> Excel column letters"""
    name = ""
    n = idx1
    while n:
        n, r = divmod(n - 1, 26)
        name = chr(65 + r) + name
    return name


def _sheet_xml(rows):
    # rows: List[List[str]]
    buf = io.StringIO()
    buf.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>')
    buf.write(
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
    )
    buf.write("<sheetData>")
    for r_idx, row in enumerate(rows, start=1):
        buf.write(f'<row r="{r_idx}">')
        for c_idx, val in enumerate(row, start=1):
            if val is None:
                continue
            s = str(val)
            if s == "":
                continue
            cref = f"{_col_name(c_idx)}{r_idx}"
            buf.write(
                f'<c r="{cref}" t="inlineStr"><is><t xml:space="preserve">{html.escape(s, quote=True)}</t></is></c>'
            )
        buf.write("</row>")
    buf.write("</sheetData></worksheet>")
    return buf.getvalue().encode("utf-8")


def write_xlsx(path, sheets):
    """
    sheets: Ordered list of tuples (sheet_name, rows)
            rows is a list of lists (first row is header)
    """
    now = datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"

    # Build override lines outside any f-string expression to avoid backslashes inside { ... }.
    overrides_lines = []
    for i, _ in enumerate(sheets, start=1):
        overrides_lines.append(
            f'  <Override PartName="/xl/worksheets/sheet{i}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    overrides_xml = "\n".join(overrides_lines)

    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
        '  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n'
        '  <Default Extension="xml" ContentType="application/xml"/>\n'
        '  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>\n'
        '  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>\n'
        '  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>\n'
        '  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>\n'
        f"{overrides_xml}\n"
        "</Types>"
    ).encode("utf-8")

    rels_root = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>""".encode(
        "utf-8"
    )

    # Workbook relationships (no backslashes in f-expression)
    rel_lines = []
    for i, _ in enumerate(sheets, start=1):
        rel_lines.append(
            f'  <Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
        )
    rel_lines.append(
        '  <Relationship Id="rIdX" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
    )
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
        + "\n".join(rel_lines)
        + "\n</Relationships>"
    ).encode("utf-8")

    # Workbook XML (build <sheet> nodes separately)
    sheet_nodes = []
    for i, (name, _) in enumerate(sheets, start=1):
        sheet_nodes.append(
            f'    <sheet name="{html.escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
        )
    wb_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"\n'
        '          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">\n'
        "  <sheets>\n" + "\n".join(sheet_nodes) + "\n  </sheets>\n</workbook>"
    ).encode("utf-8")

    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font/></fonts>
  <fills count="1"><fill/></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
</styleSheet>""".encode(
        "utf-8"
    )

    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"\n'
        '  xmlns:dc="http://purl.org/dc/elements/1.1/"\n'
        '  xmlns:dcterms="http://purl.org/dc/terms/"\n'
        '  xmlns:dcmitype="http://purl.org/dc/dcmitype/"\n'
        '  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
        f"  <dc:creator>{html.escape(TOOL_NAME)}</dc:creator>\n"
        f"  <cp:lastModifiedBy>{html.escape(TOOL_NAME)}</cp:lastModifiedBy>\n"
        f'  <dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>\n'
        f'  <dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>\n'
        "</cp:coreProperties>"
    ).encode("utf-8")

    app = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
  xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Python (stdlib)</Application>
</Properties>""".encode(
        "utf-8"
    )

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels_root)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/workbook.xml", wb_xml)
        z.writestr("xl/styles.xml", styles)
        z.writestr("docProps/core.xml", core)
        z.writestr("docProps/app.xml", app)
        for i, (_, rows) in enumerate(sheets, start=1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", _sheet_xml(rows))
def main():
    print("Starting tardigrade fetch...")

    all_rows = []
    if not SEARCH_TERMS:
        raise SystemExit("SEARCH_TERMS is empty.")

    for idx, term in enumerate(SEARCH_TERMS, 1):
        print(f"[search {idx}/{len(SEARCH_TERMS)}] {term}")
        ids = esearch_ids(term, retmax=ESARCH_PAGE, cap=MAX_RECORDS)
        print(f"  -> {len(ids)} IDs")

        if not ids:
            continue

        total_batches = (len(ids) + EFETCH_BATCH - 1) // EFETCH_BATCH
        for b, start in enumerate(range(0, len(ids), EFETCH_BATCH), 1):
            batch = ids[start:start+EFETCH_BATCH]
            print(f"  fetching batch {b}/{total_batches} ({len(batch)} IDs)")
            rows = efetch_records(batch)
            print(f"    kept {len(rows)} rows")
            all_rows.extend(rows)

    print(f"\nTotal kept after filtering: {len(all_rows)}")
    print("Writing Excel...")

    individuals = [FIELDS] + [[r.get(k, "") for k in FIELDS] for r in all_rows]
    edna = [FIELDS]
    meta = [
        ["Notes"],
        ["Dataset compiled via NCBI E-utilities (stdlib)"],
        ["Regions: Alaska → Baja California"],
        ["Data fields left blank where unavailable"],
    ]

    write_xlsx(
        OUT_FILE,
        [
            ("Marine_Tardigrades_Individuals", individuals),
            ("Marine_Tardigrades_eDNA", edna),
            ("Metadata_Notes", meta),
        ]
    )

    print(f"Saved: {OUT_FILE}")
    print("Done.")

if __name__ == "__main__":
    main()

