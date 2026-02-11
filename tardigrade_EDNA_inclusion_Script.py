#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch tardigrade nucleotide records from NCBI, write XLSX (Individuals + eDNA).

- Base Python only: urllib + xml.etree + zipfile (no Biopython/pandas/xlsxwriter).
- Uses POST for efetch (avoids long-URL issues).
- Parses both INSDSeq and GBSeq XML.
- Baseline defaults: no region filter, no eDNA split -> large Individuals sheet.
"""

import time, re, json, io, zipfile, html, datetime
from collections import defaultdict
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import xml.etree.ElementTree as ET

# ========================== CONFIGURATION ==========================

ENTREZ_EMAIL = "mab361@humboldt.edu"   # REQUIRED (NCBI policy)
TOOL_NAME = "tardigrade_fetcher_stdlib"

# Query terms
SEARCH_TERMS = [
    "Tardigrada[Organism] AND (COI OR COX1 OR 18S OR 28S OR ITS1 OR ITS2)"
    # You can add more terms here (e.g., environmental keywords), the script deduplicates IDs.
]

# NCBI / pagination
DB = "nucleotide"      # broader than "nuccore" -> more hits, closer to your earlier run
MAX_RECORDS = 5000     # cap across all search terms
ESARCH_PAGE = 500      # esearch page size
EFETCH_BATCH = 100     # ids per efetch POST
REQUEST_DELAY = 0.35   # ~3 requests/sec (NCBI-friendly)

# Diagnostics toggles (defined here to avoid NameError)
PRINT_FIRST_N_IDS = 10
WRITE_DEBUG_XML = False    # True dumps first batch XML to "debug_efetch_batch1.xml"
DEBUG_XML_LIMIT = 1

# Output
OUT_FILE = "Pacific_Marine_Tardigrades_FILLED.xlsx"

# Region filter (OFF by default for baseline big counts)
REGION_FILTER_ENABLED = True
USE_LATLON_BOX = True
LAT_MIN, LAT_MAX = 22.0, 72.0
LON_MIN, LON_MAX = -180.0, -100.0
REGION_ALLOW_PACIFIC_GENERAL = True
EDNA_REGION_RELAXED = True
DEBUG_SHOW_DROPS = 12


REGION_TEXT_TOKENS = [
  # Core jurisdictions
  "alaska","british columbia","washington","oregon","california",
  "baja california","baja california sur","baja california norte",
  "usa: alaska","usa: washington","usa: oregon","usa: california",
  "canada: british columbia","mexico: baja california","mexico: baja california sur","mexico: baja california norte",
  "us: ak","us: wa","us: or","us: ca","usa: ak","usa: wa","usa: or","usa: ca",
  "canada: bc","ca: bc","mexico: bc","mexico: bcs","mexico: bcn","bcs","bcn",

  # General Pacific signals
  "pacific ocean","north pacific","eastern pacific","northeast pacific","california current","west coast",

  # Alaska
  "gulf of alaska","aleutian","prince william sound","kodiak","kenai","kachemak",
  "southeast alaska","sitka","ketchikan","juneau",

  # British Columbia
  "vancouver island","haida gwaii","queen charlotte","inside passage","strait of georgia",
  "salish sea","juan de fuca","barkley sound","tofino","ucluelet",

  # Washington
  "puget sound","san juan islands","olympic coast","grays harbor","willapa bay",

  # Oregon
  "columbia river","tillamook","coos bay","yaquina","newport","brookings","oregon coast",

  # California (north → south)
  "humboldt","trinidad head","mendocino","bodega bay","point reyes","san francisco bay","half moon bay",
  "monterey bay","moss landing","carmel","point lobos","big sur","morro bay","avila","pismo",
  "santa barbara","goleta","santa barbara channel","channel islands","anacapa","santa cruz island",
  "point conception","santa catalina","catalina island","san pedro","long beach","palos verdes",
  "redondo","santa monica bay","malibu","ventura","la jolla","scripps pier","mission bay",
  "san diego","san diego bay","imperial beach",

  # Baja California (north & sur) + Gulf of California
  "ensenada","rosarito","san quintín","bahía san quintín","bahia san quintin",
  "bahía de los ángeles","bahia de los angeles","isla de cedros","bahía tortugas","bahia tortugas",
  "guerrero negro","laguna ojo de liebre","vizcaíno","laguna san ignacio","bahía magdalena","bahia magdalena",
  "loreto","bahía concepción","bahia concepcion","mulegé","mulege","la paz",
  "cabo san lucas","san josé del cabo","san jose del cabo","todos santos",
  "gulf of california","sea of cortez","mar de cortés","mar de cortes"
]


# eDNA split (OFF by default for baseline big counts)
CLASSIFY_EDNA = True
EDNA_KEYWORDS = [
    "edna", "environmental dna", "environmental sample", "metabarcoding",
    "metagenom", "amplicon", "bulk sample", "seawater", "sea water",
    "water filter", "plankton", "sediment", "marine sediment", "biofilm"
]

# Population representation for Individuals
COMPUTE_POPULATION_REP = True
GROUP_ROUND = 2   # decimal places for lat/lon when grouping

# Excel columns (both sheets share the same schema)
FIELDS = [
    "Species_ID", "COI", "18S", "28S", "ITS1", "ITS2",
    "Latitude", "Longitude", "Locality", "Depth",
    "Collection_Date", "Population_Representation", "Data_Type", "Citation_DOI"
]

# ==================================================================

HEADERS = {
    "User-Agent": f"{TOOL_NAME}/1.0 (+https://www.ncbi.nlm.nih.gov/; email:{ENTREZ_EMAIL})",
    "Accept": "*/*"
}

DOI_RE = re.compile(r"(10\.\d{4,9}/[^\s\"<>]+)", re.IGNORECASE)

# --------------------------- HTTP helper ---------------------------

def http_request(base, params, method="GET", retry=3, timeout=60):
    """GET/POST wrapper with basic retry and NCBI-required params."""
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

# ------------------------------ ESearch ----------------------------

def esearch_ids(term, retmax=ESARCH_PAGE, cap=MAX_RECORDS):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    ids, retstart = [], 0
    while retstart < cap:
        payload = {
            "db": DB,
            "term": term,
            "retmode": "json",
            "retmax": str(retmax),
            "retstart": str(retstart)
        }
        raw = http_request(base, payload, method="GET")
        data = json.loads(raw.decode("utf-8"))
        lst = data.get("esearchresult", {}).get("idlist", [])
        if retstart == 0:
            total_count = int(data.get("esearchresult", {}).get("count", "0") or 0)
            n_show = min(PRINT_FIRST_N_IDS, len(lst))
            print(f"  esearch total available: {total_count}")
            print("  esearch sample IDs:", lst[:n_show])
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

# -------------------------- Parse helpers --------------------------

def parse_lat_lon(s):
    if not s:
        return "", ""
    s = s.strip()
    m = re.match(r"^\s*([+-]?\d+(?:\.\d+)?)\s*[,;/]\s*([+-]?\d+(?:\.\d+)?)\s*$", s)
    if m:
        return m.group(1), m.group(2)
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*([NSns])[^0-9\-+]+([0-9]+(?:\.[0-9]+)?)\s*([EWew])", s)
    if m:
        lat = float(m.group(1)); ns = m.group(2).upper()
        lon = float(m.group(3)); ew = m.group(4).upper()
        if ns == "S": lat = -lat
        if ew == "W": lon = -lon
        return f"{lat}", f"{lon}"
    return "", ""

def pick_marker(definition, feature_text):
    text = f"{definition or ''} {feature_text or ''}".lower()
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

def extract_doi_any(ref_node):
    xp = ref_node.find("./INSDReference_xref")
    if xp is not None:
        for x in xp.findall("./INSDXref"):
            db = (x.findtext("./INSDXref_dbname") or "").lower()
            if db == "doi":
                val = x.findtext("./INSDXref_id") or ""
                if val:
                    return val.strip()
    xp = ref_node.find("./GBReference_xref")
    if xp is not None:
        for x in xp.findall("./GBXref"):
            db = (x.findtext("./GBXref_dbname") or "").lower()
            if db == "doi":
                val = x.findtext("./GBXref_id") or ""
                if val:
                    return val.strip()
    parts = []
    for tag in ("INSDReference_title","INSDReference_journal","INSDReference_remark",
                "GBReference_title","GBReference_journal","GBReference_remark"):
        t = ref_node.findtext(f"./{tag}")
        if t:
            parts.append(t)
    m = DOI_RE.search(" ".join(parts))
    return (m.group(1) if m else "")

def detect_edna(locality, definition, feature_text, qualifier_names):
    if not CLASSIFY_EDNA:
        return False, "off"
    blob = f" {(locality or '')} {(definition or '')} {(feature_text or '')} ".lower()
    if "environmental_sample" in qualifier_names or "metagenomic" in qualifier_names:
        return True, "qualifier"
    for k in EDNA_KEYWORDS:
        if k in blob:
            return True, "keyword"
    return False, ""

def region_match(locality, definition, feature_text, lat_raw):
    if not REGION_FILTER_ENABLED:
        return True, "disabled"

    blob = f" {(locality or '')} {(definition or '')} {(feature_text or '')} ".lower()

    # Pacific general pass
    if REGION_ALLOW_PACIFIC_GENERAL and any(tok in blob for tok in (
        "pacific ocean","north pacific","eastern pacific","northeast pacific","california current","west coast"
    )):
        return True, "pacific"

    # Token match
    for tok in REGION_TEXT_TOKENS:
        if tok in blob:
            return True, "text"

    # Lat/lon box fallback
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

# ------------------------- EFetch + parsers ------------------------

def efetch_records(id_batch, batch_index=1):
    """Fetch via POST; parse INSDSeq/GBSeq; return (individual_rows, edna_rows) with region applied."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    payload = {
        "db": DB,
        "id": ",".join(id_batch),
        "rettype": "gb",
        "retmode": "xml"
    }
    raw = http_request(base, payload, method="POST")
    time.sleep(REQUEST_DELAY)

    if WRITE_DEBUG_XML and batch_index <= DEBUG_XML_LIMIT:
        with open(f"debug_efetch_batch{batch_index}.xml", "wb") as f:
            f.write(raw)
        print(f"    efetch(): wrote debug_efetch_batch{batch_index}.xml ({len(raw)} bytes)")

    try:
        root = ET.fromstring(raw)
    except ET.ParseError as e:
        snippet = raw[:240].decode("utf-8", "ignore").replace("\n"," ")
        print("    efetch(): XML parse error:", e, "| first bytes:", snippet)
        return [], []

    insd_list = root.findall(".//INSDSeq")
    gbseq_list = root.findall(".//GBSeq")
    if insd_list:
        print(f"    efetch(): format=INSDSeq, nodes={len(insd_list)}")
        parsed = _parse_insdseq_nodes(insd_list)
    elif gbseq_list:
        print(f"    efetch(): format=GBSeq, nodes={len(gbseq_list)}")
        parsed = _parse_gbseq_nodes(gbseq_list)
    else:
        print("    efetch(): no INSDSeq/GBSeq nodes")
        return [], []

    ind_rows, edna_rows = _build_rows_from_parsed(parsed)
    print(f"    efetch(): kept (after region) -> individuals={len(ind_rows)}, eDNA={len(edna_rows)}")
    return ind_rows, edna_rows


def _parse_insdseq_nodes(nodes):
    parsed = []
    for insd in nodes:
        rec = {
            "accession": (insd.findtext("./INSDSeq_accession-version")
                          or insd.findtext("./INSDSeq_primary-accession") or ""),
            "organism": insd.findtext("./INSDSeq_organism") or "",
            "definition": insd.findtext("./INSDSeq_definition") or "",
            "locality": "", "lat_raw": "", "depth": "", "date": "",
            "feature_notes": [], "qual_names": set(), "doi": ""
        }
        for feat in insd.findall("./INSDSeq_feature-table/INSDFeature"):
            key = feat.findtext("./INSDFeature_key") or ""
            for q in feat.findall("./INSDFeature_quals/INSDQualifier"):
                name = (q.findtext("./INSDQualifier_name") or "").lower()
                val  = q.findtext("./INSDQualifier_value") or ""
                if key == "source":
                    rec["qual_names"].add(name)
                    if name == "country": rec["locality"] = val
                    elif name in ("lat_lon","lat-lon","lat-long","latlong"): rec["lat_raw"] = val
                    elif name == "depth": rec["depth"] = val
                    elif name == "collection_date": rec["date"] = val
                    elif name in {"isolation_source","note"} and val: rec["feature_notes"].append(val)
                else:
                    if name in {"gene","product"} and val: rec["feature_notes"].append(val)
        refs = insd.findall("./INSDSeq_references/INSDReference")
        if refs:
            rec["doi"] = extract_doi_any(refs[0])
        parsed.append(rec)
    return parsed

def _parse_gbseq_nodes(nodes):
    parsed = []
    for gb in nodes:
        rec = {
            "accession": (gb.findtext("./GBSeq_accession-version")
                          or gb.findtext("./GBSeq_primary-accession") or ""),
            "organism": gb.findtext("./GBSeq_organism") or "",
            "definition": gb.findtext("./GBSeq_definition") or "",
            "locality": "", "lat_raw": "", "depth": "", "date": "",
            "feature_notes": [], "qual_names": set(), "doi": ""
        }
        for feat in gb.findall("./GBSeq_feature-table/GBFeature"):
            key = feat.findtext("./GBFeature_key") or ""
            for q in feat.findall("./GBFeature_quals/GBQualifier"):
                name = (q.findtext("./GBQualifier_name") or "").lower()
                val  = q.findtext("./GBQualifier_value") or ""
                if key == "source":
                    rec["qual_names"].add(name)
                    if name == "country": rec["locality"] = val
                    elif name in ("lat_lon","lat-lon","lat-long","latlong"): rec["lat_raw"] = val
                    elif name == "depth": rec["depth"] = val
                    elif name == "collection_date": rec["date"] = val
                    elif name in {"isolation_source","note"} and val: rec["feature_notes"].append(val)
                else:
                    if name in {"gene","product"} and val: rec["feature_notes"].append(val)
        refs = gb.findall("./GBSeq_references/GBReference")
        if refs:
            rec["doi"] = extract_doi_any(refs[0])
        parsed.append(rec)
    return parsed

def _build_rows_from_parsed(parsed):
    """Return (individual_rows, eDNA_rows). Region filter applied here."""
    ind_rows, edna_rows = [], []

    # Counters (Part D) + drops (Part E)
    reason_counts = {"text": 0, "latlon": 0, "pacific": 0}
    drops = 0
    examples = []

    # Read flags safely
    classify_edna = bool(globals().get("CLASSIFY_EDNA", False))
    edna_relaxed  = bool(globals().get("EDNA_REGION_RELAXED", False))
    region_on     = bool(globals().get("REGION_FILTER_ENABLED", False))
    show_drops    = int(globals().get("DEBUG_SHOW_DROPS", 12))

    for rec in parsed:
        feat_text = " ".join(rec.get("feature_notes", []))

        # Region decision
        ok_region, reason = region_match(
            rec.get("locality",""), rec.get("definition",""), feat_text, rec.get("lat_raw","")
        )
        if ok_region and reason in reason_counts:
            reason_counts[reason] += 1

        # eDNA vs Individual classification (string inspect + qualifiers)
        is_edna = False
        try:
            is_edna, _ = detect_edna(rec.get("locality",""), rec.get("definition",""), feat_text, rec.get("qual_names", set()))
        except NameError:
            pass  # detect_edna not present in this variant

        # Skip early if region fails and we’re strict (handled differently for eDNA below)
        if not ok_region and not is_edna and region_on:
            drops += 1
            if len(examples) < show_drops:
                examples.append(
                    f"drop: {rec.get('accession','?')} | loc='{rec.get('locality','')}' | def='{(rec.get('definition') or '')[:80]}'"
                )
            continue

        # Build row
        markers = pick_marker(rec.get("definition",""), feat_text)
        lat, lon = parse_lat_lon(rec.get("lat_raw",""))
        row = {
            "Species_ID": rec.get("organism",""),
            "COI":  rec.get("accession","") if markers["COI"]  else "",
            "18S":  rec.get("accession","") if markers["18S"]  else "",
            "28S":  rec.get("accession","") if markers["28S"]  else "",
            "ITS1": rec.get("accession","") if markers["ITS1"] else "",
            "ITS2": rec.get("accession","") if markers["ITS2"] else "",
            "Latitude": lat,
            "Longitude": lon,
            "Locality": rec.get("locality",""),
            "Depth": rec.get("depth",""),
            "Collection_Date": rec.get("date",""),
            "Population_Representation": "",
            "Data_Type": "eDNA" if is_edna and classify_edna else "Individual",
            "Citation_DOI": rec.get("doi","")
        }

        # Route row
        if is_edna and classify_edna:
            # keep if region passes OR region filter is off OR relaxed for eDNA
            if (not region_on) or ok_region or edna_relaxed:
                edna_rows.append(row)
            else:
                drops += 1
                if len(examples) < show_drops:
                    examples.append(
                        f"drop(eDNA): {rec.get('accession','?')} | loc='{rec.get('locality','')}' | def='{(rec.get('definition') or '')[:80]}'"
                    )
            continue

        # Individuals: keep if region passes or filter is off
        if (not region_on) or ok_region:
            ind_rows.append(row)

    # Batch summary
    if region_on:
        print(
            "    region reasons this batch: "
            f"text={reason_counts.get('text',0)}, "
            f"latlon={reason_counts.get('latlon',0)}, "
            f"pacific={reason_counts.get('pacific',0)}"
        )
        print(f"    region drops in this batch: {drops}")
        for line in examples:
            print("    " + line)

    return ind_rows, edna_rows



# ------------------ Population representation ---------------------

def add_population_rep(rows):
    """For Individuals, write group size if multiple records share species+locality+date+rounded coords."""
    if not COMPUTE_POPULATION_REP or not rows:
        return rows
    counts = defaultdict(int)
    keys = []
    for r in rows:
        lat = r.get("Latitude") or ""
        lon = r.get("Longitude") or ""
        try:
            lat_f = round(float(lat), GROUP_ROUND)
            lon_f = round(float(lon), GROUP_ROUND)
        except Exception:
            lat_f, lon_f = lat, lon
        key = (r.get("Species_ID",""), r.get("Locality",""), r.get("Collection_Date",""), lat_f, lon_f)
        keys.append(key)
        counts[key] += 1
    for r, k in zip(rows, keys):
        n = counts[k]
        r["Population_Representation"] = str(n) if n > 1 else ""
    return rows

# -------------------------- Minimal XLSX ---------------------------

def _col_name(idx1):
    name = ""
    n = idx1
    while n:
        n, r = divmod(n - 1, 26)
        name = chr(65 + r) + name
    return name

def _sheet_xml(rows):
    buf = io.StringIO()
    buf.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>')
    buf.write('<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">')
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
            buf.write(f'<c r="{cref}" t="inlineStr"><is><t xml:space="preserve">{html.escape(s, quote=True)}</t></is></c>')
        buf.write("</row>")
    buf.write("</sheetData></worksheet>")
    return buf.getvalue().encode("utf-8")

def write_xlsx(path, sheets):
    """sheets: list of (sheet_name, rows) where rows is list-of-lists (row 0 is header)."""
    now = datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"

    # Build overrides for sheet parts (avoid backslashes in f-expr)
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
        '  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties"/>\n'
        '  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties"/>\n'
        f'{overrides_xml}\n'
        '</Types>'
    ).encode("utf-8")

    rels_root = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>""".encode("utf-8")

    # workbook rels
    rel_lines = []
    for i, _ in enumerate(sheets, start=1):
        rel_lines.append(
            f'  <Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
        )
    rel_lines.append('  <Relationship Id="rIdX" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>')
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
        + "\n".join(rel_lines) + '\n</Relationships>'
    ).encode("utf-8")

    # workbook xml
    sheet_nodes = []
    for i, (name, _) in enumerate(sheets, start=1):
        sheet_nodes.append(
            f'    <sheet name="{html.escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
        )
    wb_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"\n'
        '          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">\n'
        '  <sheets>\n' +
        "\n".join(sheet_nodes) +
        '\n  </sheets>\n</workbook>'
    ).encode("utf-8")

    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font/></fonts>
  <fills count="1"><fill/></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
</styleSheet>""".encode("utf-8")

    core = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"\n'
        '  xmlns:dc="http://purl.org/dc/elements/1.1/"\n'
        '  xmlns:dcterms="http://purl.org/dc/terms/"\n'
        '  xmlns:dcmitype="http://purl.org/dc/dcmitype/"\n'
        '  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
        f'  <dc:creator>{html.escape(TOOL_NAME)}</dc:creator>\n'
        f'  <cp:lastModifiedBy>{html.escape(TOOL_NAME)}</cp:lastModifiedBy>\n'
        f'  <dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>\n'
        f'  <dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>\n'
        '</cp:coreProperties>'
    ).encode("utf-8")

    app = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
  xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Python (stdlib)</Application>
</Properties>""".encode("utf-8")

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

# ------------------------------- main ------------------------------

def main():
    if not ENTREZ_EMAIL or "@" not in ENTREZ_EMAIL:
        raise SystemExit("Set ENTREZ_EMAIL at the top of the script before running.")

    print("Config:",
          f"DB={DB}; MAX_RECORDS={MAX_RECORDS}; ESARCH_PAGE={ESARCH_PAGE}; EFETCH_BATCH={EFETCH_BATCH};",
          f"RegionFilter={REGION_FILTER_ENABLED}; eDNA={CLASSIFY_EDNA}")

    all_ids, seen = [], set()
    for idx, term in enumerate(SEARCH_TERMS, 1):
        print(f"[search {idx}/{len(SEARCH_TERMS)}] {term}")
        ids = esearch_ids(term, retmax=ESARCH_PAGE, cap=MAX_RECORDS)
        print(f"  -> got {len(ids)} IDs from esearch")
        for i in ids:
            if i not in seen:
                seen.add(i); all_ids.append(i)

    print(f"Total unique IDs: {len(all_ids)}")
    if not all_ids:
        print("No IDs found. Exiting.")
        return

    individuals, edna = [], []
    total_batches = (len(all_ids) + EFETCH_BATCH - 1) // EFETCH_BATCH
    for b in range(total_batches):
        start = b * EFETCH_BATCH
        batch = all_ids[start:start+EFETCH_BATCH]
        print(f"Fetching batch {b+1}/{total_batches} ({len(batch)} IDs)…")
        ind_rows, edna_rows = efetch_records(batch, batch_index=b+1)
        individuals.extend(ind_rows)
        edna.extend(edna_rows)
        print(f"  cumulative rows so far: {len(individuals) + len(edna)}")

    if COMPUTE_POPULATION_REP:
        individuals = add_population_rep(individuals)

    print(f"Final counts | Individuals: {len(individuals)} | eDNA: {len(edna)}")

    indiv_rows = [FIELDS] + [[r.get(k, "") for k in FIELDS] for r in individuals]
    edna_rows  = [FIELDS] + [[r.get(k, "") for k in FIELDS] for r in edna]
    meta = [
        ["Notes"],
        ["Dataset compiled via NCBI E-utilities (stdlib)"],
        [f"Region filter: {'ON' if REGION_FILTER_ENABLED else 'OFF'}; lat/lon box={USE_LATLON_BOX}"],
        [f"eDNA split: {'ON' if CLASSIFY_EDNA else 'OFF'}"],
        ["Columns: " + ", ".join(FIELDS)],
        [f"Generated UTC: {datetime.datetime.utcnow().isoformat()}Z"]
    ]

    write_xlsx(
        OUT_FILE,
        [
            ("Marine_Tardigrades_Individuals", indiv_rows),
            ("Marine_Tardigrades_eDNA", edna_rows),
            ("Metadata_Notes", meta),
        ]
    )
    print(f"Saved: {OUT_FILE}")
    print("Done.")

if __name__ == "__main__":
    main()
