import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
import re
import glob

def parse_bgc_gbk_files(bgc_dir):
    """
    Fast pure-Python parser for antiSMASH region GenBank files.
    Returns a dictionary mapping locus_tag -> dict of BGC metadata.
    """
    bgc_map = {}
    bgc_path = Path(bgc_dir)
    if not bgc_path.exists():
        return bgc_map

    gbk_files = list(bgc_path.rglob("*.region*.gbk"))
    for gbk_file in gbk_files:
        bgc_id = gbk_file.stem
        strain = gbk_file.parent.name
        bgc_type = "Unknown"
        
        try:
            with open(gbk_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
        except Exception:
            continue

        # Extract BGC product/type from region feature or cand_cluster
        m_region = re.search(r'     region\s+[\d\.<>]+\n(?:      .*\n)*?\s+/product="([^"]+)"', content)
        if m_region:
            bgc_type = m_region.group(1).strip()
        else:
            m_cand = re.search(r'/product="([^"]+)"', content)
            if m_cand:
                bgc_type = m_cand.group(1).strip()

        # Extract CDS features
        cds_blocks = re.findall(r'     CDS\s+[\w\(\)\d\.<>,]+\n((?:      .*\n)+)', content)
        for block in cds_blocks:
            m_lt = re.search(r'/locus_tag="([^"]+)"', block)
            if m_lt:
                lt = m_lt.group(1).strip()
                m_gene = re.search(r'/gene="([^"]+)"', block)
                gene = m_gene.group(1).strip() if m_gene else ""
                m_prod = re.search(r'/product="([^"]+)"', block)
                prod = m_prod.group(1).strip() if m_prod else ""
                bgc_map[lt] = {
                    'bgc_id': bgc_id,
                    'bgc_type': bgc_type,
                    'strain': strain,
                    'gene_name': gene,
                    'product': prod
                }

    return bgc_map

def inject_cazyme_and_bgc_to_gexf(gexf_in, pangenome_csv, cazyme_csv, gexf_out, bgc_dir=None, overlap_csv_out=None):
    """
    Injects 4-way classification (BGC_only, CAZyme_only, BGC_and_CAZyme, Other)
    and native Gephi <viz:color> into a PPanGGOLiN GEXF file.
    Also exports an overlap report table.
    """
    gexf_in = Path(gexf_in)
    gexf_out = Path(gexf_out)
    gexf_out.parent.mkdir(parents=True, exist_ok=True)

    # Resolve BGC directory if not provided
    if not bgc_dir or not Path(bgc_dir).exists():
        # Attempt standard bgcflow locations
        possible_bgc_dirs = [
            gexf_in.parents[3] / "antismash" / "8.0.4",
            Path("data/interim/bgcs") / gexf_in.parents[3].name / "8.0.4",
            Path("data/interim/antismash/8.0.4")
        ]
        for candidate in possible_bgc_dirs:
            if candidate.exists():
                bgc_dir = candidate
                break

    # Resolve overlap CSV output path
    if not overlap_csv_out:
        overlap_csv_out = gexf_in.parents[3] / "tables" / "df_bgc_cazyme_overlap.csv"
    overlap_csv_out = Path(overlap_csv_out)
    overlap_csv_out.parent.mkdir(parents=True, exist_ok=True)

    # 1. Load CAZyme annotations
    cazyme_map = {}
    if Path(cazyme_csv).exists() and Path(cazyme_csv).stat().st_size > 0:
        try:
            df_cazyme = pd.read_csv(cazyme_csv, sep="\t" if str(cazyme_csv).endswith(".tsv") else ",")
            if "gene_id" in df_cazyme.columns and "cazyme_family" in df_cazyme.columns:
                for _, row in df_cazyme.iterrows():
                    gid = str(row["gene_id"]).strip()
                    caz = str(row["cazyme_family"]).strip()
                    cls = str(row.get("cazyme_class", "CAZyme")).strip()
                    eval_val = str(row.get("evalue", ""))
                    score_val = str(row.get("score", ""))
                    if gid not in cazyme_map:
                        cazyme_map[gid] = {
                            'families': [caz],
                            'classes': [cls],
                            'evalue': eval_val,
                            'score': score_val
                        }
                    else:
                        if caz not in cazyme_map[gid]['families']:
                            cazyme_map[gid]['families'].append(caz)
                        if cls not in cazyme_map[gid]['classes']:
                            cazyme_map[gid]['classes'].append(cls)
        except Exception as e:
            print(f"Warning reading CAZyme CSV: {e}", file=sys.stderr)

    # 2. Parse BGC region GenBanks
    bgc_map = {}
    if bgc_dir:
        bgc_map = parse_bgc_gbk_files(bgc_dir)
    print(f"Loaded {len(cazyme_map)} CAZyme gene annotations and {len(bgc_map)} BGC CDS features.")

    # 3. Map homologous gene families via PPanGGOLiN matrix
    pangenome_path = Path(pangenome_csv)
    pangenome_file = None
    if pangenome_path.is_dir():
        for candidate in ["matrix.csv", "gene_presence_absence.csv"]:
            if (pangenome_path / candidate).exists():
                pangenome_file = pangenome_path / candidate
                break
        if not pangenome_file:
            csv_files = list(pangenome_path.glob("*.csv"))
            if csv_files:
                pangenome_file = csv_files[0]
    elif pangenome_path.is_file():
        pangenome_file = pangenome_path

    family_meta = {}
    gene_to_family = {}
    if pangenome_file and pangenome_file.exists() and pangenome_file.stat().st_size > 0:
        try:
            df_pan = pd.read_csv(pangenome_file, dtype=str, low_memory=False)
            gene_col = "Gene" if "Gene" in df_pan.columns else df_pan.columns[0]
            meta_cols = {"Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", "Avg sequences per isolate", "Genome Fragment"}
            sample_cols = [c for c in df_pan.columns if c not in meta_cols]

            for _, row in df_pan.iterrows():
                fam = str(row[gene_col]).strip()
                fam_caz_fams = set()
                fam_caz_classes = set()
                fam_bgc_types = set()
                fam_bgc_ids = set()

                for sc in sample_cols:
                    val = str(row.get(sc, ""))
                    if val and val != "nan":
                        for token in val.replace("\t", " ").replace(";", " ").split():
                            t = token.strip()
                            gene_to_family[t] = fam
                            if t in cazyme_map:
                                fam_caz_fams.update(cazyme_map[t]['families'])
                                fam_caz_classes.update(cazyme_map[t]['classes'])
                            if t in bgc_map:
                                fam_bgc_types.add(bgc_map[t]['bgc_type'])
                                fam_bgc_ids.add(bgc_map[t]['bgc_id'])

                is_caz = len(fam_caz_fams) > 0
                is_bgc = len(fam_bgc_types) > 0

                if is_bgc and is_caz:
                    status = "BGC_and_CAZyme"
                elif is_bgc:
                    status = "BGC_only"
                elif is_caz:
                    status = "CAZyme_only"
                else:
                    status = "Other"

                family_meta[fam] = {
                    'status': status,
                    'is_bgc': "True" if is_bgc else "False",
                    'bgc_type': ";".join(sorted(fam_bgc_types)) if fam_bgc_types else "None",
                    'bgc_id': ";".join(sorted(fam_bgc_ids)) if fam_bgc_ids else "None",
                    'is_cazyme': "True" if is_caz else "False",
                    'cazyme_family': ";".join(sorted(fam_caz_fams)) if fam_caz_fams else "None",
                    'cazyme_class': ";".join(sorted(fam_caz_classes)) if fam_caz_classes else "Non-CAZyme"
                }
        except Exception as e:
            print(f"Warning parsing pangenome matrix: {e}", file=sys.stderr)

    # 4. Generate BGC + CAZyme overlap table
    overlap_rows = []
    for lt, binfo in bgc_map.items():
        if lt in cazyme_map:
            cinfo = cazyme_map[lt]
            fam = gene_to_family.get(lt, "Unclustered")
            overlap_rows.append({
                'strain': binfo['strain'],
                'locus_tag': lt,
                'gene_name': binfo['gene_name'],
                'product': binfo['product'],
                'bgc_id': binfo['bgc_id'],
                'bgc_type': binfo['bgc_type'],
                'cazyme_family': ";".join(cinfo['families']),
                'cazyme_class': ";".join(cinfo['classes']),
                'evalue': cinfo['evalue'],
                'score': cinfo['score'],
                'pangenome_family': fam
            })

    df_overlap = pd.DataFrame(overlap_rows)
    if not df_overlap.empty:
        df_overlap.sort_values(by=['strain', 'bgc_type', 'cazyme_family'], inplace=True)
    df_overlap.to_csv(str(overlap_csv_out), index=False)
    print(f"Saved {len(df_overlap)} overlapping genes (BGC + CAZyme) to {overlap_csv_out}")

    # 5. Parse and Enrich GEXF with <viz:color> and node attributes
    with open(str(gexf_in), 'r', encoding='utf-8') as f:
        gexf_content = f.read()

    # Sanitize unescaped XML characters (< and >) in attribute values
    gexf_content = re.sub(r'value="([^"]*)"', lambda m: 'value="' + m.group(1).replace('<', '&lt;').replace('>', '&gt;') + '"', gexf_content)
    root = ET.fromstring(gexf_content)
    tree = ET.ElementTree(root)

    uri = root.tag.split('}')[0][1:] if root.tag.startswith('{') else ""
    viz_uri = f"{uri}/viz" if uri else "http://www.gexf.net/1.2draft/viz"
    
    if uri:
        ET.register_namespace('', uri)
        ET.register_namespace('viz', viz_uri)
        ns = {'g': uri, 'viz': viz_uri}
        tag_prefix = f"{{{uri}}}"
        viz_prefix = f"{{{viz_uri}}}"
    else:
        ns = {}
        tag_prefix = ""
        viz_prefix = ""

    # Locate attributes container for nodes
    node_attributes = None
    attrs_query = './/g:attributes' if ns else './/attributes'
    for attrs in root.findall(attrs_query, ns):
        if attrs.attrib.get('class') == 'node':
            node_attributes = attrs
            break

    if node_attributes is None:
        graph = root.find('g:graph', ns) if ns else root.find('graph')
        node_attributes = ET.SubElement(graph, f"{tag_prefix}attributes", {'class': 'node', 'mode': 'static'})

    # Find highest attribute id
    max_id = 0
    attr_query = 'g:attribute' if ns else 'attribute'
    for attr in node_attributes.findall(attr_query, ns):
        try:
            max_id = max(max_id, int(attr.attrib.get('id', 0)))
        except ValueError:
            pass

    # Define new attribute IDs
    bgc_caz_overlap_attr_id = str(max_id + 1)
    overlap_attr_id = str(max_id + 2)
    is_bgc_attr_id = str(max_id + 3)
    bgc_type_attr_id = str(max_id + 4)
    is_cazyme_attr_id = str(max_id + 5)
    caz_fam_attr_id = str(max_id + 6)
    caz_cls_attr_id = str(max_id + 7)

    # Define new attribute elements ordered right after 'product'
    new_attrs = [
        ET.Element(f"{tag_prefix}attribute", {'id': bgc_caz_overlap_attr_id, 'title': 'bgc_cazyme_overlap', 'type': 'string'}),
        ET.Element(f"{tag_prefix}attribute", {'id': bgc_type_attr_id, 'title': 'bgc_type', 'type': 'string'}),
        ET.Element(f"{tag_prefix}attribute", {'id': caz_fam_attr_id, 'title': 'cazyme_family', 'type': 'string'}),
        ET.Element(f"{tag_prefix}attribute", {'id': caz_cls_attr_id, 'title': 'cazyme_class', 'type': 'string'}),
        ET.Element(f"{tag_prefix}attribute", {'id': is_bgc_attr_id, 'title': 'is_bgc', 'type': 'string'}),
        ET.Element(f"{tag_prefix}attribute", {'id': is_cazyme_attr_id, 'title': 'is_cazyme', 'type': 'string'}),
        ET.Element(f"{tag_prefix}attribute", {'id': overlap_attr_id, 'title': 'overlap_status', 'type': 'string'})
    ]

    prod_idx = -1
    for i, child in enumerate(list(node_attributes)):
        if child.attrib.get('title') == 'product':
            prod_idx = i
            break
    if prod_idx != -1:
        for j, new_attr in enumerate(new_attrs):
            node_attributes.insert(prod_idx + 1 + j, new_attr)
    else:
        for new_attr in new_attrs:
            node_attributes.append(new_attr)

    # Color Palette definitions (RGBA)
    # BGC_and_CAZyme: Vivid Purple / Magenta (#9B59B6)
    # BGC_only: Coral / Vivid Orange (#E67E22)
    # CAZyme_only: Sky Blue / Cyan (#3498DB)
    # Other: Soft Muted Gray (#D5DBDB)
    COLOR_MAP = {
        'BGC_and_CAZyme': {'r': '155', 'g': '89', 'b': '182', 'a': '1.0'},
        'BGC_only':       {'r': '230', 'g': '126', 'b': '34', 'a': '1.0'},
        'CAZyme_only':    {'r': '52',  'g': '152', 'b': '219', 'a': '1.0'},
        'Other':          {'r': '213', 'g': '219', 'b': '219', 'a': '1.0'}
    }

    nodes_query = './/g:node' if ns else './/node'
    attvalues_query = 'g:attvalues' if ns else 'attvalues'
    attvalue_query = 'g:attvalue' if ns else 'attvalue'
    viz_color_query = 'viz:color' if ns else '{*}color'

    nodes = root.findall(nodes_query, ns)
    status_counts = {'BGC_and_CAZyme': 0, 'BGC_only': 0, 'CAZyme_only': 0, 'Other': 0}

    for node in nodes:
        node_id = str(node.attrib.get('id', '')).strip()
        node_label = str(node.attrib.get('label', '')).strip()

        attvalues = node.find(attvalues_query, ns)
        if attvalues is None:
            attvalues = ET.SubElement(node, f"{tag_prefix}attvalues")

        # Resolve cluster metadata
        meta = None
        if node_label in family_meta:
            meta = family_meta[node_label]
        elif node_id in family_meta:
            meta = family_meta[node_id]
        
        # Fallback inspection across existing node attribute values
        if meta is None:
            found_caz = set()
            found_cls = set()
            found_bgc = set()
            for attval in attvalues.findall(attvalue_query, ns):
                val = str(attval.attrib.get('value', ''))
                for gid in val.replace(";", " ").replace(",", " ").split():
                    t = gid.strip()
                    if t in cazyme_map:
                        found_caz.update(cazyme_map[t]['families'])
                        found_cls.update(cazyme_map[t]['classes'])
                    if t in bgc_map:
                        found_bgc.add(bgc_map[t]['bgc_type'])
            
            is_c = len(found_caz) > 0
            is_b = len(found_bgc) > 0
            if is_b and is_c:
                st = "BGC_and_CAZyme"
            elif is_b:
                st = "BGC_only"
            elif is_c:
                st = "CAZyme_only"
            else:
                st = "Other"
            meta = {
                'status': st,
                'is_bgc': "True" if is_b else "False",
                'bgc_type': ";".join(sorted(found_bgc)) if found_bgc else "None",
                'is_cazyme': "True" if is_c else "False",
                'cazyme_family': ";".join(sorted(found_caz)) if found_caz else "None",
                'cazyme_class': ";".join(sorted(found_cls)) if found_cls else "Non-CAZyme"
            }

        cur_status = meta['status']
        status_counts[cur_status] = status_counts.get(cur_status, 0) + 1

        # Apply Gephi <viz:color>
        col = COLOR_MAP[cur_status]
        v_color = node.find(viz_color_query, ns)
        if v_color is not None:
            v_color.attrib['r'] = col['r']
            v_color.attrib['g'] = col['g']
            v_color.attrib['b'] = col['b']
            v_color.attrib['a'] = col['a']
        else:
            ET.SubElement(node, f"{viz_prefix}color", {'r': col['r'], 'g': col['g'], 'b': col['b'], 'a': col['a']})

        # Append attribute values
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': bgc_caz_overlap_attr_id, 'value': cur_status})
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': overlap_attr_id, 'value': cur_status})
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': is_bgc_attr_id, 'value': meta['is_bgc']})
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': bgc_type_attr_id, 'value': meta['bgc_type']})
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': is_cazyme_attr_id, 'value': meta['is_cazyme']})
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': caz_fam_attr_id, 'value': meta['cazyme_family']})
        ET.SubElement(attvalues, f"{tag_prefix}attvalue", {'for': caz_cls_attr_id, 'value': meta['cazyme_class']})

    tree.write(str(gexf_out), encoding='utf-8', xml_declaration=True)
    print(f"Enriched {gexf_out.name}: {len(nodes)} nodes processed.")
    print(f"Status breakdown: {status_counts}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python inject_gexf_metadata.py <gexf_in> <pangenome_csv> <cazyme_csv> <gexf_out> [bgc_dir] [overlap_csv_out]")
        sys.exit(1)
    bgc_dir_arg = sys.argv[5] if len(sys.argv) > 5 else None
    overlap_csv_arg = sys.argv[6] if len(sys.argv) > 6 else None
    inject_cazyme_and_bgc_to_gexf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], bgc_dir_arg, overlap_csv_arg)
