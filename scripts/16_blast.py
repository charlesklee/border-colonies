#!/usr/bin/env python3
"""
16_blast.py — NCBI core_nt BLAST REST Query & Entrez Source/Host Lookup (Table S2)
==================================================================================
Purpose : Submit the 39 canonical manuscript target ASVs to NCBI core_nt
          (megablast, HITLIST_SIZE=50, keep top 5 by bit score) via the NCBI
          BLAST REST API and retrieve organism / isolation_source / host from
          NCBI Entrez E-utilities. Combined with scripts/16_blast_local16s.R
          (local Biostrings alignment against RefSeq 16S type strains & COTS27)
          to produce Table S2 (31_ASV_FULL_BLAST_TABLE_v6.csv & blast_summary_v6).
Inputs  : output/intermediate/block19_target_ASVs.csv
Outputs : output/intermediate/block19_blast_cache.json
          output/tables/TableS2_full_blast.csv
          output/tables/TableS2_blast_summary.csv
Supports: Table S2; §3.5 / §4 (ASV_8, ASV_17, COTS27, Yarza et al. 2014 flags)
Note    : Replace NCBI_EMAIL below before running against NCBI APIs.
          Reported manuscript results were queried on 2026-11-23.
"""

import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

NCBI_EMAIL  = os.environ.get("NCBI_EMAIL", "your.email@institution.edu")
TARGETS_CSV = os.path.join("output", "intermediate", "block19_target_ASVs.csv")
CACHE_JSON  = os.path.join("output", "intermediate", "block19_blast_cache.json")
OUT_LONG    = os.path.join("output", "tables", "TableS2_full_blast.csv")
OUT_SUMMARY = os.path.join("output", "tables", "TableS2_blast_summary.csv")

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS   = {"User-Agent": f"AscidianMicrobiomePipeline/1.0 (mailto:{NCBI_EMAIL})"}


def http_post(url, data_dict, timeout=120):
    encoded = urllib.parse.urlencode(data_dict).encode("utf-8")
    req = urllib.request.Request(url, data=encoded, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def http_get(url, timeout=120):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def submit_blast(fasta_str, database="core_nt", megablast=True, hitlist_size=50, expect="1e-10", entrez_query=None):
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": database,
        "QUERY": fasta_str,
        "HITLIST_SIZE": str(hitlist_size),
        "EXPECT": str(expect),
        "FORMAT_TYPE": "XML",
        "TOOL": "AscidianMicrobiomePipeline",
    }
    if entrez_query:
        params["ENTREZ_QUERY"] = entrez_query
    if megablast:
        params["MEGABLAST"] = "on"
    else:
        params["WORD_SIZE"] = "11"
    text = http_post(BLAST_URL, params)
    m_rid = re.search(r"RID = ([A-Z0-9]+)", text)
    m_rtoe = re.search(r"RTOE = (\d+)", text)
    if not m_rid:
        raise RuntimeError(f"Failed to obtain RID:\n{text[:600]}")
    return m_rid.group(1), (int(m_rtoe.group(1)) if m_rtoe else 15)


def poll_and_fetch(rid, rtoe, max_wait=300):
    wait0 = max(10, min(rtoe, 25))
    time.sleep(wait0)
    elapsed = wait0
    while elapsed < max_wait:
        st = http_get(f"{BLAST_URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}")
        if "Status=READY" in st:
            if "ThereAreHits=no" in st:
                return None
            return http_get(f"{BLAST_URL}?CMD=Get&FORMAT_TYPE=XML&RID={rid}", timeout=180)
        if "Status=FAILED" in st or "Status=UNKNOWN" in st:
            raise RuntimeError(f"BLAST RID {rid} failed/expired")
        time.sleep(15)
        elapsed += 15
    raise TimeoutError(f"BLAST RID {rid} timed out after {max_wait}s")


def parse_xml(xml_txt, top_n=5):
    res = {}
    if not xml_txt:
        return res
    root = ET.fromstring(xml_txt)
    for it in root.findall(".//Iteration"):
        qid = it.findtext("Iteration_query-def", "").strip().split()[0]
        qlen = int(it.findtext("Iteration_query-len", "0"))
        hits = []
        for hit in it.findall(".//Hit"):
            hit_id = hit.findtext("Hit_id", "").strip()
            acc = hit.findtext("Hit_accession", "").strip()
            m_ver = re.search(r"\|([A-Z]{1,4}_?[0-9]+\.[0-9]+)\|", hit_id)
            acc_ver = m_ver.group(1) if m_ver else acc
            hdef = hit.findtext("Hit_def", "").strip()
            hsps = hit.findall(".//Hsp")
            if not hsps:
                continue
            bh = max(hsps, key=lambda h: float(h.findtext("Hsp_bit-score", "0")))
            al = int(bh.findtext("Hsp_align-len", "1"))
            ident = int(bh.findtext("Hsp_identity", "0"))
            qf, qt = int(bh.findtext("Hsp_query-from", "1")), int(bh.findtext("Hsp_query-to", "1"))
            hits.append({
                "Accession": acc_ver,
                "Description": hdef,
                "Pct_identity": round(100.0 * ident / al, 2),
                "Query_cover": round(100.0 * (abs(qt - qf) + 1) / qlen, 1) if qlen else 0.0,
                "Align_len": al,
                "E_value": float(bh.findtext("Hsp_evalue", "1")),
                "Bit_score": round(float(bh.findtext("Hsp_bit-score", "0")), 1),
            })
        hits.sort(key=lambda x: (-x["Bit_score"], x["E_value"]))
        res[qid] = hits[:top_n]
    return res


def fetch_nuccore(accessions):
    meta = {}
    accs = sorted(set(a for a in accessions if a))
    for i in range(0, len(accs), 20):
        ids_str = ",".join(accs[i:i + 20])
        url = f"{EUTILS}/efetch.fcgi?db=nuccore&id={urllib.parse.quote(ids_str)}&rettype=gb&retmode=xml"
        try:
            root = ET.fromstring(http_get(url, timeout=120))
            for seq in root.findall(".//GBSeq"):
                acc_v = seq.findtext("GBSeq_accession-version", "").strip()
                org = seq.findtext("GBSeq_organism", "").strip()
                iso, host, country = "", "", ""
                for feat in seq.findall(".//GBFeature"):
                    if feat.findtext("GBFeature_key", "") == "source":
                        for q in feat.findall(".//GBQualifier"):
                            qn, qv = q.findtext("GBQualifier_name", ""), q.findtext("GBQualifier_value", "").strip()
                            if qn == "isolation_source" and not iso:
                                iso = qv
                            elif qn == "host" and not host:
                                host = qv
                            elif qn in ("geo_loc_name", "country") and not country:
                                country = qv
                parts = ([f"host: {host}"] if host else []) + ([f"source: {iso}"] if iso else []) + ([country] if country else [])
                meta[acc_v] = {"Scientific_name": org, "Isolation_source_or_host": "; ".join(parts) if parts else "Not specified"}
        except Exception as e:
            print(f"Warning: Entrez fetch failed: {e}", flush=True)
        time.sleep(0.5)
    return meta


def main():
    os.makedirs(os.path.dirname(OUT_LONG), exist_ok=True)
    with open(TARGETS_CSV, newline="", encoding="utf-8") as f:
        targets = list(csv.DictReader(f))
    cache = json.load(open(CACHE_JSON, encoding="utf-8")) if os.path.exists(CACHE_JSON) else {"core_nt": {}, "16S_type": {}, "nuccore_meta": {}, "cots27": {}}

    missing = [t for t in targets if t["ASV"] not in cache["core_nt"]]
    for i in range(0, len(missing), 3):
        batch = missing[i:i + 3]
        fasta = "\n".join(f">{t['ASV']}\n{t['Sequence']}" for t in batch)
        try:
            rid, rtoe = submit_blast(fasta, "core_nt", True, 50)
            parsed = parse_xml(poll_and_fetch(rid, rtoe), 5)
            for t in batch:
                cache["core_nt"][t["ASV"]] = parsed.get(t["ASV"], [])
        except Exception:
            for t in batch:
                f1 = f">{t['ASV']}\n{t['Sequence']}"
                rid1, rtoe1 = submit_blast(f1, "core_nt", True, 50, entrez_query="16S[Title] OR rRNA[Title] OR ribosomal[Title]")
                cache["core_nt"][t["ASV"]] = parse_xml(poll_and_fetch(rid1, rtoe1), 5).get(t["ASV"], [])
        json.dump(cache, open(CACHE_JSON, "w", encoding="utf-8"), indent=2)
        time.sleep(3)

    accs = {h["Accession"] for db in ("core_nt", "16S_type") for hl in cache[db].values() for h in hl}
    miss_a = [a for a in accs if a not in cache["nuccore_meta"]]
    if miss_a:
        cache["nuccore_meta"].update(fetch_nuccore(miss_a))
        json.dump(cache, open(CACHE_JSON, "w", encoding="utf-8"), indent=2)

    long_rows, sum_rows = [], []
    for t in targets:
        aid = t["ASV"]
        for db in ("core_nt", "16S_type"):
            for rk, h in enumerate(cache[db].get(aid, []), 1):
                m = cache["nuccore_meta"].get(h["Accession"], {})
                long_rows.append({
                    "ASV": aid, "Old_IDs": t["Old_IDs_from_31table"] or t["Old_IDs_crosswalk"],
                    "Length_bp": t["Length_bp"], "Total_reads": t["Total_reads"],
                    "Dominant_species_tissue": t["Dominant_species_tissue"],
                    "SILVA_138.2_taxonomy": t["SILVA_138.2_taxonomy"],
                    "Database": db, "Hit_rank": rk, "Accession": h["Accession"],
                    "Description": h["Description"], "Scientific_name": m.get("Scientific_name", ""),
                    "Pct_identity": h["Pct_identity"], "Query_cover": h["Query_cover"],
                    "Align_len": h["Align_len"], "E_value": f"{h['E_value']:.2e}",
                    "Bit_score": h["Bit_score"], "Isolation_source_or_host": m.get("Isolation_source_or_host", ""),
                })
        bc = cache["core_nt"].get(aid, [None])[0]
        bt = cache["16S_type"].get(aid, [None])[0]
        ts_id = bt["Pct_identity"] if bt else 0.0
        note = f"Potential novel family (<86.5% type-strain id; {ts_id:.2f}%)" if ts_id < 86.5 else (f"Potential novel genus (<94.5% type-strain id; {ts_id:.2f}%)" if ts_id < 94.5 else f"Known genus-level match (>=94.5%; {ts_id:.2f}%)")
        sum_rows.append({
            "ASV": aid, "SILVA_138.2_taxonomy": t["SILVA_138.2_taxonomy"],
            "Best_core_nt_accession": bc["Accession"] if bc else "None",
            "Best_core_nt_pct_id": bc["Pct_identity"] if bc else "",
            "Best_core_nt_qcov": bc["Query_cover"] if bc else "",
            "Best_core_nt_source_host": cache["nuccore_meta"].get(bc["Accession"], {}).get("Isolation_source_or_host", "") if bc else "",
            "Best_type_strain_accession": bt["Accession"] if bt else "None",
            "Best_type_strain_name": cache["nuccore_meta"].get(bt["Accession"], {}).get("Scientific_name", "") if bt else "",
            "Best_type_strain_pct_id": bt["Pct_identity"] if bt else "",
            "Note_Yarza2014": note,
        })

    if long_rows:
        with open(OUT_LONG, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(long_rows[0].keys())); w.writeheader(); w.writerows(long_rows)
    if sum_rows:
        with open(OUT_SUMMARY, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(sum_rows[0].keys())); w.writeheader(); w.writerows(sum_rows)
    print(f"Wrote {len(long_rows)} rows to {OUT_LONG} and {len(sum_rows)} rows to {OUT_SUMMARY}.")


if __name__ == "__main__":
    main()
