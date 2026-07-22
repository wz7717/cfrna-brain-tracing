#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "braintrace_source_tracing.db"
DEFAULT_OUTDIR = ROOT / "reports"
DEFAULT_ORTHOLOGY = ROOT / "data" / "orthology" / "ensembl_mfascicularis_hsapiens_homology.tsv"


CURATION_ROWS = [
    {
        "old_region_id": "44563",
        "new_region_id": "1-2",
        "region_name": "somatosensory areas 1 and 2",
        "region_acronym": "1-2",
        "parent_region_id": "Parietal",
        "lobe": "Parietal",
        "roi173": "1-2",
        "regional_map": "S1",
        "saleem_network": "Parietal, and Parieto-occipital region",
        "status": "corrected_id",
        "evidence": "Bo2023 joined annotation has Region=1-2, Full_name=somatosensory areas 1 and 2; current 44563 came from the 1-2 samples and is consistent with an Excel serial/date conversion artifact.",
    },
    {
        "old_region_id": "AI",
        "new_region_id": "A1",
        "region_name": "auditory area I, core region of the auditory cortex",
        "region_acronym": "A1",
        "parent_region_id": "Temporal",
        "lobe": "Temporal",
        "roi173": "A1",
        "regional_map": "A1",
        "saleem_network": "Temporal",
        "status": "corrected_id",
        "evidence": "Original sample rows have roi173=A1 and rm94=A1; Bo2023 joined annotation uses Region=A1 with this full name.",
    },
    {
        "old_region_id": "PrCo",
        "new_region_id": "PrCO",
        "region_name": "precentral opercular area",
        "region_acronym": "PrCO",
        "parent_region_id": "Frontal",
        "lobe": "Frontal",
        "roi173": "PrCO",
        "regional_map": "PFCol",
        "saleem_network": "Operculum/Insula",
        "status": "case_corrected_id",
        "evidence": "Bo2023 abbreviation/joined annotation uses PrCO and full name precentral opercular area.",
    },
    {
        "old_region_id": "NARegion",
        "new_region_id": "NAcc",
        "region_name": "nucleus accumbens",
        "region_acronym": "NAcc",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "Striatum",
        "regional_map": "NAcc",
        "saleem_network": "Subcortical",
        "status": "corrected_id",
        "evidence": "Samples labelled NARegion have rm94=NAcc; Bo2023 dictionary maps NAcc to nucleus accumbens.",
    },
    {
        "old_region_id": "GPeGPi",
        "new_region_id": "GP",
        "region_name": "globus pallidus",
        "region_acronym": "GPeGPi",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "GPeGPi",
        "regional_map": "GP",
        "saleem_network": "Subcortical",
        "status": "corrected_id",
        "evidence": "Bo2023 joined annotation uses Region=GP, roi173=GPeGPi, Full_name=globus pallidus.",
    },
    {
        "old_region_id": "amy",
        "new_region_id": "Amy",
        "region_name": "amygdala",
        "region_acronym": "Amy",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "amy",
        "regional_map": "Amyg",
        "saleem_network": "Subcortical",
        "status": "case_corrected_id",
        "evidence": "Bo2023 dictionary uses Amy for amygdala.",
    },
    {
        "old_region_id": "cd",
        "new_region_id": "Cd",
        "region_name": "caudate",
        "region_acronym": "Cd",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "Striatum",
        "regional_map": "Cau",
        "saleem_network": "Subcortical",
        "status": "case_corrected_id",
        "evidence": "Bo2023 dictionary uses Cd for caudate.",
    },
    {
        "old_region_id": "thalamus",
        "new_region_id": "Tha",
        "region_name": "thalamus",
        "region_acronym": "Tha",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "TH",
        "regional_map": "TH",
        "saleem_network": "Subcortical",
        "status": "corrected_id",
        "evidence": "Bo2023 dictionary uses Tha for thalamus.",
    },
    {
        "old_region_id": "pu",
        "new_region_id": "pu",
        "region_name": "putamen",
        "region_acronym": "Pu",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "Striatum",
        "regional_map": "Put",
        "saleem_network": "Subcortical",
        "status": "display_name_completed",
        "evidence": "Bo2023 dictionary lists Pu as putamen; joined annotation keeps Region=pu with roi173=Striatum.",
    },
    {
        "old_region_id": "PL",
        "new_region_id": "PL",
        "region_name": "PL subcortical region",
        "region_acronym": "PL",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "PL",
        "regional_map": "PL",
        "saleem_network": "Subcortical",
        "status": "display_name_completed_unresolved",
        "evidence": "Bo2023 joined annotation contains Region=PL, roi173=PL, but the local dictionary lacks a full anatomical expansion; keep ID and flag as unresolved for publication wording.",
    },
    {
        "old_region_id": "MT",
        "new_region_id": "MT",
        "region_name": "middle temporal visual area",
        "region_acronym": "MT",
        "parent_region_id": "Temporal",
        "lobe": "Temporal",
        "roi173": "MT",
        "regional_map": "MT",
        "saleem_network": "Occipital/Temporal",
        "status": "display_name_completed",
        "evidence": "Bo2023 joined annotation contains Region=MT, roi173=MT; standard macaque abbreviation MT is middle temporal visual area.",
    },
    {
        "old_region_id": "23a",
        "new_region_id": "23a",
        "region_name": "posterior cingulate cortex, area 23a",
        "region_acronym": "23a",
        "parent_region_id": "Cingulate",
        "lobe": "Cingulate",
        "roi173": "23a",
        "regional_map": "23a",
        "saleem_network": "Cingulate gyrus",
        "status": "display_name_completed",
        "evidence": "Bo2023 joined annotation contains Region=23a; companion regions 23b/23c are labelled posterior cingulate cortex in the current atlas.",
    },
    {
        "old_region_id": "31",
        "new_region_id": "31",
        "region_name": "posterior cingulate cortex, area 31",
        "region_acronym": "31",
        "parent_region_id": "Cingulate",
        "lobe": "Cingulate",
        "roi173": "31",
        "regional_map": "31",
        "saleem_network": "Cingulate gyrus",
        "status": "display_name_completed",
        "evidence": "Bo2023 joined annotation contains Region=31; assigned cingulate display name consistent with parent/lobe labels.",
    },
    {
        "old_region_id": "cla",
        "new_region_id": "Cla",
        "region_name": "claustrum",
        "region_acronym": "Cla",
        "parent_region_id": "Subcortical",
        "lobe": "Subcortical",
        "roi173": "Claustrum",
        "regional_map": "Claustrum",
        "saleem_network": "Subcortical",
        "status": "case_corrected_id",
        "evidence": "Samples labelled cla have roi173=Claustrum; Bo2023 joined annotation/dictionary identify the region as claustrum.",
    },
]


def merge_coordinates(existing: str | None, row: dict[str, str]) -> str:
    try:
        data = json.loads(existing or "{}")
    except json.JSONDecodeError:
        data = {}
    data.update(
        {
            "lobe": row["lobe"],
            "saleem_network": row["saleem_network"],
            "roi173": row["roi173"],
            "regional_map": row["regional_map"],
            "label_curation_status": row["status"],
            "label_curation_evidence": row["evidence"],
            "curated_region_id": row["new_region_id"],
        }
    )
    return json.dumps(data, ensure_ascii=False)


def curate_region_labels(db_path: Path, outdir: Path, dry_run: bool) -> pd.DataFrame:
    curation = pd.DataFrame(CURATION_ROWS)
    outdir.mkdir(parents=True, exist_ok=True)
    curation_path = outdir / "bo2023_publication_label_curation_map_20260704.csv"
    curation.to_csv(curation_path, index=False, encoding="utf-8-sig")
    if dry_run:
        return curation

    backup = db_path.with_suffix(f".pre_bo2023_label_curation_20260704{db_path.suffix}")
    if not backup.exists():
        shutil.copy2(db_path, backup)

    con = sqlite3.connect(db_path)
    try:
        cur = con.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS bo2023_publication_label_curation (
                old_region_id TEXT PRIMARY KEY,
                new_region_id TEXT NOT NULL,
                region_name TEXT NOT NULL,
                region_acronym TEXT,
                parent_region_id TEXT,
                lobe TEXT,
                roi173 TEXT,
                regional_map TEXT,
                saleem_network TEXT,
                status TEXT,
                evidence TEXT,
                updated_at TEXT
            )
            """
        )
        updated_at = datetime.now(timezone.utc).isoformat()
        for row in CURATION_ROWS:
            old_id = row["old_region_id"]
            new_id = row["new_region_id"]
            existing = cur.execute(
                "SELECT coordinates FROM macaque_brain_atlas WHERE atlas_id=4 AND region_id=?",
                (old_id,),
            ).fetchone()
            coordinates = merge_coordinates(existing[0] if existing else None, row)
            if old_id != new_id:
                old_exists = cur.execute(
                    "SELECT COUNT(*) FROM macaque_brain_atlas WHERE atlas_id=4 AND region_id=?",
                    (old_id,),
                ).fetchone()[0]
                target = cur.execute(
                    "SELECT COUNT(*) FROM macaque_brain_atlas WHERE atlas_id=4 AND region_id=?",
                    (new_id,),
                ).fetchone()[0]
                if target and old_exists:
                    raise ValueError(f"Cannot rename {old_id} -> {new_id}: target already exists")
                if old_exists:
                    cur.execute(
                        """
                        UPDATE signature_genes
                        SET region_id=?
                        WHERE sigset_id IN (SELECT sigset_id FROM signature_sets WHERE atlas_id=4)
                          AND region_id=?
                        """,
                        (new_id, old_id),
                    )
                    cur.execute(
                        """
                        UPDATE reference_expression
                        SET region_id=?, region_name=?
                        WHERE atlas_id=4 AND region_id=?
                        """,
                        (new_id, row["region_name"], old_id),
                    )
                    cur.execute(
                        """
                        UPDATE macaque_brain_atlas
                        SET region_id=?, region_name=?, region_acronym=?, parent_region_id=?, coordinates=?
                        WHERE atlas_id=4 AND region_id=?
                        """,
                        (
                            new_id,
                            row["region_name"],
                            row["region_acronym"],
                            row["parent_region_id"],
                            coordinates,
                            old_id,
                        ),
                    )
                else:
                    current = cur.execute(
                        "SELECT coordinates FROM macaque_brain_atlas WHERE atlas_id=4 AND region_id=?",
                        (new_id,),
                    ).fetchone()
                    refreshed_coordinates = merge_coordinates(current[0] if current else coordinates, row)
                    cur.execute(
                        """
                        UPDATE reference_expression
                        SET region_name=?
                        WHERE atlas_id=4 AND region_id=?
                        """,
                        (row["region_name"], new_id),
                    )
                    cur.execute(
                        """
                        UPDATE macaque_brain_atlas
                        SET region_name=?, region_acronym=?, parent_region_id=?, coordinates=?
                        WHERE atlas_id=4 AND region_id=?
                        """,
                        (
                            row["region_name"],
                            row["region_acronym"],
                            row["parent_region_id"],
                            refreshed_coordinates,
                            new_id,
                        ),
                    )
            else:
                cur.execute(
                    """
                    UPDATE reference_expression
                    SET region_name=?
                    WHERE atlas_id=4 AND region_id=?
                    """,
                    (row["region_name"], old_id),
                )
                cur.execute(
                    """
                    UPDATE macaque_brain_atlas
                    SET region_name=?, region_acronym=?, parent_region_id=?, coordinates=?
                    WHERE atlas_id=4 AND region_id=?
                    """,
                    (
                        row["region_name"],
                        row["region_acronym"],
                        row["parent_region_id"],
                        coordinates,
                        old_id,
                    ),
                )
            cur.execute(
                """
                INSERT OR REPLACE INTO bo2023_publication_label_curation
                (old_region_id, new_region_id, region_name, region_acronym, parent_region_id,
                 lobe, roi173, regional_map, saleem_network, status, evidence, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    old_id,
                    new_id,
                    row["region_name"],
                    row["region_acronym"],
                    row["parent_region_id"],
                    row["lobe"],
                    row["roi173"],
                    row["regional_map"],
                    row["saleem_network"],
                    row["status"],
                    row["evidence"],
                    updated_at,
                ),
            )

        notes_text = cur.execute("SELECT notes FROM atlas_versions WHERE atlas_id=4").fetchone()[0]
        notes = json.loads(notes_text)
        notes["publication_label_curation"] = {
            "updated_at_utc": updated_at,
            "map_file": str(curation_path),
            "policy": "Publication-facing Bo2023 region labels were curated from the local Bo2023 joined annotation/dictionary; expression values and signature weights were not recomputed.",
        }
        cur.execute(
            "UPDATE atlas_versions SET notes=? WHERE atlas_id=4",
            (json.dumps(notes, ensure_ascii=False),),
        )
        con.commit()
    finally:
        con.close()
    return curation


def build_humanized_signature(db_path: Path, orthology_path: Path, outdir: Path) -> pd.DataFrame:
    sig = pd.read_sql_query(
        """
        SELECT sg.sigset_id, sg.region_id, mba.region_name, mba.parent_region_id,
               sg.gene_symbol AS macaque_signature_gene, sg.weight
        FROM signature_genes sg
        JOIN signature_sets ss ON ss.sigset_id = sg.sigset_id
        LEFT JOIN macaque_brain_atlas mba ON mba.atlas_id = ss.atlas_id AND mba.region_id = sg.region_id
        WHERE ss.atlas_id=4
        """,
        sqlite3.connect(db_path),
    )
    ortho = pd.read_csv(orthology_path, sep="\t")
    ortho = ortho.rename(
        columns={
            "Gene stable ID": "macaque_gene_id",
            "Gene name": "macaque_gene_name",
            "Human gene stable ID": "human_gene_id",
            "Human gene name": "human_gene_name",
            "Human homology type": "human_homology_type",
            "Human orthology confidence [0 low, 1 high]": "human_orthology_confidence",
        }
    )
    for col in ["macaque_gene_id", "macaque_gene_name", "human_gene_name"]:
        ortho[col] = ortho[col].fillna("").astype(str).str.strip()
    one2one = ortho[
        (ortho["human_homology_type"].fillna("").astype(str).str.contains("ortholog", case=False))
        & (pd.to_numeric(ortho["human_orthology_confidence"], errors="coerce").fillna(0) >= 1)
        & (ortho["human_gene_name"] != "")
    ].copy()
    id_map = dict(zip(one2one["macaque_gene_id"], one2one["human_gene_name"]))
    name_map = dict(zip(one2one["macaque_gene_name"], one2one["human_gene_name"]))

    def humanize(gene: str) -> tuple[str, str]:
        gene = str(gene).strip()
        if gene.startswith("ENSMFAG"):
            return id_map.get(gene, ""), "orthology_id_map" if gene in id_map else "unmapped_ensmfag"
        return name_map.get(gene, gene), "orthology_name_map" if gene in name_map else "as_symbol"

    mapped = sig["macaque_signature_gene"].map(humanize)
    sig["human_gene_symbol"] = [x[0] for x in mapped]
    sig["humanization_status"] = [x[1] for x in mapped]
    out = outdir / "bo2023_signature_genes_humanized_20260704.csv"
    sig.to_csv(out, index=False, encoding="utf-8-sig")
    return sig


def write_report(outdir: Path, curation: pd.DataFrame, humanized: pd.DataFrame, db_path: Path) -> None:
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    signature_regions = cur.execute(
        "SELECT COUNT(DISTINCT region_id) FROM signature_genes WHERE sigset_id IN (SELECT sigset_id FROM signature_sets WHERE atlas_id=4)"
    ).fetchone()[0]
    missing_labels = cur.execute(
        """
        SELECT COUNT(*) FROM signature_genes sg
        WHERE sg.sigset_id IN (SELECT sigset_id FROM signature_sets WHERE atlas_id=4)
          AND NOT EXISTS (
            SELECT 1 FROM macaque_brain_atlas mba
            WHERE mba.atlas_id=4 AND mba.region_id=sg.region_id
          )
        """
    ).fetchone()[0]
    con.close()

    status_counts = humanized["humanization_status"].value_counts().to_dict()
    mapped_rows = int((humanized["human_gene_symbol"].fillna("") != "").sum())
    lines = [
        "# Bo2023 publication label curation report",
        "",
        f"Database: `{db_path}`",
        "",
        "## Region-label curation",
        "",
        f"- Curated rows: {len(curation)}",
        f"- Signature regions after curation: {signature_regions}",
        f"- Signature rows with missing atlas label after curation: {missing_labels}",
        "",
        "Corrected IDs:",
    ]
    for row in curation.itertuples(index=False):
        if row.old_region_id != row.new_region_id:
            lines.append(f"- `{row.old_region_id}` -> `{row.new_region_id}`: {row.region_name} ({row.status})")
    lines.extend(["", "Display-name completions:"])
    for row in curation.itertuples(index=False):
        if row.old_region_id == row.new_region_id:
            lines.append(f"- `{row.old_region_id}`: {row.region_name} ({row.status})")
    lines.extend(
        [
            "",
            "## Humanized signature mapping",
            "",
            f"- Humanized signature rows with a non-empty human gene symbol: {mapped_rows}/{len(humanized)}",
            f"- Humanization status counts: {status_counts}",
            "",
            "The original macaque signature genes are preserved. The humanized CSV is intended for human cfRNA overlap/enrichment analyses.",
        ]
    )
    (outdir / "bo2023_publication_label_curation_report_20260704.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--orthology", type=Path, default=DEFAULT_ORTHOLOGY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    curation = curate_region_labels(args.db, args.outdir, args.dry_run)
    humanized = build_humanized_signature(args.db, args.orthology, args.outdir)
    write_report(args.outdir, curation, humanized, args.db)
    print(args.outdir / "bo2023_publication_label_curation_report_20260704.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
