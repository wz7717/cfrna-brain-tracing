#!/usr/bin/env python3
"""Round-4 P2 calculations: ablation bounds, matched RF, F1/IQR, enrichment, confusion."""
from __future__ import annotations

import argparse, json, math, csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "raw_datasets_v0.1.9_20260728/01_Bo2023"
OUT = ROOT / "manuscript/calculations/p2"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260730
ROUTE = "hybrid_projected_network_logcpm_exact"

parser = argparse.ArgumentParser()
parser.add_argument(
    "--rf-only",
    action="store_true",
    help="Recompute only the frozen-truth 200-gene RF outputs and update P2_audit.json.",
)
args = parser.parse_args()

panel = pd.read_csv(ROOT / "code/data/models/bo2023_saleem_network_top200_model_genes.csv")
panel_idx = set(panel["gene_index"].astype(int))
vsd_path = RAW / "primary/mfas5_819samples_23605genes_vsd4_rmbatch.xls"
with vsd_path.open(encoding="utf-8") as f:
    sample_ids = f.readline().rstrip("\r\n").split("\t")
    selected = {}
    for idx, line in enumerate(f):
        if idx in panel_idx:
            parts = line.rstrip("\r\n").split("\t")
            selected[idx] = np.asarray(parts[1:], dtype=np.float32)
X = np.vstack([selected[i] for i in panel["gene_index"].astype(int)]).T

frozen_truth_path = (
    ROOT
    / "code/reproducibility/p2_rf200_frozen_truth/formal_lomo_network_detail.csv"
)
frozen = pd.read_csv(
    frozen_truth_path,
    dtype={"sample_id": str, "monkey_id": str},
)
frozen = frozen.loc[frozen["route_family"].eq(ROUTE)].copy()
frozen["sample_id"] = frozen["sample_id"].str.strip()
frozen["monkey_id"] = frozen["monkey_id"].str.strip()
if len(frozen) != 819 or frozen["sample_id"].nunique() != 819:
    raise ValueError("Frozen Network truth is not the expected 819-sample route")
if frozen["monkey_id"].nunique() != 9:
    raise ValueError("Frozen Network truth does not contain nine donor folds")
frozen = frozen.drop_duplicates("sample_id").set_index("sample_id")
if set(frozen.index) != set(sample_ids):
    raise ValueError("Frozen Network truth and VSD sample sets differ")
frozen = frozen.loc[sample_ids]
y = frozen["label"].astype(str).to_numpy()
donor = frozen["monkey_id"].astype(str).to_numpy()
classes = sorted(set(y))

details = []
matrices = {}
for heldout in sorted(set(donor)):
    train, test = donor != heldout, donor == heldout
    model = RandomForestClassifier(
        n_estimators=300, random_state=SEED, class_weight="balanced_subsample",
        n_jobs=-1, min_samples_leaf=2,
    )
    model.fit(X[train], y[train])
    prob = model.predict_proba(X[test])
    order = np.argsort(prob, axis=1)[:, ::-1]
    pred = model.classes_[order[:, 0]]
    top3 = model.classes_[order[:, :3]]
    truths = y[test]
    ids = np.asarray(sample_ids)[test]
    for sid, truth, p1, p3 in zip(ids, truths, pred, top3):
        details.append({
            "sample_id": sid, "donor": heldout, "truth": truth, "pred_top1": p1,
            "pred_top3": " | ".join(p3), "hit1": int(p1 == truth),
            "hit3": int(truth in p3), "panel_genes": 200,
        })
    matrices[heldout] = confusion_matrix(truths, pred, labels=classes)
detail = pd.DataFrame(details)
detail.to_csv(OUT / "P2_RF200_lomo_detail.csv", index=False)
rf_summary = {
    "n": len(detail), "donors": detail["donor"].nunique(),
    "top1_hits": int(detail.hit1.sum()), "top1": float(detail.hit1.mean()),
    "top3_hits": int(detail.hit3.sum()), "top3": float(detail.hit3.mean()),
    "specification": "locked 200-gene panel; 300 trees; balanced_subsample; min_samples_leaf=2; frozen-truth nine-fold LOMO",
    "truth_source": str(frozen_truth_path.relative_to(ROOT)).replace("\\", "/"),
    "truth_route_family": ROUTE,
    "random_state": SEED,
}

# donor matrices: long-form CSV and a dependency-light SVG small-multiple figure
long = []
for d, cm in matrices.items():
    for i, truth in enumerate(classes):
        for j, pred in enumerate(classes):
            long.append({"donor": d, "truth": truth, "predicted": pred, "count": int(cm[i, j])})
pd.DataFrame(long).to_csv(OUT / "P2_RF200_donor_confusion_long.csv", index=False)

def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;")
cell, gap, left, top = 12, 70, 85, 42
panel_w, panel_h = cell*10, cell*10
W, H = 3*(panel_w+gap)+left, 3*(panel_h+55)+top
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
       '<rect width="100%" height="100%" fill="white"/>',
       '<style>text{font-family:Arial;font-size:8px}.title{font-size:11px;font-weight:bold}</style>']
for k,(d,cm) in enumerate(sorted(matrices.items())):
    ox=left+(k%3)*(panel_w+gap); oy=top+(k//3)*(panel_h+55)
    mx=max(1,int(cm.max()))
    svg.append(f'<text class="title" x="{ox}" y="{oy-10}">Donor {esc(d)} (n={int(cm.sum())})</text>')
    for i in range(10):
        for j in range(10):
            v=int(cm[i,j]); shade=255-round(190*v/mx)
            svg.append(f'<rect x="{ox+j*cell}" y="{oy+i*cell}" width="{cell}" height="{cell}" fill="rgb({shade},{shade},{255})" stroke="#ddd" stroke-width=".3"/>')
            if v: svg.append(f'<text x="{ox+j*cell+cell/2}" y="{oy+i*cell+9}" text-anchor="middle">{v}</text>')
    for i,label in enumerate(classes):
        ab=esc(label[:4])
        svg.append(f'<text x="{ox-3}" y="{oy+i*cell+9}" text-anchor="end">{ab}</text>')
        svg.append(f'<text x="{ox+i*cell+6}" y="{oy+panel_h+10}" text-anchor="middle" transform="rotate(90 {ox+i*cell+6},{oy+panel_h+10})">{ab}</text>')
svg.append('</svg>')
(OUT/"P2_RF200_donor_confusion.svg").write_text("\n".join(svg),encoding="utf-8")
try:
    from PIL import Image, ImageDraw, ImageFont
    scale = 3
    im = Image.new("RGB", (W*scale, H*scale), "white")
    dr = ImageDraw.Draw(im)
    for k,(d,cm) in enumerate(sorted(matrices.items())):
        ox=(left+(k%3)*(panel_w+gap))*scale
        oy=(top+(k//3)*(panel_h+55))*scale
        mx=max(1,int(cm.max()))
        dr.text((ox,oy-30),f"Donor {d} (n={int(cm.sum())})",fill="black")
        for i in range(10):
            for j in range(10):
                v=int(cm[i,j]); shade=255-round(190*v/mx)
                box=(ox+j*cell*scale,oy+i*cell*scale,ox+(j+1)*cell*scale,oy+(i+1)*cell*scale)
                dr.rectangle(box,fill=(shade,shade,255),outline=(220,220,220))
                if v: dr.text((box[0]+5,box[1]+4),str(v),fill="black")
        for i,label in enumerate(classes):
            dr.text((ox-65,oy+i*cell*scale+5),label[:4],fill="black")
            dr.text((ox+i*cell*scale,oy+panel_h*scale+5),label[:4],fill="black")
    im.save(OUT/"P2_RF200_donor_confusion.png")
except Exception:
    pass

if args.rf_only:
    audit_path = OUT / "P2_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8")) if audit_path.exists() else {}
    audit["rf200"] = rf_summary
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(rf_summary, ensure_ascii=False, indent=2))
    raise SystemExit(0)

# Exact-region F1 median/IQR.
macro = json.loads((ROOT/"code/reproducibility/macro_f1_class_data.json").read_text(encoding="utf-8"))
f1_rows = pd.DataFrame(macro["data"])
f1_rows["f1"] = pd.to_numeric(f1_rows["f1"], errors="coerce")
f1_stats = {}
for endpoint in ["LOSO_Exact", "LOMO_Exact", "LOSO_Network", "LOMO_Network"]:
    x=f1_rows.loc[f1_rows.endpoint.eq(endpoint),"f1"].dropna().to_numpy()
    if len(x):
        q=np.quantile(x,[.25,.5,.75])
        f1_stats[endpoint]={"n_classes":len(x),"q1":float(q[0]),"median":float(q[1]),"q3":float(q[2]),"iqr":float(q[2]-q[0])}

# Bo2023 label to Saleem-style nomenclature crosswalk.
ann = pd.read_csv(RAW/"auxiliary_buildkit/33_region_annotation_joined.csv", low_memory=False)
cross = (pd.DataFrame({"Region": info.reset_index().drop_duplicates("sample_id")["sample_id"]}).iloc[:0])
meta_full = pd.read_excel(
    RAW/"primary/Information of sequenced samples_update_full878_filter819.xlsx",
    sheet_name="mfas5_819samples_phenSet4",
    usecols=["Region","SaleemNetworks"],
).dropna().drop_duplicates()
cross = meta_full.merge(ann[["Region","Full_name","Dictionary_lobe","Regional_map"]].drop_duplicates("Region"),
                        on="Region",how="left")
cross.rename(columns={"Region":"bo2023_region_id","SaleemNetworks":"locked_network",
                      "Full_name":"saleem_style_full_name","Dictionary_lobe":"dictionary_lobe",
                      "Regional_map":"broad_regional_map"}).sort_values(
                          ["locked_network","bo2023_region_id"]).to_csv(
                              OUT/"P2_Bo2023_Saleem_crosswalk.csv",index=False)

# Frozen-tier diagnostic bounds, not retrained models.
ablation = {
 "full_route":{"network_top3":91.94,"group_top3":72.48,"exact_top3":45.21},
 "remove_exact_output":{"network_top3":91.94,"group_top3":72.48,"exact_top3":None},
 "remove_group_output":{"network_top3":91.94,"group_top3":None,"exact_top3":45.21},
 "remove_all_fine_outputs":{"network_top3":91.94,"group_top3":None,"exact_top3":None},
 "beam_error_attribution":{
   "exact_misses":446,"network_beam_misses_within_exact_denominator_approx":66,
   "beam_share_of_exact_misses_approx":66/446,
   "exact_top3_given_beam_hit":49.07,
   "group_top3_given_beam_hit":78.32,
   "recovery_after_beam_miss":0.0,
   "note":"post hoc frozen-route diagnostic; not a retrained no-beam model"
 }
}

# Cell-type overlap using the Bo2023-provided marker list and panel symbols.
markers=pd.read_csv(RAW/"auxiliary_buildkit/12_celltype_high_expression_markers.csv")
panel_genes=set(panel.gene_symbol.astype(str))
background=set(pd.read_csv(ROOT/"code/reproducibility/p0_bio3_projector/vsd_gene_symbol_mapping_audit.csv").gene_symbol.astype(str))
cell_rows=[]
types=sorted(markers["Cell type"].dropna().astype(str).unique())
pvals=[]
for ct in types:
    g=set(markers.loc[markers["Cell type"].astype(str).eq(ct),"Gene"].dropna().astype(str)) & background
    k=len(panel_genes & g); M=len(background); n=len(g); N=len(panel_genes & background)
    p=float(hypergeom.sf(k-1,M,n,N)) if M and n else 1.0
    pvals.append(p); cell_rows.append({"cell_type":ct,"overlap":k,"marker_genes_in_background":n,"panel_in_background":N,"background":M,"p_raw":p})
order=np.argsort(pvals); q=np.ones(len(pvals)); running=1.0
for rankpos in range(len(order)-1,-1,-1):
    idx=order[rankpos]; running=min(running,pvals[idx]*len(pvals)/(rankpos+1)); q[idx]=running
for r,qq in zip(cell_rows,q): r["bh_q"]=float(qq)
pd.DataFrame(cell_rows).sort_values("bh_q").to_csv(OUT/"P2_panel_celltype_overlap.csv",index=False)

# Attempt exploratory GO/KEGG via g:Profiler with an explicit custom background.
enrichment_status={"status":"not_run"}
try:
    import requests
    payload={"organism":"hsapiens","query":sorted(panel_genes),"background":sorted(background),
             "sources":["GO:BP","KEGG"],"user_threshold":0.05,
             "significance_threshold_method":"fdr","no_evidences":True,"domain_scope":"custom"}
    resp=requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",json=payload,timeout=60)
    resp.raise_for_status()
    result=resp.json().get("result",[])
    pd.DataFrame(result).to_csv(OUT/"P2_panel_GO_KEGG_gprofiler.csv",index=False)
    enrichment_status={"status":"completed","query_n":len(panel_genes),"background_n":len(background),
                       "significant_terms":len(result),"retrieved_on":"2026-07-30","sources":["GO:BP","KEGG"]}
except Exception as e:
    enrichment_status={"status":"failed","error":str(e)}

audit={"scope":{"report_claimed_p2":28,"explicitly_listed":15},
       "rf200":rf_summary,"f1_stats":f1_stats,"ablation":ablation,
       "cell_type_significant_bh05":[r for r in cell_rows if r["bh_q"]<.05],
       "go_kegg":enrichment_status,
       "biomart_recheck":{"frozen_query_release":"not archived","frozen_query_date":"not archived",
          "current_recheck_release":"Ensembl 116 (June 2026)","recheck_date":"2026-07-30",
          "mapping_replaced":False},
       "requirements_roles":{"requirements.txt":"minimal app","requirements_reproducible.txt":"full manuscript reproduction",
          "requirements-repro.txt":"testing extras that includes requirements_reproducible.txt"},
       "live_demo":{"application_rate_limiter":"none found","application_upload_limit":"none set in code",
          "deployment_note":"host-level Streamlit limits apply; benchmarked raw-count input 153 samples x 28,415 genes, peak working set 222.0 MiB"}}
(OUT/"P2_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
print(json.dumps(audit,indent=2))
