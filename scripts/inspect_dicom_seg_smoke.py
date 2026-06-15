#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor" / "dicom_tools_py314"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import pydicom  # noqa: E402


def main() -> int:
    root = ROOT / "data" / "tcia_tcga_glioma_mri" / "idc_seg_smoke"
    for path in sorted(root.rglob("*.dcm")):
        ds = pydicom.dcmread(path, stop_before_pixels=False)
        print("\nFILE", path)
        print("Modality", getattr(ds, "Modality", ""))
        print("SeriesDescription", getattr(ds, "SeriesDescription", ""))
        print("SOPClassUID", getattr(ds, "SOPClassUID", ""))
        print("Rows/Columns", getattr(ds, "Rows", ""), getattr(ds, "Columns", ""))
        print("NumberOfFrames", getattr(ds, "NumberOfFrames", ""))
        print("FrameOfReferenceUID", getattr(ds, "FrameOfReferenceUID", ""))
        if hasattr(ds, "SegmentSequence"):
            print("Segments", len(ds.SegmentSequence))
            for seg in ds.SegmentSequence:
                print(
                    "  ",
                    getattr(seg, "SegmentNumber", ""),
                    getattr(seg, "SegmentLabel", ""),
                    getattr(seg, "SegmentDescription", ""),
                    getattr(seg, "SegmentAlgorithmName", ""),
                )
        print("HasPixelData", hasattr(ds, "PixelData"), "PixelDataBytes", len(getattr(ds, "PixelData", b"")))
        if hasattr(ds, "PerFrameFunctionalGroupsSequence"):
            print("PerFrameFunctionalGroups", len(ds.PerFrameFunctionalGroupsSequence))
        try:
            arr = ds.pixel_array
            print("pixel_array", arr.shape, arr.dtype, "nonzero", int((arr > 0).sum()))
        except Exception as exc:
            print("pixel_array_error", type(exc).__name__, str(exc)[:300])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
