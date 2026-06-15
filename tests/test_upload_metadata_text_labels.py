from io import BytesIO

from app.pages.upload_page import _read_uploaded_file
from data_processor import DataProcessor


class NamedUpload(BytesIO):
    name = "sample.csv"


def test_uploaded_numeric_region_label_is_preserved_as_text() -> None:
    uploaded = NamedUpload(
        b"gene_symbol,tpm_value,sample_id,ground_truth_region\n"
        b"GENE1,1.2,SAMPLE001,44563\n"
        b"GENE2,1.4,,\n"
    )

    df = _read_uploaded_file(uploaded)
    metadata = DataProcessor().extract_embedded_metadata(df)

    assert metadata["sample_id"] == "SAMPLE001"
    assert metadata["ground_truth_region"] == "44563"
