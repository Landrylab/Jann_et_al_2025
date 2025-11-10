## Barcode Design

Barcodes were generated using the **DNAbarcodes** package  
<https://github.com/lzamparo/DNAbarcodes>.

1. The script `generate_barcodes.py` was used to enumerate every possible **12-nucleotide** barcode sequence.

2. The script `validate_barcodes.py` was used to select **30,000** barcodes meeting the following criteria:
   - GC content between **40% and 60%**
   - Minimum **Hamming distance ≥ 3** between any two barcodes
   - Barcode length: **12 bp**

**Example command:**
```bash
validate_barcodes.py \
    -input all_possible_barcodes_12bp.txt \
    -output 30000_barcode_3mismatch_gc_04_06.txt \
    30000 12 0.4 0.6 3
