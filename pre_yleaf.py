import pysam

inp = "filt_Y_ad.vcf.gz"
out = "filt_Y_ad_final.vcf.gz"

with pysam.VariantFile(inp, "r") as vcf_in:
    out_header = vcf_in.header.copy()
    with pysam.VariantFile(out, "w", header=out_header) as vcf_out:
        for rec in vcf_in:
            for s in rec.samples.values():
                gt = s.get('GT')
                if gt is None or None in gt:
                    s['AD'] = (0, 0)
                else:
                    ref = sum(1 for a in gt if a == 0)
                    alt = sum(1 for a in gt if a == 1)
                    s['AD'] = (ref * 5, alt * 5)
            vcf_out.write(rec)