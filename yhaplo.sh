bcftools view -r Y,chrY cohort.vcf.gz -Ov - \
  | bcftools norm -m -any -f hg19.no_alt.fa - \
  | bcftools view -f PASS -Oz - \
  | bcftools annotate -h <(echo '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depths">') -Oz - \
  > filt_Y_ad.vcf.gz

bcftools index filt_Y_ad.vcf.gz

python3 pre_yleaf.py

bcftools index filt_Y_ad_final.vcf.gz

mkdir -p temp_yleaf

Yleaf -vcf filt_Y_ad_final.vcf.gz \
      -o temp_yleaf \
      --reference_genome hg19 \
      -pq 0.8 \
      -t 4 \
      --force

rm filt_Y_ad*

mv temp_yleaf/hg_prediction.hg .
rm -rf temp_yleaf

python3 yhaplo_visual.py
