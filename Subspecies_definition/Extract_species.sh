ml BCFtools/1.10.2-GCC-9.3.0
ml VCFtools/0.1.16-GCC-10.2.0
VCF=/data/umcg-tifn/SNV/merged_raw/lld1.lld2.300ob.500fg.ibd.vcf.gz


Species_id=$1 #102293
temp_bcf=temp/$Species_id\.bcf
Genotypes=Genotypes/genotype_$Species_id\.tsv
Allele_frequencies=Genotypes/AF_$Species_id\.tsv

echo "Extracting species"
if [ ! -f $temp_bcf ]; then
	bcftools view $VCF --regions $Species_id  -O b -o $temp_bcf
fi
echo "Encoding genotypes"
if [ ! -f $Genotypes\.012 ]; then
	vcftools --bcf $temp_bcf --012  --out $Genotypes
fi
echo "Generate table with per sample allele frequencies"
bcftools query -f '%CHROM\t%POS\t[ %AF_ALT\t]\n' $temp_bcf  > $Allele_frequencies
echo "Cleaning temporal files"


