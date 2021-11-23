# Generate_Epitope_Kmer
Generate Epitopes for Neoantigen from VCF

# Purpose
For neoantigen, there must be epitope candidates which have mutant peptide sequence from tumor cell. So using somatic vcf file, first get nonsynonymous variant like missnese, frameshift, indel and then From the mutation position, generate N-mer epitope peptide sequence from both sides. At final, it returns epitope candidate results with manufacturability.

# Usage
python Generate_EpitopeCandi.py {vcf annotated from VEP} {length of epitope to get}
