#!usr/bin/python
# This script is to get epitope candidate from vep-vcf.
# Author : KJM
# Usage : python this.script vep.vcf epitope_length > output.file


from vaxrank.manufacturability import ManufacturabilityScores
from vaxrank.vaccine_peptide import *
import sys

class PvacpeptideVaccinePeptide(VaccinePeptide):
    def __new__(cls, peptide):
        return VaccinePeptideBase.__new__(
            cls,
            mutant_protein_fragment="",
            mutant_epitope_predictions="",
            wildtype_epitope_predictions="",
            mutant_epitope_score="",
            wildtype_epitope_score="",
            num_mutant_epitopes_to_keep="",
            manufacturability_scores=ManufacturabilityScores.from_amino_acids(peptide)
        )

def manufacturability_headers(self):
    return[
            'cterm_7mer_gravy_score',
            'max_7mer_gravy_score',
            'difficult_n_terminal_residue',
            'c_terminal_cysteine',
            'c_terminal_proline',
            'cysteine_count',
            'n_terminal_asparagine',
            'asparagine_proline_bond_count']

def append_manufacturability_metrics(self, line, peptide):
    metrics = peptide.manufacturability_scores
    line['cterm_7mer_gravy_score'] = metrics.cterm_7mer_gravy_score
    line['max_7mer_gravy_score'] = metrics.max_7mer_gravy_score
    line['difficult_n_terminal_residue'] = metrics.difficult_n_terminal_residue
    line['c_terminal_cysteine'] = metrics.c_terminal_cysteine
    line['c_terminal_proline'] = metrics.c_terminal_proline
    line['cysteine_count'] = metrics.cysteine_count
    line['n_terminal_asparagine'] = metrics.n_terminal_asparagine
    line['asparagine_proline_bond_count'] = metrics.asparagine_proline_bond_count

    return line

anno_info='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|TSL|APPRIS|SIFT|PolyPhen|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|DownstreamProtein|ProteinLengthChange|WildtypeProtein'
anno_info=anno_info.split('|')

protein_pos=anno_info.index('Protein_position')
amino_index=anno_info.index('Amino_acids')
HGVSp=anno_info.index('HGVSp')
DS=anno_info.index('DownstreamProtein')
WT=anno_info.index('WildtypeProtein')
Gene=anno_info.index('SYMBOL')


amino_acid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



fr=open(sys.argv[1],'r')
epitope_len=int(sys.argv[2])

text='chrNum\tcoordinate\tRef\tAlt\tConsequence\tGene\tPeptidePos\tChange\tWT\tMT\n'

for i in fr.readlines():
	if i.startswith('#'):
		pass
	else:
		ii=i.strip().split('\t')
		for info in ii[7].split(';'):
			if info.startswith('CSQ'):
				vep=info.split('|')
				Consequence=vep[1]
				WT_seq=vep[WT]
				if 'missense' in Consequence or 'inframe_insertion' in Consequence:
					wildtype_amino_acid, mutant_amino_acid = vep[amino_index].split('/')
					if '*' in wildtype_amino_acid:
						wildtype_amino_acid = wildtype_amino_acid.split('*')[0]
					elif 'X' in wildtype_amino_acid:
						wildtype_amino_acid = wildtype_amino_acid.split('X')[0]
					if '*' in mutant_amino_acid:
						mutant_amino_acid = mutant_amino_acid.split('*')[0]
						stop_codon_added = True
					elif 'X' in mutant_amino_acid:
						mutant_amino_acid = mutant_amino_acid.split('X')[0]
						stop_codon_added = True
					else:
						stop_codon_added = False
					if wildtype_amino_acid == '-':
						position = int(vep[protein_pos].split('-', 1)[0])
						wildtype_amino_acid_length = 0
					else:
						if '-' in vep[protein_pos]:
							position = int(vep[protein_pos].split('-', 1)[0]) - 1
							wildtype_amino_acid_length = len(wildtype_amino_acid)
						else:
							position = int(vep[protein_pos]) - 1
							wildtype_amino_acid_length = len(wildtype_amino_acid)

					if len(WT_seq) + len(mutant_amino_acid) - len(wildtype_amino_acid) < epitope_len:
						pass
					else:
						already=[]
						for p in range(position-epitope_len+1,position): # +1 is to avoid wildtype mutation in front of mutation position
							if p < 0:
								pass
							else:
								MT_last=epitope_len-(position-p)
								WT_last=MT_last-len(mutant_amino_acid[:MT_last])
								if WT_last <= 0 and MT_last <= len(mutant_amino_acid):
									peptide=WT_seq[p:position]+mutant_amino_acid[:MT_last]
									if WT_seq[p:p+epitope_len] == peptide: # becase of repeated base ex) PALAAAAAAAAAAAAA  -> insert 3 position from L to LA and epitope lengtg 5 => WT and MT is same!
										pass
									else:
										if peptide in already:
											pass
										else:
											text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
											if len(WT_seq[p:p+epitope_len]) < epitope_len:
												text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+'No_WT'+'\t'+peptide+'\n'
											else:	
												text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WT_seq[p:p+epitope_len]+'\t'+peptide+'\n'
											already.append(peptide)
								elif WT_last > 0 and position+wildtype_amino_acid_length+WT_last <= len(WT_seq):
									peptide=WT_seq[p:position]+mutant_amino_acid[:MT_last]+WT_seq[position+wildtype_amino_acid_length:position+wildtype_amino_acid_length+WT_last] # +wildtype_amino_acid_length cause of there are two type '-/GA or G/GA', if '-' wildtype_amino_acid_length = 0, or wildtype_amino_acid_length = len('G/GA'.split('/')[0])
									if WT_seq[p:p+epitope_len] == peptide:
										pass
									else:
										if peptide in already:
											pass
										else:
											text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
											if len(WT_seq[p:p+epitope_len]) < epitope_len:
												text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+'No_WT'+'\t'+peptide+'\n'
											else:
												text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WT_seq[p:p+epitope_len]+'\t'+peptide+'\n'
											already.append(peptide)
						for p in range(len(mutant_amino_acid)): # if len(mutant_amino_acid) ==1, p is only '0'
							last=epitope_len-(len(mutant_amino_acid[p:p+epitope_len]))
							WTpep=WT_seq[position+p:position+p+epitope_len]
							if last <= 0:
								peptide=mutant_amino_acid[p:p+epitope_len]
								if WTpep == peptide:
									pass
								else:
									if peptide in already:
										pass
									else:
										text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
										if len(WTpep) < epitope_len:
											text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+'No_WT'+'\t'+peptide+'\n'
										else:
											text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WTpep+'\t'+peptide+'\n'
										already.append(peptide)
							elif position+wildtype_amino_acid_length+last <= len(WT_seq):
								peptide=mutant_amino_acid[p:p+epitope_len]+WT_seq[position+wildtype_amino_acid_length:position+wildtype_amino_acid_length+last]
								if WTpep == peptide:
									pass
								else:
									if peptide in already:
										pass
									else:
										text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
										if len(WTpep) < epitope_len:
											text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+'No_WT'+'\t'+peptide+'\n'
										else:
											text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WTpep+'\t'+peptide+'\n'

				elif 'inframe_del' in Consequence:
					wildtype_amino_acid, mutant_amino_acid = vep[amino_index].split('/')
					if '*' in wildtype_amino_acid:
						wildtype_amino_acid = wildtype_amino_acid.split('*')[0]
					elif 'X' in wildtype_amino_acid:
						wildtype_amino_acid = wildtype_amino_acid.split('X')[0]
					if '*' in mutant_amino_acid:
						mutant_amino_acid = mutant_amino_acid.split('*')[0]
						stop_codon_added = True
					elif 'X' in mutant_amino_acid:
						mutant_amino_acid = mutant_amino_acid.split('X')[0]
						stop_codon_added = True
					else:
						stop_codon_added = False
					position = int(vep[protein_pos].split('-', 1)[0]) - 1
					wildtype_amino_acid_length = len(wildtype_amino_acid)
					if mutant_amino_acid == '-':
						mutant_amino_acid = ''
					if len(WT_seq) - wildtype_amino_acid_length + len(mutant_amino_acid) < epitope_len:
						pass
					elif position < epitope_len-1:
						for p in range(position+1):
							last=p+epitope_len+wildtype_amino_acid_length
							if last > len(WT_seq):
								pass
							else:
								peptide=WT_seq[p:position]+mutant_amino_acid+WT_seq[position+wildtype_amino_acid_length:last]
					else:
						for p in range(position-epitope_len+1,position+1):
							last=p+epitope_len+wildtype_amino_acid_length-len(mutant_amino_acid)
							if last > len(WT_seq):
								pass
							else:
								peptide=WT_seq[p:position]+mutant_amino_acid+WT_seq[position+wildtype_amino_acid_length:last]
								if WT_seq[p:p+epitope_len] == peptide:
									pass
								else:
									text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
									text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WT_seq[p:p+epitope_len]+'\t'+peptide+'\n'
				elif 'frameshift_variant' in Consequence:
					already=[]
					wildtype_amino_acid, mutant_amino_acid = vep[amino_index].split('/')
					DW_seq=vep[DS]
					position = int(vep[protein_pos].split('-', 1)[0]) - 1
					if len(WT_seq) + len(DW_seq) < epitope_len:
						pass
					elif position < epitope_len:
						for p in range(position):
							last=epitope_len-(position-p)
							if last > len(DW_seq):
								pass
							else:
								peptide=WT_seq[p:position]+DW_seq[:last]
								if WT_seq[p:p+epitope_len] == peptide:
									pass
								else:
									if peptide in already:
										pass
									else:
										text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
										text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WT_seq[p:p+epitope_len]+'\t'+peptide+'\n'
										already.append(peptide)
					else:
						for p in range(position-epitope_len+1,position):
							last=epitope_len-(position-p)
							if last > len(DW_seq):
								pass
							else:
								peptide=WT_seq[p:position]+DW_seq[:last]
								if WT_seq[p:p+epitope_len] == peptide:
									pass
								else:
									if peptide in already:
										pass
									else:
										text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
										text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+WT_seq[p:p+epitope_len]+'\t'+peptide+'\n'
										already.append(peptide)
					if len(DW_seq) >= epitope_len:
						for p in range(len(DW_seq)):
							if int(p+epitope_len) > int(len(DW_seq)):
								pass
							else:
								peptide=DW_seq[p:p+epitope_len]
								if peptide in already:
									pass
								else:
									text+=ii[0]+'\t'+ii[1]+'\t'+ii[3]+'\t'+ii[4]+'\t'
									text+=Consequence+'\t'+vep[Gene]+'\t'+vep[protein_pos]+'\t'+vep[amino_index]+'\t'+'No_WT'+'\t'+peptide+'\n'
									already.append(peptide)



count=0
final_text=''
for i in text.strip().split('\n'):
	if count==0:
		final_text+=i+'\t'+'cterm_7mer_gravy_score'+'\t'+'max_7mer_gravy_score'+'\t'+'difficult_n_terminal_residue'+'\t'+'c_terminal_cysteine'+'\t'+'c_terminal_proline'+'\t'+'cysteine_count'+'\t'+'n_terminal_asparagine'+'\t'+'asparagine_proline_bond_count'+'\n'
		count+=1
	else:
		ii=i.split('\t')
		sequence=ii[-1]
		seq_num=ii[0]+'-'+ii[1]+'-'+ii[2]+'-'+ii[3]
		line = {'id': seq_num,'peptide_sequence': sequence}
		manufacturability_scores=ManufacturabilityScores.from_amino_acids(sequence)
		final_text+=i+'\t'+str(round(manufacturability_scores.cterm_7mer_gravy_score,3))+'\t'+str(round(manufacturability_scores.max_7mer_gravy_score,3))+'\t'+str(manufacturability_scores.difficult_n_terminal_residue)+'\t'+str(manufacturability_scores.c_terminal_cysteine)+'\t'+str(manufacturability_scores.c_terminal_proline)+'\t'+str(manufacturability_scores.cysteine_count)+'\t'+str(manufacturability_scores.n_terminal_asparagine)+'\t'+str(manufacturability_scores.asparagine_proline_bond_count)+'\n'

fw=open(sys.argv[1].replace('vcf',sys.argv[2]+'mer.epitope.tsv'),'w')
fw.write(final_text)
fw.close()
