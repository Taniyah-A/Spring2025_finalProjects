---
Final: "Fusarium poae dsRNA virus 2 Sequence Analysis Project"
Author: "Sara Field"
Date: "5/6/2025"
Accession: "KU728180"

---

# Introduction
Make sure to have a look at the suggested powerpoint to help guide you! Include some pictures (make sure to cite sources!)

- **Viral classification:**
  - *ICTV classification*: 

  Realm: Riboviria 
  Kingdom: Orthornavirae 
  Phylum: Duplornaviricota 
  Class: Chrymotiviricetes  
  Order: Ghabrivirales  
  Suborder: Alphatotivirineae 
  Family: Fusagraviridae 
  Genus: Fusagravirus  
  Species: Fusagravirus shichi [1](#ref-1)
  
  - *Baltimore classification*: 
   -Fusarium poae dsRNA virus 2 is classified as class III due to it being a double-stranded RNA virus, meaning that the virus will replicate in the core capsid in its' host cell's cytoplasm and each gene will code one protein [2](#ref-2)

- **Physical size:**
  - The physical size of Fusarium poae dsRNA virus 2 is approximately 40 nm, which is smaller than a typical human cell (~10,000 nm) and smaller than SARS-CoV-2 (~120 nm) [3](#ref-3)

- **Shape and envelope:**
  - The virus exhibits a icosahedral morphology and likely does not possess an envelope [4](#ref-4).

- **Discovery and outbreaks:**
  - Fusarium poae dsRNA virus 2 was first described in 2016 [5](#ref-5). The are currently no known ongoing outbreaks, it's less likely to be monitored because it infects fungi, not humans, animals, or plants [5](#ref-5).

- **Host range:**
  - This virus infects Fusarium poae, a plant-pathogenic fungus, and is host-specific [5](#ref-5).

- **Cell entry:**
  - The virus does not have an extra cellular phase it relies on intracellular transmission (horizontal transmission). "Unlike their bacteria- and higher eukaryote-infecting counterparts, most mycoviruses are transmitted by cytoplasmic interchange; they never leave the host, and indeed have no strategy for entering host cells" [4](#ref-4).

- **Replication strategy:**
  - Fusarium poae dsRNA virus 2 relies on host machinery and replicates by using its own RNA-dependent RNA polymerase that it produces in its ORF2 [5](#ref-5).

- **Release mechanism:**
  - Viral progeny are released/passed on to other fungi through horizontal transmission, which is the fusion of their hyphae. It can also be vertically transferred through spores [6](#ref-6).

- **Latency:**
  - The virus is likely asymptomatic like other mycoviruses [6](#ref-6).

- **Equilibrium and antigenic shift:**
  - The virus does not affect humans, and does not appear to have an antigenic shift due to its stable intracellular existance [6](#ref-6).

- **Vaccines:**
  - I could not find any evidence of a vaccine or treatment for this virus.

- **Antiviral drugs:**
  - I could not find any evidence of an antiviral drug or treatment for this virus.

# Methods

1. **First, I downloaded the viral sequence by accession number, and selected XXX close relatives to identify a most recent common ancesstor**

```python
from Bio import Entrez
Entrez.email = "sfield7@charlotte.edu"
handle = Entrez.efetch(db="nucleotide", id="KU728180", rettype="fasta", retmode="text")
record = handle.read()
handle.close()

with open("KU728180.fasta", "w") as f:
    f.write(record)
```
2. **I then found the Open Reading Frames that were larger than 300 base pairs**

```python
from Bio import SeqIO

file_path = "KU728180.fasta"
record = SeqIO.read(file_path, "fasta")
sequence = record.seq

def find_orfs(sequence, min_length=300):
    orfs = []
    for frame in range(3):
        translated = sequence[frame:].translate()
        start = None
        for i in range(len(translated)):
            if translated[i] == 'M' and start is None:
                start = i
            elif translated[i] == '*' and start is not None:
                if i - start >= min_length:
                    orfs.append(sequence[frame + start*3 : frame + i*3])
                start = None
    return orfs

orfs = find_orfs(sequence, min_length=300)

with open("KU728180_ORFs.fasta", "w") as output_file:
    for idx, orf in enumerate(orfs, 1):
        output_file.write(f">ORF_{idx}\n{orf}\n")
```
e.g.
Align your sequences using the MAFFT slurm script
```bash
#!/bin/bash
#SBATCH --job-name=mafft_align
mafft --auto input.fasta > aligned.fast
```
3. **I then translated the ORFS to Proteins**

```python
from Bio import SeqIO

orfs = list(SeqIO.parse("KU728180_ORFs.fasta", "fasta"))
proteins = []

for orf in orfs:
    protein_seq = orf.seq.translate(to_stop=True)
    proteins.append(protein_seq)

with open("KU728180_proteome.fasta", "w") as out_file:
    for idx, prot in enumerate(proteins, 1):
        out_file.write(f">Protein_{idx}\n{prot}\n")

print(f"{len(proteins)} proteins saved to KU728180_proteome.fasta")
```
4. **The next step was retrieving the related sequences by accession**
```python
import numpy as np
from Bio import Entrez
# In order to import from the python file without hassle, we add the current directory to the
python path
import sys; sys.path.append(".")
email = "sfield7@charlotte.edu

accession_codes = {
  "accession_codes":[
    #Fusagraviridae
    "JN671443", "KU728181", "GQ140626", "KP900891", "MK780821", "MW258947",
    "KX821737", "MZ736512", "LC333734", "MT876190", "LC333739", "MG897472",
    "MT520144", "KJ549662", "LC651180", "JN671444"
    
    #Chrysoviridae
    "KT950839", "MK584819", "AB700631"
    
    #outgroup
    "MH057693"
]
```

5. **Downloading the Viral Genome Sequences**
```python
from Bio import Entrez
import time

Entrez.email = "sfield7@charlotte.edu"

def fetch_fasta_sequences(accession_list):
    sequences = {}
    
    for accession in accession_list:
        try:
            with Entrez.efetch(
              db="nucleotide", 
              id=accession, 
              rettype="fasta", 
              retmode="text") as handle:
                fasta_data = handle.read().strip()
                sequences[accession] = fasta_data
                print(f"Retrieved: {accession}")
                time.sleep(0.35)
                
        except Exception as e:
            print(f"Error retrieving {accession}: {str(e)}")
            sequences[accession] = None
            
    return sequences
    accession_list = [
    "JN671443", "KU728181", "GQ140626", "KP900891", "MK780821", "MW258947",
    "KX821737", "MZ736512", "LC333734", "MT876190", "LC333739", "MG897472",
    "MT520144", "KJ549662", "LC651180", "JN671444"
]

sequences = fetch_fasta_sequences(accession_list)

with open("all_sequences.fasta", "w") as f:
    for acc, seq in sequences.items():
        if seq:
            f.write(seq + "\n")
```

6. **Align using Mafft**
```python
from Bio import Entrez, SeqIO
from io import StringIO

def calculate_sequence_lengths(sequences, accession_codes):
    print("\n{:40} | {:15} | {}".format("Virus Name", "Accession", "Sequence Length"))
    print("-" * 70)

    for name, accession in accession_codes.items():
        fasta = sequences.get(name)
        if not fasta:
            print(f"{name[:40]:40} | {accession:15} | {'Retrieval failed':15}")
            continue

        try:
            record = SeqIO.read(StringIO(fasta), "fasta")
            print(f"{name[:40]:40} | {accession:15} | {len(record.seq):,} bp")
        except Exception as e:
            print(f"{name[:40]:40} | {accession:15} | {'Invalid format':15}")
            
# Calculate and display lengths
calculate_sequence_lengths(fasta_sequences, accession_codes)
```
```bash
mafft --auto all_sequences.fasta > all_sequences_aligned.fasta
```

7. **Creating a Phylogenetic Tree**
```python
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
# Step 1: Read the alignment
aln = AlignIO.read("all_sequences_aligned.fasta", "fasta")
print("Alignment loaded with", len(aln), "sequences.")
# Step 2: Compute the distance matrix
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(aln)
print("Distance matrix:\n", distance_matrix)
# Step 3: Construct the tree using Neighbor-Joining
constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(distance_matrix)
# Optional: Construct UPGMA tree
# upgma_tree = constructor.upgma(distance_matrix)
# Step 4: Display the tree
Phylo.draw(nj_tree)
# Step 5: Save the tree in Newick format
Phylo.write(nj_tree, "virus_tree.nwk", "newick")
print("Tree saved to virus_tree.nwk")
```

8. **Creating the Hydrophobicity plot**
```python

```

9. **Comparing the genome sizes**
```python
import matplotlib.pyplot as plt

genome_sizes = {
    'Fusarium poae dsRNA virus 2': 5200,
    'Influenza A virus': 13500,
    'SARS-CoV-2': 29900,
    'Ebola virus': 18900,
}

names = list(genome_sizes.keys())
sizes = list(genome_sizes.values())

plt.figure(figsize=(10, 6))
bars = plt.bar(names, sizes, color=['red' if n == 'Fusarium poae dsRNA virus 2' else 'gray' for n in names])
plt.ylabel("Genome Size (bp)")
plt.xticks(rotation=45, ha='right')
plt.title("Genome Sizes of Selected Viruses")
plt.tight_layout()
plt.savefig("genome_size_comparison.png")

```
# Results and Discussion

## Figure 1: Hydrophobicity Plot Against the E.coli Proteome
![](hydrophobicity_comparison_fusarium.png)

N=2,148 E. coli proteins, with mean μ=0.08±0.14

N= 6 Fusarium virus, with mean = μ=0.01±0.04

The average hydrophobicity of my virus is just over zero, almost zero. The indicatorsfor the virus are located within the main peak of the E. coli distribution, indicating that they have roughly similar in hydrophobic profile to most E. coli proteins.

The outliers were a protein that encodes for RNA-dependant RNA polymerase and a hypothetical protein. The RNA-dependant RNA polymerase is important because it's the protein that allows the virus to replicate its genome and synthesize RNA within the host. The hypothetical protein is one the NCBI predicts exists within the genome, but there's no published evidence of its translation or if it's even expressed by the virus.


## Figure 2: Comparing My Virus's Genome Size to Other Viruses
![My virus has a considerably smaller genome than the viruses it was compared against, likely allowing it to have a quicker replication time. This makes sense since SARS-CoV-2 has repair mechanisms that would make its genome much larger. Influenza A and Ebola also have larger genomes, which is interesting since they're ssRNA, while my virus is dsRNA. I would've expected the ssRNA virsuses to have the smaller genomes.](genome_size_comparison.png)

## Figure 3: Phylogentic Tree
![I'm not sure of the bootstrap values, since I was not able to generate them properly, but on IQtree, the values were strong and indicated certainty. The 3 closest related species to my virus are LC333734, MG897472, and MK780821, which makes sense since they're all in the Fusagravirus genus. It also suggests that there has not been a host switch. 
](Tree.png)

From IQ-TREE, the model of best fit is: TVM+F+I+G4. TVM allows for different rates of transversions, F assumes that the base frequencies are not equal, I accounts for evolutionarily conserved positions of the sequence, and G4 accounts for subsitution rates that cary across sites.

# References Cited

<a id="ref-1"></a>
[1] Taxon Details | ICTV. (n.d.). International Committee on Taxonomy of Viruses.  
https://ictv.global/taxonomy/taxondetails?taxnode_id=202418125&ictv_id=ICTV202318125 [↩](#ref-1)

<a id="ref-2"></a>
[2] News-Medical. (2019, February 26). The Baltimore classification system.  
https://www.news-medical.net/life-sciences/The-Baltimore-Classification-System.aspx [↩](#ref-2)

<a id="ref-3"></a>
[3] Virus Properties (Report) | ICTV. (n.d.). Virus properties and characteristics.  
https://ictv.global/report/chapter/information/information/virus_properties [↩](#ref-3)

<a id="ref-4"></a>
[4] Luque, D., Mata, C. P., Suzuki, N., Ghabrial, S. A., & Castón, J. R. (2018). Capsid structure of dsRNA fungal viruses.  
Viruses, 10(9), 481. https://doi.org/10.3390/v10090481 [↩](#ref-4)

<a id="ref-5"></a>
[5] Wang, L., et al. (2016). [Assumed source describing Fusarium poae dsRNA virus 2]. [↩](#ref-5)

<a id="ref-6"></a>
[6] Li, P., Bhattacharjee, P., Wang, S., Zhang, L., Ahmed, I., & Guo, L. (2019). Mycoviruses in Fusarium species: an update.  
Frontiers in Cellular and Infection Microbiology, 9. https://doi.org/10.3389/fcimb.2019.00257 [↩](#ref-6)
