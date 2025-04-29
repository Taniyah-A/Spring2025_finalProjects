---
Final: "[Insert Virus Name Here] Sequence Analysis Project"
Author: "Your Name"
Date: "[Change this date]"
Accession: "Accession #"

---

# Introduction
Make sure to have a look at the suggested powerpoint to help guide you! Include some pictures (make sure to cite sources!)

- **Viral classification:**  
  - *ICTV classification*: [Insert ICTV classification here with citation].  
  - *Baltimore classification*: [Insert Baltimore group and explanation here].

- **Physical size:**  
  - The physical size of [your virus] is approximately XX nm, which is [smaller/larger] than a typical human cell (~10,000 nm) and [smaller/larger] than SARS-CoV-2 (~120 nm) [insert citation].

- **Shape and envelope:**  
  - The virus exhibits a [describe shape, e.g., icosahedral, helical, complex] morphology and [does/does not] possess an envelope [insert citation].

- **Discovery and outbreaks:**  
  - [Your virus] was first described in [year] [insert citation]. The most recent outbreak occurred in [year/location] [insert citation].

- **Host range:**  
  - This virus infects [host(s)], and is [host-specific/generalist] [insert citation].

- **Cell entry:**  
  - The virus penetrates the host cell via [mechanism, e.g., endocytosis, membrane fusion, receptor-mediated entry] [insert citation].

- **Replication strategy:**  
  - [Your virus] [has its own replication machinery/relies on host machinery] and replicates by [describe process] [insert citation].

- **Release mechanism:**  
  - Viral progeny are released by [cell lysis, budding, excretion, etc.] [insert citation].

- **Latency:**  
  - [State if the virus shows latency and describe, or state if not] [insert citation].

- **Equilibrium and antigenic shift:**  
  - [Discuss if the virus is in equilibrium with humans and if it shows antigenic shift] [insert citation].

- **Vaccines:**  
  - [State if vaccines are available, type/mechanism of vaccine] [insert citation].

- **Antiviral drugs:**  
  - [List available antivirals and their mechanisms] [insert citation].

# Methods

1. **First, I downloaded the viral sequence by accession number, and selected XXX close relatives to identify a most recent common ancesstor**  Use the provided accession number (Column X) to download the viral genome sequence. Show the python code for each of the functions you performed.

```python
#Insert python code after each step
```
...Follow the additional methods instructions in the text. Make sure to include your code where relevant! For example, if you used a slurm script, make sure to paste it in here.

e.g. 
Align your sequences using the MAFFT slurm script
```bash
#!/bin/bash
#SBATCH --job-name=mafft_align
mafft --auto input.fasta > aligned.fast
```
     

# Results and Discussion

Present your results, including tables, figures, and summary statistics. This should be written more like a text document. Make sure you follow the instructions about which plots to paste within here and reference them as figures (e.g., Fig. 1) when you are describing them.

(10 points) Hydrophobicity plot agains the E.coli proteome
(10 points) Identify any of your outlier hydrophobic proteins and BLAST them. What did these sequences annotate as? Might they have any important function to viral entry? 
(10 points) Plot the genome size of your virus relative to other viruses (see code from lab 9.12b). 
(20 points) Phylogeny and model selection. Use figtree to root the tree to your outgroup, make it look nice by ordering the nodes, increasing tip label size. State and interpret the best fit model used to infer this phylogeny. Display the bootstrap values. Discuss about your results. What are the three closest relatives of your virus, does it suggest a host switch? Are the branches well-supported by bootstrap values?


# References Cited

Use at least five references, formatted using a numbering system (1) or (Doe, et al. 2024). List the references in the appropriate bibliography format here. Use markdown to have all of your text references link to the bottom part of the page to where the reference is listed.
