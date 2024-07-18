# Chemogenomic Library using Open-Access Data

# Chemogenomics Project

Chemogenomics employs well-characterized tool compounds to functionally annotate proteins within complex cellular systems, facilitating target discovery and validation. Unlike highly selective chemical probes, the small molecule modulators (e.g., agonists or antagonists) used in chemogenomic studies do not require exclusive selectivity. Due to the limited availability of high-quality chemical probes for a small fraction of targets, the broader criteria for defining small molecules in chemogenomics allow for a more extensive target coverage.

Open-access data is critical in deciphering small molecules and targets used in drug discovery. By utilizing trusted sources that employ rigid criteria in defining small molecules, one may have confidence that the resulting chemogenomic library comprises high-quality compounds. These compounds can then be leveraged to explore various biological targets, providing valuable insights into protein functions and interactions. This approach accelerates the identification and validation of drug targets and enhances the reproducibility and reliability of research findings, ultimately contributing to more effective drug development processes.

## Project Overview

The project housed in this repository is a fraction of the work completed over six months, during which I created a chemogenomic library from overlapping molecules of credible datasets. The datasets used were:
- SGC Donated Chemical Probes
- Chemical Probes Portal
- Drug Central 2023
- EUbOpen

Using RDKit, I ensured that the compound representations of the datasets were consistent and then created a merged dataset of the overlaps. One difficulty I encountered was the different naming conventions (Canonical SMILES, InChI string, InChI Key) and ensuring that the same compounds were recognized as equivalent. I implemented a standardization protocol to address this, converting all representations to a unified InChI string format to facilitate accurate comparisons and merging. Additionally, I incorporated cross-referencing techniques and manual verification steps to ensure the integrity of the compound matches.

## Data Query and Normalization

Using SQL, I queried the dataset in GOSTAR, a third-party application similar to ChEMBL, known for its diverse bioactivity data. This diversity is beneficial from a chemogenomic library standpoint, as it broadens the range of biological activities. However, the data's varied sources present challenges in ensuring consistency and reliability. Despite GOSTAR's efforts to reconcile these differences, additional work was needed. To further harmonize the data, I applied several normalization and validation techniques. This included standardizing measurement units (uM) and activity types (IC50, AC50, EC50, Ki, Kd) and ensuring that experimental conditions were comparable. I also employed statistical methods to identify and address outliers or inconsistencies within the dataset. By doing so, I aimed to enhance the quality and coherence of the integrated chemogenomic library.

## Compound Acquisition

Using the library, I determined which compounds were not in UCB's possession and created a rationale for acquiring each, seeking second opinions from chemists. This step involved assessing the novelty, potential therapeutic relevance, and synthetic feasibility of each compound. We also considered factors such as purity, stability, and cost-effectiveness/commercial availability to acquire high-quality and practical candidates for further research. The combination of expert opinions and commercial viability checks highlighted the importance of a multidisciplinary approach in chemogenomics.

## Similarity Scoring and Visualization

I determined the Morgan fingerprints for each of UCB's compounds and those suggested for acquisition. Morgan fingerprints are binary bit representations of chemical structures. I employed Tanimoto similarity scoring to determine the similarity among the groups of compounds and then used KNN to determine, for every candidate, the five most similar UCB compounds. I visualized this relationship using NetworkX.

After determining which compounds had a similarity score over 0.8, I manually reviewed the targets against which the candidates were tested. Based on the assumption that structure relates to function, I suggested that UCB compounds with high similarity scores be tested against those same targets. This approach leverages the structural similarity to infer potential bioactivity, aiming to uncover new therapeutic potentials for UCB's existing compounds.

## Conclusion

The availability of extensive public datasets allows for the implementation of advanced machine learning and data mining techniques, further enhancing the capability to predict and identify novel drug candidates. This project demonstrates that, despite data integration and standardization challenges, the resulting comprehensive and well-curated chemogenomic library can significantly streamline the drug discovery process. By providing a centralized and reliable source of annotated compounds, this library facilitates the identification of new therapeutic targets, understanding protein functions, and exploring complex biological pathways. The insights gained from this resource can accelerate the development of effective and targeted therapies, ultimately contributing to improved patient outcomes and addressing unmet medical needs.
