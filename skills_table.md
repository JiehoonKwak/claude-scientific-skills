# Scientific Skills Overview

This document contains tables of scientific skills organized by category.

## Scientific Databases

| Name | Description |
|------|-------------|
| alphafold-database | Access AlphaFold's 200M+ AI-predicted protein structures. Retrieve structures by UniProt ID, download PDB/mmCIF files, analyze confidence metrics (pLDDT, PAE), for drug discovery and structural biology. |
| biorxiv-database | Efficient database search tool for bioRxiv preprint server. Use this skill when searching for life sciences preprints by keywords, authors, date ranges, or categories, retrieving paper metadata, downloading PDFs, or conducting literature reviews. |
| chembl-database | Query ChEMBL's bioactive molecules and drug discovery data. Search compounds by structure/properties, retrieve bioactivity data (IC50, Ki), find inhibitors, perform SAR studies, for medicinal chemistry. |
| clinicaltrials-database | Query ClinicalTrials.gov via API v2. Search trials by condition, drug, location, status, or phase. Retrieve trial details by NCT ID, export data, for clinical research and patient matching. |
| clinpgx-database | Access ClinPGx pharmacogenomics data (successor to PharmGKB). Query gene-drug interactions, CPIC guidelines, allele functions, for precision medicine and genotype-guided dosing decisions. |
| clinvar-database | Query NCBI ClinVar for variant clinical significance. Search by gene/position, interpret pathogenicity classifications, access via E-utilities API or FTP, annotate VCFs, for genomic medicine. |
| cosmic-database | Access COSMIC cancer mutation database. Query somatic mutations, Cancer Gene Census, mutational signatures, gene fusions, for cancer research and precision oncology. Requires authentication. |
| ena-database | Access European Nucleotide Archive via API/FTP. Retrieve DNA/RNA sequences, raw reads (FASTQ), genome assemblies by accession, for genomics and bioinformatics pipelines. Supports multiple formats. |
| ensembl-database | Query Ensembl genome database REST API for 250+ species. Gene lookups, sequence retrieval, variant analysis, comparative genomics, orthologs, VEP predictions, for genomic research. |
| fda-database | Query openFDA API for drugs, devices, adverse events, recalls, regulatory submissions (510k, PMA), substance identification (UNII), for FDA regulatory data analysis and safety research. |
| gene-database | Query NCBI Gene via E-utilities/Datasets API. Search by symbol/ID, retrieve gene info (RefSeqs, GO, locations, phenotypes), batch lookups, for gene annotation and functional analysis. |
| geo-database | Access NCBI GEO for gene expression/genomics data. Search/download microarray and RNA-seq datasets (GSE, GSM, GPL), retrieve SOFT/Matrix files, for transcriptomics and expression analysis. |
| gwas-database | Query NHGRI-EBI GWAS Catalog for SNP-trait associations. Search variants by rs ID, disease/trait, gene, retrieve p-values and summary statistics, for genetic epidemiology and polygenic risk scores. |
| hmdb-database | Access Human Metabolome Database (220K+ metabolites). Search by name/ID/structure, retrieve chemical properties, biomarker data, NMR/MS spectra, pathways, for metabolomics and identification. |
| kegg-database | Query KEGG REST API for pathway analysis, gene-pathway mapping, metabolic pathways, drug interactions, ID conversion, compound structures, for systems biology and enrichment analysis. |
| metabolomics-workbench-database | Access NIH Metabolomics Workbench via REST API (4,200+ studies). Query metabolites, RefMet nomenclature, MS/NMR data, m/z searches, study metadata, for metabolomics and biomarker discovery. |
| opentargets-database | Query Open Targets Platform for target-disease associations, drug target discovery, tractability/safety data, genetics/omics evidence, known drugs, for therapeutic target identification. |
| pdb-database | Access RCSB PDB for 3D protein/nucleic acid structures. Search by text/sequence/structure, download coordinates (PDB/mmCIF), retrieve metadata, for structural biology and drug discovery. |
| pubchem-database | Query PubChem via PUG-REST API/PubChemPy (110M+ compounds). Search by name/CID/SMILES, retrieve properties, similarity/substructure searches, bioactivity, for cheminformatics. |
| pubmed-database | Search PubMed biomedical literature. Advanced queries with Boolean/MeSH/field tags, E-utilities API access, batch processing, citation management, for literature review and systematic analysis. |
| reactome-database | Query Reactome REST API for pathway analysis, enrichment, gene-pathway mapping, disease pathways, molecular interactions, expression analysis, for systems biology studies. |
| string-database | Query STRING API for protein-protein interactions (59M proteins, 20B interactions). Network analysis, GO/KEGG enrichment, interaction discovery, 5000+ species, for systems biology. |
| uniprot-database | Query UniProt protein database REST API. Search proteins by name/gene/accession, retrieve sequences (FASTA), annotations, ID mapping, Swiss-Prot/TrEMBL, GO terms, for protein analysis. |
| uspto-database | Access USPTO APIs for patent/trademark searches, examination history (PEDS), assignments, citations, office actions, TSDR, for IP analysis and prior art searches. |
| zinc-database | Access ZINC (230M+ purchasable compounds). Search by ZINC ID/SMILES, similarity searches, 3D-ready structures for docking, analog discovery, for virtual screening and drug discovery. |

## Scientific Helpers

| Name | Description |
|------|-------------|
| context-initialization | Always Auto-invoked skill that creates/updates workspace AGENT.md to instruct the agent to always search for existing skills before attempting any scientific task. |
| get-available-resources | This skill should be used at the start of any computationally intensive scientific task to detect and report available system resources (CPU cores, GPUs, memory, disk space). It creates a JSON file with resource information and strategic recommendations that inform computational approach decisions such as whether to use parallel processing (joblib, multiprocessing), out-of-core computing (Dask, Zarr), GPU acceleration (PyTorch, JAX), or memory-efficient strategies. Use this skill before running analyses, training models, processing large datasets, or any task where resource constraints matter. |

## Scientific Integrations

| Name | Description |
|------|-------------|
| benchling-integration | Benchling R&D platform integration. Access registry (DNA, proteins), inventory, ELN entries, workflows via API, build Benchling Apps, query Data Warehouse, for lab data management automation. |
| dnanexus-integration | DNAnexus cloud genomics platform. Build apps/applets, manage data (upload/download), dxpy Python SDK, run workflows, FASTQ/BAM/VCF, for genomics pipeline development and execution. |
| labarchive-integration | Electronic lab notebook API integration. Access notebooks, manage entries/attachments, backup notebooks, integrate with Protocols.io/Jupyter/REDCap, for programmatic ELN workflows. |
| latchbio-integration | Latch platform for bioinformatics workflows. Build pipelines with Latch SDK, @workflow/@task decorators, deploy serverless workflows, LatchFile/LatchDir, Nextflow/Snakemake integration. |
| omero-integration | Microscopy data management platform. Access images via Python, retrieve datasets, analyze pixels, manage ROIs/annotations, batch processing, for high-content screening and microscopy workflows. |
| opentrons-integration | Lab automation platform for Flex/OT-2 robots. Write Protocol API v2 protocols, liquid handling, hardware modules (heater-shaker, thermocycler), labware management, for automated pipetting workflows. |

## Scientific Packages

| Name | Description |
|------|-------------|
| anndata | Manipulate AnnData objects for single-cell genomics. Load/save .h5ad files, manage obs/var metadata, layers, embeddings (PCA/UMAP), concatenate datasets, for scRNA-seq workflows. |
| arboreto | Gene regulatory network inference with GRNBoost2/GENIE3 algorithms. Infer TF-target relationships from expression data, scalable with Dask, for scRNA-seq and GRN analysis. |
| astropy | Astronomy toolkit. FITS I/O, celestial coordinate transforms, cosmology calculations, time systems, WCS, units, astronomical tables, for astronomical data analysis and imaging. |
| biomni | Autonomous biomedical AI agent framework for executing complex research tasks across genomics, drug discovery, molecular biology, and clinical analysis. Use this skill when conducting multi-step biomedical research including CRISPR screening design, single-cell RNA-seq analysis, ADMET prediction, GWAS interpretation, rare disease diagnosis, or lab protocol optimization. Leverages LLM reasoning with code execution and integrated biomedical databases. |
| biopython | Toolkit for computational molecular biology. Manipulate sequences, parse biological file formats (FASTA, GenBank, FASTQ, PDB), access NCBI databases, run BLAST, analyze structures, build phylogenetic trees, for bioinformatics workflows. |
| bioservices | Access 40+ bio web services (UniProt, KEGG, ChEMBL, PubChem, BLAST, QuickGO, Reactome). Retrieve data, map IDs, pathway analysis, for multi-database integration workflows. |
| cellxgene-census | Query CZ CELLxGENE Census (61M+ cells). Filter by cell type/tissue/disease, retrieve expression data, integrate with scanpy/PyTorch, for population-scale single-cell analysis. |
| cobrapy | Constraint-based metabolic modeling (COBRA). FBA, FVA, gene knockouts, flux sampling, SBML models, for systems biology and metabolic engineering analysis. |
| dask | Parallel/distributed computing. Scale pandas/NumPy beyond memory, parallel DataFrames/Arrays, multi-file processing, task graphs, for larger-than-RAM datasets and parallel workflows. |
| datamol | Pythonic RDKit wrapper for cheminformatics. SMILES parsing, standardization, descriptors, fingerprints, clustering, 3D conformers, parallel processing, for drug discovery workflows. |
| deepchem | Molecular machine learning toolkit. Property prediction (ADMET, toxicity), GNNs (GCN, MPNN), MoleculeNet benchmarks, pretrained models, featurization, for drug discovery ML. |
| deeptools | NGS analysis toolkit. BAM to bigWig conversion, QC (correlation, PCA, fingerprints), heatmaps/profiles (TSS, peaks), for ChIP-seq, RNA-seq, ATAC-seq visualization. |
| diffdock | Diffusion-based molecular docking. Predict protein-ligand binding poses from PDB/SMILES, confidence scores, virtual screening, for structure-based drug design. Not for affinity prediction. |
| esm | Comprehensive toolkit for protein language models including ESM3 (generative multimodal protein design across sequence, structure, and function) and ESM C (efficient protein embeddings and representations). Use this skill when working with protein sequences, structures, or function prediction; designing novel proteins; generating protein embeddings; performing inverse folding; or conducting protein engineering tasks. Supports both local model usage and cloud-based Forge API for scalable inference. |
| etetoolkit | Phylogenetic tree toolkit (ETE). Tree manipulation (Newick/NHX), evolutionary event detection, orthology/paralogy, NCBI taxonomy, visualization (PDF/SVG), for phylogenomics. |
| flowio | Parse FCS (Flow Cytometry Standard) files v2.0-3.1. Extract events as NumPy arrays, read metadata/channels, convert to CSV/DataFrame, for flow cytometry data preprocessing. |
| gget | Bioinformatics query toolkit. Gene info (Ensembl/UniProt), BLAST, AlphaFold, ARCHS4 expression, Enrichr, OpenTargets, COSMIC, genome downloads, CLI/Python, for rapid database queries. |
| matchms | Mass spectrometry analysis. Process mzML/MGF/MSP, spectral similarity (cosine, modified cosine), metadata harmonization, compound ID, for metabolomics and MS data processing. |
| matplotlib | Foundational plotting library. Create line plots, scatter, bar, histograms, heatmaps, 3D, subplots, export PNG/PDF/SVG, for scientific visualization and publication figures. |
| medchem | Medicinal chemistry filters. Apply drug-likeness rules (Lipinski, Veber), PAINS filters, structural alerts, complexity metrics, for compound prioritization and library filtering. |
| molfeat | Molecular featurization for ML (100+ featurizers). ECFP, MACCS, descriptors, pretrained models (ChemBERTa), convert SMILES to features, for QSAR and molecular ML. |
| paper-2-web | This skill should be used when converting academic papers into promotional and presentation formats including interactive websites (Paper2Web), presentation videos (Paper2Video), and conference posters (Paper2Poster). Use this skill for tasks involving paper dissemination, conference preparation, creating explorable academic homepages, generating video abstracts, or producing print-ready posters from LaTeX or PDF sources. |
| polars | Fast DataFrame library (Apache Arrow). Select, filter, group_by, joins, lazy evaluation, CSV/Parquet I/O, expression API, for high-performance data analysis workflows. |
| pydeseq2 | Differential gene expression analysis (Python DESeq2). Identify DE genes from bulk RNA-seq counts, Wald tests, FDR correction, volcano/MA plots, for RNA-seq analysis. |
| pydicom | Python library for working with DICOM (Digital Imaging and Communications in Medicine) files. Use this skill when reading, writing, or modifying medical imaging data in DICOM format, extracting pixel data from medical images (CT, MRI, X-ray, ultrasound), anonymizing DICOM files, working with DICOM metadata and tags, converting DICOM images to other formats, handling compressed DICOM data, or processing medical imaging datasets. Applies to tasks involving medical image analysis, PACS systems, radiology workflows, and healthcare imaging applications. |
| pyhealth | Comprehensive healthcare AI toolkit for developing, testing, and deploying machine learning models with clinical data. This skill should be used when working with electronic health records (EHR), clinical prediction tasks (mortality, readmission, drug recommendation), medical coding systems (ICD, NDC, ATC), physiological signals (EEG, ECG), healthcare datasets (MIMIC-III/IV, eICU, OMOP), or implementing deep learning models for healthcare applications (RETAIN, SafeDrug, Transformer, GNN). |
| pymatgen | Materials science toolkit. Crystal structures (CIF, POSCAR), phase diagrams, band structure, DOS, Materials Project integration, format conversion, for computational materials science. |
| pymc-bayesian-modeling | Bayesian modeling with PyMC. Build hierarchical models, MCMC (NUTS), variational inference, LOO/WAIC comparison, posterior checks, for probabilistic programming and inference. |
| pymoo | Multi-objective optimization framework. NSGA-II, NSGA-III, MOEA/D, Pareto fronts, constraint handling, benchmarks (ZDT, DTLZ), for engineering design and optimization problems. |
| pyopenms | Mass spectrometry toolkit (OpenMS Python). Process mzML/mzXML, peak picking, feature detection, peptide ID, proteomics/metabolomics workflows, for LC-MS/MS analysis. |
| pysam | Genomic file toolkit. Read/write SAM/BAM/CRAM alignments, VCF/BCF variants, FASTA/FASTQ sequences, extract regions, calculate coverage, for NGS data processing pipelines. |
| pytdc | Therapeutics Data Commons. AI-ready drug discovery datasets (ADME, toxicity, DTI), benchmarks, scaffold splits, molecular oracles, for therapeutic ML and pharmacological prediction. |
| pytorch-lightning | Deep learning framework (PyTorch Lightning). Organize PyTorch code into LightningModules, configure Trainers for multi-GPU/TPU, implement data pipelines, callbacks, logging (W&B, TensorBoard), distributed training (DDP, FSDP, DeepSpeed), for scalable neural network training. |
| rdkit | Cheminformatics toolkit for molecules. Parse SMILES/SDF, calculate descriptors (MW, LogP, TPSA), fingerprints, substructure search, 2D/3D generation, similarity, for drug discovery. |
| reportlab | PDF generation toolkit. Create invoices, reports, certificates, forms, charts, tables, barcodes, QR codes, Canvas/Platypus APIs, for professional document automation. |
| scanpy | Single-cell RNA-seq analysis. Load .h5ad/10X data, QC, normalization, PCA/UMAP/t-SNE, Leiden clustering, marker genes, cell type annotation, trajectory, for scRNA-seq analysis. |
| scikit-bio | Biological data toolkit. Sequence analysis, alignments, phylogenetic trees, diversity metrics (alpha/beta, UniFrac), ordination (PCoA), PERMANOVA, FASTA/Newick I/O, for microbiome analysis. |
| scikit-learn | ML toolkit. Classification, regression, clustering, PCA, preprocessing, pipelines, GridSearch, cross-validation, RandomForest, SVM, for general machine learning workflows. |
| scikit-survival | Comprehensive toolkit for survival analysis and time-to-event modeling in Python using scikit-survival. Use this skill when working with censored survival data, performing time-to-event analysis, fitting Cox models, Random Survival Forests, Gradient Boosting models, or Survival SVMs, evaluating survival predictions with concordance index or Brier score, handling competing risks, or implementing any survival analysis workflow with the scikit-survival library. |
| scvi-tools | This skill should be used when working with single-cell omics data analysis using scvi-tools, including scRNA-seq, scATAC-seq, CITE-seq, spatial transcriptomics, and other single-cell modalities. Use this skill for probabilistic modeling, batch correction, dimensionality reduction, differential expression, cell type annotation, multimodal integration, and spatial analysis tasks. |
| seaborn | Statistical visualization. Scatter, box, violin, heatmaps, pair plots, regression, correlation matrices, KDE, faceted plots, for exploratory analysis and publication figures. |
| shap | Model interpretability and explainability using SHAP (SHapley Additive exPlanations). Use this skill when explaining machine learning model predictions, computing feature importance, generating SHAP plots (waterfall, beeswarm, bar, scatter, force, heatmap), debugging models, analyzing model bias or fairness, comparing models, or implementing explainable AI. Works with tree-based models (XGBoost, LightGBM, Random Forest), deep learning (TensorFlow, PyTorch), linear models, and any black-box model. |
| simpy | Process-based discrete-event simulation framework in Python. Use this skill when building simulations of systems with processes, queues, resources, and time-based events such as manufacturing systems, service operations, network traffic, logistics, or any system where entities interact with shared resources over time. |
| statsmodels | Statistical modeling toolkit. OLS, GLM, logistic, ARIMA, time series, hypothesis tests, diagnostics, AIC/BIC, for rigorous statistical inference and econometric analysis. |
| tooluniverse | Use this skill when working with scientific research tools and workflows across bioinformatics, cheminformatics, genomics, structural biology, proteomics, and drug discovery. This skill provides access to 600+ scientific tools including machine learning models, datasets, APIs, and analysis packages. Use when searching for scientific tools, executing computational biology workflows, composing multi-step research pipelines, accessing databases like OpenTargets/PubChem/UniProt/PDB/ChEMBL, performing tool discovery for research tasks, or integrating scientific computational resources into LLM workflows. |
| torch-geometric | Graph Neural Networks (PyG). Node/graph classification, link prediction, GCN, GAT, GraphSAGE, heterogeneous graphs, molecular property prediction, for geometric deep learning. |
| torchdrug | Graph-based drug discovery toolkit. Molecular property prediction (ADMET), protein modeling, knowledge graph reasoning, molecular generation, retrosynthesis, GNNs (GIN, GAT, SchNet), 40+ datasets, for PyTorch-based ML on molecules, proteins, and biomedical graphs. |
| transformers | Work with state-of-the-art machine learning models for NLP, computer vision, audio, and multimodal tasks using HuggingFace Transformers. This skill should be used when fine-tuning pre-trained models, performing inference with pipelines, generating text, training sequence models, or working with BERT, GPT, T5, ViT, and other transformer architectures. Covers model loading, tokenization, training with Trainer API, text generation strategies, and task-specific patterns for classification, NER, QA, summarization, translation, and image tasks. (plugin:scientific-packages@claude-scientific-skills) |
| umap-learn | UMAP dimensionality reduction. Fast nonlinear manifold learning for 2D/3D visualization, clustering preprocessing (HDBSCAN), supervised/parametric UMAP, for high-dimensional data. |
| zarr-python | Chunked N-D arrays for cloud storage. Compressed arrays, parallel I/O, S3/GCS integration, NumPy/Dask/Xarray compatible, for large-scale scientific computing pipelines. |

## Scientific Thinking

| Name | Description |
|------|-------------|
| docx | Document toolkit (.docx). Create/edit documents, tracked changes, comments, formatting preservation, text extraction, for professional document processing. |
| exploratory-data-analysis | EDA toolkit. Analyze CSV/Excel/JSON/Parquet files, statistical summaries, distributions, correlations, outliers, missing data, visualizations, markdown reports, for data profiling and insights. |
| hypothesis-generation | Generate testable hypotheses. Formulate from observations, design experiments, explore competing explanations, develop predictions, propose mechanisms, for scientific inquiry across domains. |
| pdf | PDF manipulation toolkit. Extract text/tables, create PDFs, merge/split, fill forms, for programmatic document processing and analysis. |
| peer-review | Systematic peer review toolkit. Evaluate methodology, statistics, design, reproducibility, ethics, figure integrity, reporting standards, for manuscript and grant review across disciplines. |
| pptx | Presentation toolkit (.pptx). Create/edit slides, layouts, content, speaker notes, comments, for programmatic presentation creation and modification. |
| scientific-brainstorming | Research ideation partner. Generate hypotheses, explore interdisciplinary connections, challenge assumptions, develop methodologies, identify research gaps, for creative scientific problem-solving. |
| scientific-critical-thinking | Evaluate research rigor. Assess methodology, experimental design, statistical validity, biases, confounding, evidence quality (GRADE, Cochrane ROB), for critical analysis of scientific claims. |
| scientific-visualization | Create publication figures with matplotlib/seaborn/plotly. Multi-panel layouts, error bars, significance markers, colorblind-safe, export PDF/EPS/TIFF, for journal-ready scientific plots. |
| scientific-writing | Write scientific manuscripts. IMRAD structure, citations (APA/AMA/Vancouver), figures/tables, reporting guidelines (CONSORT/STROBE/PRISMA), abstracts, for research papers and journal submissions. |
| statistical-analysis | Statistical analysis toolkit. Hypothesis tests (t-test, ANOVA, chi-square), regression, correlation, Bayesian stats, power analysis, assumption checks, APA reporting, for academic research. |
| xlsx | Spreadsheet toolkit (.xlsx/.csv). Create/edit with formulas/formatting, analyze data, visualization, recalculate formulas, for spreadsheet processing and analysis. |

