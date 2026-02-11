# Trends in Computational Metabolomics (2021–2025)

[![Awesome](https://cdn.jsdelivr.net/gh/sindresorhus/awesome@d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)
A curated list of computational tools, databases, and software for metabolomics research published between 2021-2025.

## Citation

> Domingo-Fernández, D., Healy, D., Kind, T., Allen, A., Colluru, V., and Misra, B. (2025). Trends in computational metabolomics in the last five years (2021–2025). chemRxiv. <https://doi.org/10.26434/chemrxiv-2025-f6lvq>

## Methodology

We compiled this comprehensive list of computational tools and databases through a systematic approach. We began with manual curation of published metabolomics resources starting in 2021, building upon our prior review efforts ([Misra, 2021](https://doi.org/10.1007/s11306-021-01796-1); [Misra and van der Hooft, 2015](https://doi.org/10.1002/elps.201500417)). To ensure complete coverage and minimize gaps, we subsequently employed three complementary strategies:

1. **Computational screening**: Used AI-assisted tools and Google Scholar to identify emerging metabolomics software and databases across the 2021–2025 period
2. **Literature review**: Systematically examined citations and references in recent comprehensive metabolomics reviews for each functional category
3. **Field experts**: Reviewed our compilation with domain experts across different metabolomics subdisciplines

## Collaborating

Please make a PR editing the [list](tools_list.tsv) to add new tools.

## 📚 Table of Contents

- [🏷️ Annotation](#️-annotation)
- [🧪 Benchmark/Dataset](#-benchmarkdataset)
- [🧬 Biosynthetic Gene Clusters](#-biosynthetic-gene-clusters)
- [💧 CE-MS](#-ce-ms)
- [⚡ DIMS](#-dims)
- [🗃️ Database](#️-database)
- [💊 Drug Discovery](#-drug-discovery)
- [🌱 Exposomics](#-exposomics)
- [🌀 FT-ICR MS](#-ft-icr-ms)
- [🗂️ Formats](#️-formats)
- [🔥 GC-MS](#-gc-ms)
- [🔋 Holistic/Standalone Tools](#-holisticstandalone-tools)
- [🌈 IR](#-ir)
- [🖼️ Imaging MS](#️-imaging-ms)
- [🌬️ Ion Mobility MS](#️-ion-mobility-ms)
- [🧮 Isotopic](#-isotopic)
- [🔭 Large Scale](#-large-scale)
- [🧬 Lipidomics](#-lipidomics)
- [🕸️ Metabolic Networks](#️-metabolic-networks)
- [🗃️ Metadata](#️-metadata)
- [🔀 Multifunctional](#-multifunctional)
- [🔗 Multiomics](#-multiomics)
- [🔬 NMR](#-nmr)
- [🌱 Organism Specific](#-organism-specific)
- [🛤️ Pathway, Enrichment and Ontology Tools](#️-pathway-enrichment-and-ontology-tools)
- [🔎 Patterns](#-patterns)
- [⚙️ Pre-processing](#️-pre-processing)
- [✅ Quality Control](#-quality-control)
- [⏱️ RT (Retention Time)](#️-rt-retention-time)
- [🦠 Single Cell Metabolomics](#-single-cell-metabolomics)
- [🗺️ Spatial Metabolomics](#️-spatial-metabolomics)
- [🧩 Specialized](#-specialized)
- [📕 Spectral Library](#-spectral-library)
- [📈 Statistical](#-statistical)
- [🎯 Targeted](#-targeted)
- [📊 Visualization](#-visualization)

---

## 🏷️ Annotation

### 🧪 Class | Property Prediction

- [CANOPUS](https://bio.informatik.uni-jena.de/software/canopus/) - Class annotation tool
- [Mass Spectrum Transformer](https://github.com/chensaian/TransG-Net) - Multi-modal fusion of molecular graphs and mass spectra for property prediction
- [MWFormer](https://github.com/zhanghailiangcsu/MWFormer) - Molecular weight estimation from EI-MS spectra
- [Spec2Class](https://huggingface.co/VickiPol/binary_models) - Spectrum to class prediction

### 🔍 De Novo Generation

- [DiffMS](https://github.com/coleygroup/DiffMS) - Diffusion-based de novo structure generation from MS/MS spectra
- [MASSISTANT](https://github.com/BartaLazar/MASSISTANT/) - De novo molecular structure prediction from EI-MS spectra via SELFIES encoding
- [MS-BART](https://github.com/OpenDFM/MS-BART) - Unified encoder-decoder model for de novo structure elucidation from mass spectra
- [MS2Mol](https://doi.org/10.26434/chemrxiv-2023-vsmpx-v3) - MS/MS to molecule prediction
- [MSGo](http://github.com/aaronma2020/MSGO) - Generates molecular structures directly from mass spectra for exposomics
- [MSNovelist](https://github.com/meowcat/MSNovelist) - De novo structure elucidation
- [Mass2SMILES](https://github.com/volvox292/mass2smiles) - Mass to SMILES conversion
- [MassGenie](https://github.com/neilswainston/FragGenie) - De novo molecular structure prediction from mass spectra using transformers
- [OMG](https://github.com/HassounLab/OMGG) - Optimal molecule generation
- [SEISMiQ](https://github.com/Boehringer-Ingelheim/seismiq) - Structure elucidation
- [Spec2Mol](https://github.com/KavrakiLab/Spec2Mol) - Spectrum to molecule
- [TeFT](https://github.com/thumingo/TeFT) - Transformer for fragmentation

### 🧠 Learned Spectrum Representations

- [ChemEmbed](https://github.com/massspecdl/ChemEmbed) - Metabolite identification using enhanced MS/MS data and multidimensional molecular embeddings
- [CLERMS](https://github.com/HaldamirS/CLERMS) - Contrastive learning-based embeddings for tandem mass spectra
- [DeepMASS](https://github.com/hcji/DeepMASS2_Data_Processing) - Compound annotation via semantic similarity of mass spectral language
- [DreaMS](https://github.com/pluskal-lab/DreaMS) - Self-supervised transformer pre-trained on millions of unannotated MS/MS spectra
- [LSM1-MS2](https://github.com/matterworksbio/LSM1-MS2) - Self-supervised foundation model for tandem mass spectrometry
- [MS2DeepScore](https://github.com/matchms/ms2deepscore) - Deep learning similarity scoring
- [MS2DeepScore 2.0](https://github.com/matchms/MS2DeepScore) - Updated deep learning similarity
- [MSBERT](https://github.com/zhanghailiangcsu/MSBERT) - BERT for mass spectrometry
- [NaFM](https://github.com/TomAIDD/NaFM-Official) - Pre-trained foundation model for small-molecule natural products
- [Spec2Vec](https://github.com/iomega/spec2vec) - Word2Vec for spectra
- [SpecEmbedding](https://huggingface.co/spaces/xp113280/SpecEmbedding) - Spectrum embedding

### 🏎️ Molecular Formula Prediction

- [CRB-FCC](https://github.com/iconSS/FCC/releases/tag/v1) - Fragmental chain characterization
- [FIDDLE](https://github.com/JosieHong/FIDDLE) - Formula identification
- [FSA](https://github.com/kslynn128171/FSA) - Formula subset analysis
- [IDSL.UFA](https://cran.r-project.org/package=IDSL.UFA) - United formula annotation
- [MIST-CF](https://github.com/samgoldman97/mist-cf) - Chemical formula prediction
- [RASSP](https://github.com/thejonaslab/rassp-public) - Rapid approximate subset-based EI-MS spectra prediction
- [BUDDY](https://github.com/HuanLab/BUDDY) - Formula prediction and annotation

### 📚 Molecular Library Retrieval

- [CMSSP](https://huggingface.co/OliXio/CMSSP) - Cross-modal spectrum-structure prediction
- [COSMIC](https://bio.informatik.uni-jena.de/cosmic/) - Confidence scoring for structure annotation
- [CSU-MS2](https://github.com/tingxiecsu/CSU-MS2) - MS/MS library search
- [IDSL_MINT](https://github.com/idslme/IDSL_MINT) - Metabolite identification
- [JESTR](https://github.com/HassounLab/JESTR1/) - Joint embedding for structure retrieval
- [LC-MS2Struct](https://github.com/aalto-ics-kepaco/msms_rt_ssvm) - LC-MS/MS structure annotation
- [MIST](https://github.com/samgoldman97/mist) - Mass spectrum identification
- [MVP](https://github.com/HassounLab/MVP) - Metabolite virtual profiling
- [VInSMoC](https://github.com/mohimanilab/VInSMoC) - Variable interpretation of spectrum-molecule couples

### 🕸️ Molecular Networking

- [3D-MPEA](https://github.com/ZibianFan/3D-MPEA) - Graph attention model for multidimensional omics annotation
- [BAM](https://github.com/HassounLab/BAM) - Biotransformation-based annotation method
- [ChemWalker](https://github.com/computational-chemical-biology/ChemWalker) - Annotation propagation via random walks
- [ConCISE](http://github.com/zquinlan/concise) - Consensus annotation propagation fusing molecular networking and in silico predictions
- [E-SGMN](https://sourceforge.net/projects/e-sgmn/) - Enhanced structure-guided molecular networking
- [IIMN](https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/) - Ion identity molecular networking
- [IMN4NPD](https://github.com/mwang87/NP-Classifier) - Integrated molecular networking for natural products
- [MCN](https://github.com/Alexander0/molecular_communities) - Optimizes per-metabolite connectivity for molecular networking
- [MMSA](https://multiplealignment.gnps2.org/setscreation) - Interactive visualization of fragmentation patterns across molecular families
- [ModiFinder](https://github.com/Wang-Bioinformatics-Lab/ModiFinder_base) - Modification site localization
- [MolNotator](https://github.com/ZzakB/MolNotator) - Reduces molecular network redundancy by grouping ion species into neutral molecules
- [MS-Net](https://zenodo.org/records/17669288) - Multi-similarity network-based metabolite annotation
- [SGMNS](https://github.com/DLUT-datas/SGMNS) - Structure-guided molecular network strategy
- [SIMILE](https://github.com/biorack/simile) - Spectral alignment with statistical significance

### Others

- [AnnoMe](https://github.com/chrboku/AnnoMe) - MS/MS spectra classification
- [AnnoSM](https://github.com/YangHuaLab/AnnoSM) - Substituent mode annotation
- [BioTransformer4.0](https://github.com/Wishartlab-openscience/Biotransformer) - Biotransformation prediction
- [IDSL.CSA](https://cran.r-project.org/package=IDSL.CSA) - Composite spectra analysis
- [Inventa](https://luigiquiros.github.io/inventa/) - Structural novelty discovery
- [ipaPy2](https://github.com/francescodc87/ipaPy2) - Integrated probabilistic annotation
- [MADGEN](https://github.com/HassounLab/MADGEN) - Scaffold-based de novo molecular generation from MS/MS spectra
- [MCheM](https://github.com/sirius-ms/sirius) - Multiplexed post-column derivatization for improved metabolite annotation
- [MetaboAnnotatoR](https://github.com/gggraca/MetaboAnnotatoR) - Automated annotation of all-ion fragmentation LC-MS data
- [MMST](https://github.com/mpriessner/MultiModalSpectralTransformer) - Predicts chemical structures from multimodal spectroscopic data
- [MS2DECIDE](https://github.com/MejriY/MS2DECIDE) - Decision theory for annotation
- [mWISE](https://github.com/b2slab/mWISE) - Context-based annotation of LC-MS features through diffusion in graphs
- [NetID](https://github.com/LiChenPU/NetID/releases/tag/v1.0) - Network-based identification
- [OrbiFragsNets](https://github.com/EdwinChingate/OrbiFragsNets) - Orbitrap fragmentation networks

### 🔍 Spectral Similarity Retrieval

- [BLINK](https://github.com/biorack/blink) - Fast spectral similarity
- [compareMS2 2.0](https://github.com/524D/compareMS2) - MS/MS comparison
- [falcon](https://github.com/bittremieux-lab/falcon) - Fast spectral clustering
- [FastEI](https://github.com/Qiong-Yang/FastEI) - Fast EI spectral matching
- [Flash entropy search](https://github.com/YuanyueLi/EntropySearch) - Entropy-based search
- [MASST+](https://github.com/mohimanilab/MASSTplus) - Enhanced mass spectrum search
- [metID](https://jaspershen.github.io/metID/) - Metabolite identification
- [mineMS2](https://github.com/odisce/mineMS2) - MS/MS mining
- [MolDiscovery](https://github.com/mohimanilab/molDiscovery) - Molecule discovery
- [MS/MS spectral entropy similarity](https://github.com/YuanyueLi/MSEntropy) - Entropy-based similarity
- [Reverse Spectral Search](https://github.com/Philipbear/reverse_search) - Reverse search algorithm
- [SimMS](https://github.com/PangeAI/simms) - Similarity for mass spectra
- [TransExION](https://github.com/banhdzui/TransExION) - Transformer for ion extraction
- [Weighted spectral similarity for library search](https://github.com/enveda/weighting-spectral-similarity) - Weighted similarity

### ⚗️ Spectrum Prediction

- [3DMolMS](https://github.com/JosieHong/3DMolMS) - 3D molecular MS prediction
- [CFM-ID 4.0](https://cfmid.wishartlab.com/) - Competitive fragmentation modeling
- [CIDMD](https://github.com/jesilee/) - Collision-induced dissociation via molecular dynamics
- [DeepCDM](https://github.com/ADNLab-SCU/DeepCDMs) - Deep learning prediction of MS/MS spectra for chemically derived molecules
- [ESP](https://github.com/HassounLab/ESP) - Ensemble spectral prediction model for metabolite annotation
- [FIORA](https://github.com/BAMeScience/fiora) - Fragmentation prediction
- [fragnnet](https://github.com/FraGNNet/fragnnet) - Fragmentation neural network
- [GrAFF-MS](https://www.envedabio.com/posts/a-new-scalable-ml-model-for-accurate-prediction-of-small-molecule-mass-spectra) - Graph attention for fragmentation
- [HDSE-MS](https://github.com/lzjforyou/HDSE-MS) - MS/MS spectrum prediction via hierarchical distance structural encoding
- [ICEBERG](https://github.com/coleygroup/ms-pred) - Geometric deep learning model for CID mass spectrum simulation
- [MARASON](https://github.com/coleygroup/ms-pred) - Mass spectrum prediction
- [MassFormer](https://github.com/Roestlab/massformer/) - Transformer for mass spectra
- [MassKG](https://xomics.com.cn/masskg) - Knowledge-based fragmentation and deep learning for natural product MS annotation
- [MS2Compound](https://github.com/beherasan/MS2Compound) - Generates custom predicted spectra databases for metabolite identification
- [ms-pred](https://github.com/coleygroup/ms-pred) - Prefix-tree decoding for predicting mass spectra from molecules
- [NPS-MS](https://nps-ms.ca/) - Deep learning MS/MS prediction for novel psychoactive substance identification
- [PPGB-MS2](https://github.com/zhengfj1994/PPGB_MS2) - MS/MS prediction using precursor-product ion pair graph bag embeddings
- [QC-GN2oMS2](https://github.com/PNNL-m-q/qcgnoms) - Quantum chemistry MS prediction
- [SingleFrag](https://github.com/MaribelPR/SingleFrag) - Single fragmentation prediction

### 🧰 Workflow Tools

- [ENTAiLS Toolkit](https://whitehead-heather.github.io/ENTAiLSToolkit/) - Annotation toolkit
- [MetaboCoreUtils, MetaboAnnotation and CompoundDb](https://github.com/jorainer/MetaboAnnotationTutorials) - R packages for annotation
- [MetaFetcheR](https://github.com/komorowskilab/MetaFetcheR/) - Metabolite fetching
- [MS2Query](https://github.com/iomega/ms2query) - Query-based annotation
- [mscompiler](https://github.com/QizhiSu/mspcompiler) - R package for compiling and searching MS/MS spectral libraries
- [mspepsearchr](https://github.com/AndreySamokhin/mssearchr) - R package for high-throughput mass spectral library searching
- [PS2MS](https://github.com/jhhung/PS2MS) - Deep learning-based identification of new psychoactive substances from mass spectra
- [RapidMass](https://github.com/Katherine00689/RapidMass) - Automated species authentication for high-throughput mass spectrometry
- [Scannotation](https://github.com/scannotation/Scannotation_software) - Scan annotation
- [ShinyMetID](https://github.com/jjs3098/ShinyMetID) - Shiny app for metabolite ID
- [SiMD](https://si-simd.com/) - Curated metabolite reference spectral library for LC-MS identification

---

## 🧪 Benchmark/Dataset

- [MassSpecGym](https://github.com/pluskal-lab/MassSpecGym) - Benchmark for MS/MS molecule discovery and identification tasks
- [MetaBench](https://github.com/metabench/metabench) - Benchmark for evaluating LLM capabilities in metabolomics

---

## 🧬 Biosynthetic Gene Clusters

- [antiSMASH 8.0](https://antismash.secondarymetabolites.org/) - Biosynthetic gene cluster detection
- [iPRESTO](https://pypi.org/project/ipresto/) - Prediction of secondary metabolites
- [MariClus](https://www.mariclus.com) - Marine cluster analysis
- [MIBiG 3.0](https://mibig.secondarymetabolites.org/) - Minimum information about BGCs
- [NPClassScore](https://github.com/NPLinker/nplinker) - Natural product classification scoring
- [NPLinker](https://github.com/NPLinker/nplinker) - Natural products linking
- [NPOmix](https://github.com/tiagolbiotech/NPOmix_python) - Natural products omics

---

## 💧 CE-MS

- [MobilityTransformR](https://github.com/LiesaSalzer/MobilityTransformR) - Mobility transformation
- [PeakMeister](https://github.com/PBM-Group/PeakMeister) - Peak detection for CE-MS

---

## ⚡ DIMS

- [EASY-FIA](https://github.com/AMrbt20/EASY-FIA/) - Flow injection analysis
- [rIDIMS](https://github.com/BioinovarLab/rIDIMS) - R package for DIMS
- [Tidy-Direct-to-MS](https://github.com/griquelme/tidyms) - Tidy workflow for direct MS

---

## 🗃️ Database

- [AMDB](https://amdb.online/) - Quantitative animal metabolite database with statistical analysis tools
- [BinDiscover database](https://bindiscover.metabolomics.us/) - Standardized GC-MS metabolomics data from 156,000+ samples across 2,000 studies
- [CCDB](https://github.com/idslme/CCDB/tree/main) - Chemical correlation database cataloguing inter-chemical correlations in metabolomics
- [COCONUT](https://coconut.naturalproducts.net/) - Natural products database
- [DNA adduct database](https://gitlab.nexs-metabolomics/projects/dna_adductomics_database) - LC-HRMS database and bioinformatics tool for DNA adductomics screening
- [EnzyMine](http://www.rxnfinder.org/enzymine/) - Knowledge base linking enzymatic reaction features with enzyme sequences
- [foodMASST](https://masst.ucsd.edu/foodmasst) - MS/MS search tool against a food and beverage metabolomics reference database
- [FragHub](https://github.com/eMetaboHUB/FragHub) - Standardized, merged open mass spectral libraries for metabolite annotation
- [GNPS Dashboard](https://dashboard.gnps2.org/) - GNPS data visualization
- [HMDB 5.0](https://hmdb.ca) - Human metabolome database
- [HREI-MSDB](https://github.com/mtshn/gchrmsexplain) - High-resolution EI-MS database
- [Human Hair Atlas](https://metabolomics.cloud/hair/) - Hair metabolome atlas
- [LEAFBot](https://ccms-ucsd.github.io/GNPSDocumentation/gnpslibraries/) - MS/MS spectral library of 300+ botanical natural product standards
- [MarkerDB 2.0](https://markerdb.ca/) - Biomarker database
- [MassBase](http://webs2.kazusa.or.jp/massbase/) - Multi-platform MS depository of metabolites detected across organisms
- [MassSpecBlocks](https://github.com/privrja/MassSpecBlocks) - Database of nonribosomal peptide and polyketide building blocks for MS
- [MCID database](http://www.mycompoundid.org/) - LC-MS metabolite identification using mass search with ppm-level accuracy
- [MedMeta](https://cbcb.cdutcm.edu.cn/MedMeta) - Integrative database connecting medicinal secondary metabolites and biosynthetic pathways
- [MetaboLights](https://www.ebi.ac.uk/metabolights) - Metabolomics data repository
- [Metabolome atlas of the aging mouse brain](https://mouse.atlas.metabolomics.us/) - Mouse brain metabolome
- [MetalinksDB](https://metalinks.omnipathdb.org/) - Prior knowledge resource for metabolite-mediated intercellular signaling
- [MetaNetX 2025 update](https://www.metanetx.org/) - Reconciliation and integration of metabolic data across biochemical databases
- [MetHoS](https://methos.cebitec.uni-bielefeld.de/welcomeview) - Automated web platform for large-scale multi-experiment metabolomics processing
- [METLIN-CCS](https://metlin.scripps.edu/) - METLIN with CCS values
- [microbeMASST](https://github.com/robinschmid/microbe_masst) - MS/MS search tool linking spectra to microbial producers taxonomically
- [MiMeDB](https://mimedb.org/) - Microbiome metabolome database
- [MiMeDB 2.0](https://mimedb.org/) - Updated microbiome metabolome database
- [MSCAT](https://mscat.metabolomicsworkbench.org/) - Catalog and guidance tool for metabolomics software selection
- [MSnLib](https://github.com/corinnabrungs/msn_tree_library) - Open multi-stage fragmentation (MSn) spectral library for ML
- [Natural Products Atlas 2.0](https://www.npatlas.org/) - Natural products database
- [NGlycDB](https://metaspace2020.org/project/velickovic_Nglyc_2021) - N-glycan annotation and visualization for mass spectrometry imaging
- [NIST23](https://chemdata.nist.gov/dokuwiki/doku.php?id=chemdata:nist23-msms) - NIST mass spectral library
- [NMRlipids Databank](https://nmrlipids.github.io/) - Quality-evaluated MD simulation database for lipid membranes
- [NPMine](https://github.com/computational-chemical-biology/npmine) - Text-mining tool extracting natural product structures from literature
- [PharmMet DB](https://github.com/ClinicalBiomarkersLaborabory/PharmMet) - In vitro human liver S9-generated drug metabolite reference library
- [plantMASST](https://github.com/helenamrusso/plantmasst) - MS/MS search tool for plant metabolites across 2,793 species
- [Pyrfume](https://pyrfume.org/) - Unified framework aggregating stimulus-linked olfactory datasets
- [RepoRT](https://github.com/michaelwitting/RepoRT) - Community-curated retention time repository
- [Spectraverse](https://github.com/skinniderlab/spectraverse-analysis) - Preprocessed, quality-filtered, merged public MS/MS spectral library
- [TOMATOMET](http://metabolites.in/tomato‐fruits/) - Tomato metabolome

---

## 💊 Drug Discovery

- [BitBIRCH-Lean](https://github.com/mqcomplab/bblean) - Memory-efficient clustering of billion-scale drug-like molecular libraries
- [ChromaQuant](https://hplcfdu.shinyapps.io/ChromaQuant/) - Identifies drug-perturbed protein complexes via co-fractionation MS
- [DMetFinder](https://github.com/KeWHs/dmetdata) - Drug metabolite finder
- [Limelight](https://limelight-ms.org) - Web-based visualization and analysis of drug-protein adduct MS data
- [SynFrag](https://synfrag.simm.ac.cn/) - Evaluates synthetic accessibility of AI-designed drug molecules

---

## 🌱 Exposomics

- [CMDN](https://github.com/LeaveMeNotTonight/CMDN) - Untargeted metabolomics network for high-throughput xenobiotic metabolite annotation
- [EISA-EXPOSOME](https://github.com/Lab-XUE/EISA-EXPOSOME) - Enhanced in-source fragmentation annotation for sensitive chemical identification
- [ExposomeX](http://www.exposomex.cn/) - Integrated platform for standardized exposome-wide data analysis
- [FeatureHunter](https://msomics.abrc.sinica.edu.tw/FeatureHunter/) - Software for nucleic acid adduct detection from HR-MS/MS data
- [FluoroMatch IM](https://fluoromatch.com) - PFAS identification integrating ion mobility with nontargeted MS
- [HalogenFinder](https://halogenfinder.com) - ML-based workflow for chlorine/bromine disinfection byproduct identification
- [HExpMetDB](https://github.com/FangLabNTU/HExpMetDB) - Human exposome metabolite database
- [MDRB](https://github.com/933ZhangDD/MDRB) - MS-based identification of in vivo traditional Chinese medicine components
- [MetabFlow](https://bddg.hznu.edu.cn/metabflow/) - Database for metabolism of exogenous natural products in vivo and in vitro
- [MSThunder](https://github.com/LQZ0123/MSThunder.git) - Nontargeted analysis framework for rapid unknown organic pollutant identification

---

## 🌀 FT-ICR MS

- [MetaboDirect](https://github.com/Coayala/MetaboDirect) - Automated pipeline for FT-ICR MS organic matter characterization
- [MoleTrans](https://github.com/JibaoLiu/MoleTrans) - Browser-based tool for FT-ICR MS dissolved organic matter formula analysis

---

## 🗂️ Formats

- [Aird](https://github.com/CSi-Studio/AirdPro) - Computation-oriented MS format with high compression and fast decoding
- [AlphaTims](https://github.com/MannLabs/alphatims) - Python package for efficient indexing, access, and visualization of LC-TIMS-MS data
- [CloMet](https://github.com/rmallol/clomet) - Converts raw metabolomics data from public repositories into formats for analysis platforms
- [Dear-OMG](https://github.com/jianweishuai/Dear-OMG) - Unified compressed storage format for proteomics, genomics, and metabolomics metadata
- [mspack](https://github.com/fhanau/mspack) - Lossless and lossy MS data compression exploiting cross-scan redundancy in mzML/mzXML files
- [mwtab Python Library for RESTful Access](https://github.com/MoseleyBioinformaticsLab/mwtab) - Python library for accessing and validating Metabolomics Workbench data via REST API
- [MZA](https://github.com/PNNL-m-q/mza) - HDF5-based storage and access tool for multidimensional MS data including ion mobility
- [mzapy](https://github.com/PNNL-m-q/mzapy) - Python package for extraction and processing of MS data in MZA format
- [mzPeak](https://github.com/mobiusklein/mzpeak_prototyping) - Community data format designed for high-throughput multidimensional MS workflows
- [pyOpenMS-viz](https://github.com/OpenMS/pyopenms_viz) - OpenMS visualization
- [rawrr R](https://github.com/fgcz/rawrr) - Raw file reading in R
- [SpectriPy](https://github.com/rformassspectrometry/SpectriPy) - R package enabling translation between R and Python mass spectrometry data structures
- [TIMSCONVERT](https://github.com/gtluu/timsconvert) - Converts Bruker timsTOF raw data to open mzML and imzML formats

---

## 🔥 GC-MS

- [ADAP-KDB](https://adap.cloud/) - Spectral knowledgebase for tracking and prioritizing unknown GC-MS spectra
- [CRISP](https://github.com/vivekmathema/GCxGC-CRISP) - GCxGC data processing
- [eROI](https://github.com/Gaozhuang-big/eROI) - Equidistant region-of-interest strategy for deconvolving co-eluted peaks in GC-Orbitrap HRMS
- [gc-ims-tools](https://github.com/Charisma-Mannheim/gc-ims-tools) - Open-source Python package for customizable GC-IMS data handling and analysis
- [GcDUO](https://github.com/mariallr/GcDuo) - Open-source R software for GC×GC-MS data processing using PARAFAC/PARAFAC2 deconvolution
- [GCMS-ID](https://gcms-id.ca) - GC-MS identification
- [isoSCAN](https://github.com/jcapelladesto/isoSCAN) - R tool for automated isotopologue quantification in GC-chemical-ionization-MS flux experiments
- [MACE](http://www.oc.tu-bs.de/schulz/html/MACE.html) - Open-access EI-MS spectral database of natural products for chemical ecology
- [MetaboPAC](https://github.com/gtStyLab/MetaboPAC) - Infers absolute concentrations from GC-MS relative abundances using kinetic models
- [mzrtsim](https://github.com/yufree/mzrtsim) - R package for simulating GC/LC full-scan raw data in mzML for software benchmarking
- [RIpred](https://ripred.ca/) - Retention index prediction
- [SERDA](https://github.com/slfan2013/SERDA) - Denoising autoencoder for normalizing systematic signal drift in large-scale GC-MS studies
- [SpecTUS](https://github.com/hejjack/SpecTUS) - Deep learning model for de novo molecular structure annotation from GC-EI-MS spectra
- [uafR](https://github.com/castratton/uafR) - R package automating GC-MS data retrieval, compound identification, and similarity matching
- [UniqPy](https://github.com/DNKonanov/uni_cli) - Thermodynamic model for quantifying short-chain fatty acids from headspace GC-MS

---

## 🔋 Holistic/Standalone Tools

- [FluoroMatch 2.0](http://innovativeomics.com/software/fluoromatch-flow-covers-entire-pfas-workflow/) - PFAS workflow
- [GNPS2](https://gnps2.org) - Global natural products social networking
- [MAVEN2](https://github.com/eugenemel/maven/releases/tag/2.11.8) - Open-source LC-MS/MS metabolomics and lipidomics analysis software
- [MAW](https://github.com/zmahnoor14/MAW) - Metabolomics analysis workflow
- [MassCube](https://github.com/huaxuyu/masscube) - Python framework for LC-MS peak detection, annotation, and quality control
- [MeRgeION](https://github.com/daniellyz/meRgeION2) - Transforms LC-MS/MS data into information-rich spectral databases
- [MetaboReport](https://metaboreport.com/) - Interactive web app for metabolomics analysis and reporting
- [MetaboScape 2025](https://www.bruker.com/en/products-and-solutions/mass-spectrometry/ms-software/metaboscape.html) - Proprietary Bruker software for metabolomics with automated annotation and ion mobility integration
- [Metabox 2.0](http://metsysbio.com/metabox) - Metabolomics toolbox
- [MetEx](https://mo.princeton.edu/MetEx/) - Web app for exploring microbial secondary metabolome datasets
- [PlantMetSuite](https://github.com/YaonanLab/PlantMetSuite) - Web-based platform for plant-specific metabolomics analysis and multi-omics integration
- [Workflow4Metabolomics](https://workflow4metabolomics.org/) - Galaxy workflows

---

## 🌈 IR

- [Chemprop-IR](https://github.com/gfm-collab/chemprop-IR) - ML package for predicting infrared spectra from molecular structures
- [Graphormer-IR](https://github.com/HopkinsLaboratory/Graphormer-IR) - IR prediction with Graphormer
- [Spectro](https://github.com/ChemAI-Lab/spectro) - Multi-modal approach combining IR and NMR data for molecular structure elucidation

---

## 🖼️ Imaging MS

- [13C-SpaceM](https://github.com/Buglakova/13C-SpaceM) - 13C spatial metabolomics
- [Aird-MSI](https://github.com/CSi-Studio/mzmine3/tree/Aird2.0) - Optimized Aird compression for MSI data achieving 70% storage reduction over imzML
- [Cardinal v3](https://github.com/kuwisdelu/Cardinal) - MSI analysis in R
- [DeepION](https://github.com/BioNet-XMU/DeepION) - Deep learning ion image representation for mass spectrometry imaging
- [deepPseudoMSI](https://www.deeppseudomsi.org/) - Deep learning pseudo MSI
- [DIMPLE](https://github.com/dickinsonlab/DIMPLE-code) - Pipeline for clustering metabolite enrichment patterns along developmental gradients in DESI-MSI
- [HT SpaceM](https://github.com/delafior/HT-SpaceM) - High-throughput single-cell metabolomics combining MALDI imaging MS with batch processing
- [i2i](https://github.com/LanekoffLab/i2i) - Processes mass spectrometry imaging data into ion images
- [imzML Writer](https://github.com/VIU-Metabolomics/imzML_Writer) - Python tool converting continuously acquired MSI data into pixel-aligned imzML files
- [LipidQMap](https://github.com/swinnenteam/LipidQMap) - Open-source platform for omics-scale quantitative lipidomics in MSI
- [M2aia](https://github.com/m2aia/m2aia) - Interactive tool for 3D reconstruction, visualization, and analysis of MSI datasets
- [macroMS](https://github.com/macroMS/macroMS_web) - Web-based suite for image-guided MS of macroscopic samples like microbial colonies
- [mass2adduct](https://github.com/kbseah/mass2adduct) - Annotates adducts and background ions in MSI data
- [massNet](https://github.com/wabdelmoula/massNet) - Deep learning model for classifying MSI data directly from raw spectra
- [MassVision](https://SlicerMassVision.readthedocs.io) - MSI platform integrating visualization, segmentation, AI training, and statistical analysis
- [Met-ID](https://github.com/pbjarterot/Met-ID) - GUI software for metabolite identification from MALDI-MSI using MS1/MS2
- [METASPACE-ML](https://github.com/metaspace2020/metaspace) - ML-based metabolite annotation for MSI with improved low-intensity detection
- [MetaVision3D](https://metavision3d.rc.ufl.edu) - Computer-vision pipeline transforming serial 2D MALDI-MSI into 3D spatial metabolome models
- [MSI-Explorer](https://github.com/MMV-Lab/MSI-Explorer) - Napari plugin for interactive MSI data processing, spectral analysis, and annotation
- [MSIannotator](https://github.com/Yingzhu96/MSIannotator) - Organ-specific, database-driven automated metabolite annotation for MSI
- [MSIGen](https://github.com/LabLaskin/MSIGen) - Python package for visualizing line-wise MSI data across MS1, MS2, MRM, and IMS modes
- [MSight](https://github.com/laurenfields/MSIght) - Automated integration of MSI with histological stains and LC-MS/MS peptide IDs
- [MSIpixel](https://github.com/gmartano/MSIpixel) - Automated pipeline for compound annotation and quantitation in DDA MSI
- [msiFlow](https://github.com/Immunodynamics-Engel-Lab/msiflow) - End-to-end software for reproducible multimodal MALDI MSI and microscopy analysis
- [MSroi](http://www.mcrals.info/) - Software for compressing and preprocessing MS data from direct infusion, LC-MS, or imaging
- [Multi‑MSI Processor](https://github.com/multimsiproc/MMP) - Software for batch MSI data analysis, biomarker discovery, and spatial workflows
- [Pew2](https://github.com/djdt/pewpew) - Open-source LA-ICP-MS image processing software
- [rMSIannotation](https://github.com/prafols/rMSIproc) - Algorithm for annotating carbon isotopes and adducts in MSI without compound libraries
- [S2IsoMEr](https://github.com/alexandrovteam/S2IsoMEr) - R package propagating isomeric/isobaric ambiguities into enrichment analysis for spatial metabolomics
- [SagMSI](https://github.com/BioNet-XMU/SagMSI) - Unsupervised graph convolution network for spatial segmentation of MSI tissue data
- [SmartGate](https://github.com/zhanglabtools/SmartGate) - AI tool for iterative peak selection and spatial tissue segmentation in imaging MS
- [SONAR-MSI](https://github.com/GrantBellic/SONAR-MSI) - Workflow integrating synchronized ion acquisition with pseudo-MSI and deep learning
- [SpaceM](https://github.com/alexandrovteam/SpaceM) - In situ single-cell metabolomics via microscopy and MALDI imaging MS
- [SpatialMETA](https://github.com/WanluLiuLab/SpatialMETA) - Variational autoencoder for cross-modal integration of spatial transcriptomics and metabolomics
- [SSN](https://github.com/LanekoffLab/i2i) - Tandem MS imaging annotation using spatial correlations among product ions
- [subMALDI](https://github.com/wesleyburr/subMaldi) - Open-framework R package for organizing, preprocessing, and normalizing MS spectral data

---

## 🌬️ Ion Mobility MS

- [AllCCS2](http://allccs.zhulab.cn/) - CCS prediction database
- [AutoCCS](https://github.com/PNNL-Comp-Mass-Spec/AutoCCS) - Automated CCS calculation
- [AutonoMS](https://autonoms.readthedocs.io) - Integrates analytical equipment into automated laboratory workflows
- [CCS Predictor 2.0](https://github.com/facundof2016/CCSP2.0) - Open-source SVR models to predict CCS values and filter false positives
- [CCSfind](https://github.com/sangeeta97/ccs_find) - Tool for building CCS databases from experimental LC-IM-MS measurements
- [GCIMS](https://github.com/sipss/GCIMS) - R package for untargeted GC-IMS data processing including denoising, alignment, and clustering
- [HyperCCS](https://github.com/NeoNexusX/HyperCCS) - CCS prediction framework powered by chemical large language models
- [Met4DX](https://met4dx.zhulab.cn/) - LC-IM-MS/MS data processing for 4D metabolomics
- [MobiLipid](https://github.com/FelinaHildebrand/MobiLipid) - CCS quality control for IM-MS lipidomics using ¹³C-labeled lipid standards
- [MOCCal](https://github.com/HinesLab/MOCCal) - Automated multiomics CCS calibration for traveling wave ion mobility
- [Mol2CCS](https://github.com/enveda/ccs-prediction) - GNN-based CCS prediction benchmarking with molecular fingerprints and confidence models
- [OpenTIMS, TimsPy, and TimsR](www.github.com/michalsta/opentims) - TIMS data tools
- [PACCS](https://github.com/yuxuanliao/PACCS) - Voxel-projected area and GNN-based CCS prediction for diverse molecules
- [PNNL PreProcessor](https://github.com/PNNL-Comp-Mass-Spec/PNNL-PreProcessor/releases/tag/v6.0) - Converts LC-IM-MS data to LC-MS format for existing analysis software
- [POMICS](https://www.pomics.org/) - ML-assisted ab initio CCS computation using DFT and conformational modeling
- [PubChemLite plus CCS](https://pubchemlite.lcsb.uni.lu) - Curated chemical collection for exposomics with predicted CCS values from PubChem
- [SigmaCCS](https://github.com/zmzhang/SigmaCCS) - GNN-based CCS prediction from 3D conformers with 282M-value in-silico database
- [Snakemake CCS](https://github.com/DasSusanta/snakemake_ccs) - Automated ab initio CCS prediction from SMILES using QM on HPC via Snakemake
- [SSN-CCSPIF](https://github.com/Anqi-Guo0105/Trendline-based-IM-MS-Feature-Filtering-Software-V1.0) - Structural similarity networking with CCS prediction interval filtering for IM-MS

---

## 🧮 Isotopic

- [Aerith](https://github.com/thepanlab/Aerith) - Visualization of isotopic patterns in stable isotope probing MS data
- [CIL-ExPMRM](http://www.exposomemrm.com/) - Ultrasensitive chemical isotope labeling pseudo-MRM platform for exposomic suspect screening
- [FAMetA](https://www.fameta.es) - Estimates fatty acid import, de novo lipogenesis, elongation, and desaturation from ¹³C profiles
- [IsoPairFinder](https://github.com/DoddLab/IsoPairFinder) - Identifies pathway intermediates from paired ¹²C/¹³C metabolomics data
- [isopair](https://github.com/kilgain/MassSpec/tree/main) - Suite of tools for LC-MS peak picking using isotopologue pairs
- [Isoreader](https://github.com/isoverse/isoreader) - R package to read stable isotope ratio MS data from common instrument formats
- [IsoSolve](https://pypi.org/project/isosolve/) - Integrative framework consolidating MS and NMR isotopic measurements
- [khipu](https://github.com/shuzhao-li-lab/khipu) - Annotates adducts, fragments, and isotopes to infer neutral mass in untargeted metabolomics
- [Khipu-web](https://metabolomics.cloud/khipu) - Web app for annotating isotope-labeled metabolites in stable isotope tracing
- [SIMPEL](https://github.com/SIMPELmetabolism/SIMPEL) - R package for automated isotopologue mining from HRMS data for metabolic flux analysis

---

## 🔭 Large Scale

- [LargeMetabo](https://github.com/LargeMetabo/LargeMetabo) - R package integrating multiple experiments for biomarker identification at scale
- [peakPantheR](https://bioconductor.org/packages/peakPantheR/) - Peak integration at scale

---

## 🧬 Lipidomics

- [ADViSELipidomics](https://github.com/ShinyFabio/ADViSELipidomics) - Shiny app for preprocessing, analyzing, and visualizing lipidomics data with batch correction
- [BATL](https://complimet.ca/batl/) - Bayesian naïve Bayes classifier annotating targeted lipidomics peaks using RT and intensity features
- [DBLiPro](http://lipid.cloudna.cn/home) - Webserver for human lipid metabolism knowledge base and integrative lipidomic analysis
- [DIATAGe](https://github.com/Velenosi-Lab/DIATAGeR) - R package for triacylglycerol identification in DIA-based lipidomics using target-decoy scoring
- [Doxlipid](https://figshare.com/articles/journal_contribution/_em_strong_In_silico_strong_em_strong_-predicted_dynamic_oxlipidomics_MS_MS_library_High-throughput_discovery_and_characterization_of_unknown_oxidized_lipids_strong_/23700177/4) - In silico MS/MS library for annotation of 32 subclasses of oxidized lipids
- [LINEX2](https://gitlab.lrz.de/lipitum-projects/linex) - Lipid network analysis framework inferring enzymatic activity from lipidomics via graph algorithms
- [LipiDetective](https://github.com/LipiTUM/lipidetective) - Deep learning transformer for molecular lipid species identification from tandem MS
- [LipiDex 2](https://github.com/coongroup/LipiDex-2) - Lipidomics processing with MSn tree-based fragmentation, spectral matching, and automated QC
- [LipiDisease](http://cbdm-01.zdv.uni-mainz.de:3838/piyusmor/LipiDisease/) - Literature-mining web server associating lipid sets with diseases via PubMed analysis
- [Lipid Spectrum Generator](https://github.com/98104781/LSG) - Open-source in silico spectral library generator for lipid identification in LC-MS/MS
- [Lipid Wizard](https://github.com/RaoboXu/Lipidwizard) - Lipid assignment and isotopic peak stripping in 2D-LC-MS using ECN-based RT prediction
- [Lipid4Danalyzer](http://lipid4danalyzer.zhulab.cn/) - Lipid identification via 4D alignment of RT, MS1, MS/MS, and CCS
- [LipidA-IDER](https://github.com/Systems-Biology-Of-Lipid-Metabolism-Lab/LipidA-IDER) - Automated lipid A structure annotation from HR-MS/MS for Gram-negative bacteria
- [LipidCruncher](https://lipidcruncher.org/) - Open-source web platform for processing, visualizing, and analyzing lipidomic data
- [Lipidepifind](https://github.com/Xinguang-Liu/Lipidepifind) - Identifies structurally modified lipid epi-metabolites using metabolic network expansion
- [LipidFinder 2.0](https://github.com/ODonnell-Lipidomics/LipidFinder) - Lipidomics pipeline with artefact filtering, in-source fragment removal, and FDR assessment
- [LipidIN](https://github.com/LinShuhaiLAB/LipidIN) - Platform-independent lipid annotation framework with 168.5M fragmentation library
- [LipidOne 2.0](https://lipidone.eu) - Web tool analyzing lipidomic data at class, species, and building block levels
- [LipidOracle](https://hub.docker.com/r/zambonilab/lipidoracle) - Infers lipid double bond and sn-positions from electron-induced dissociation MS
- [LipidOz](https://github.com/PNNL-m-q/lipidoz) - Automated determination of lipid C=C double bond positions from OzID IMS-MS data
- [LipidQuant 1.0](https://holcapek.upce.cz/lipidquant.php) - Automated lipidomic quantitation using lipid class separation with high-resolution MS
- [LipidQuant 2.1](https://zenodo.org/record/10364585) - MATLAB software for high-throughput lipid identification and quantification
- [LipidSig](chenglab.cmu.edu.tw/lipidsig) - Web server for lipidomic profiling, differential expression, correlation, and network analysis
- [LipidSig 2.0](https://lipidsig.bioinfomics.org/) - Auto-identifies lipid species and assigns 29 characteristics for integrated lipidomic analysis
- [LipidSigR](https://github.com/BioinfOMICS/LipidSigR) - R package companion to LipidSig for customizable lipidomics analysis and visualization
- [LipidSpace](https://lifs-tools.org/tools/lipidspace.html) - Tool comparing lipidomes using graph-based structural similarity models and ML analysis
- [LipidSuite](https://suite.lipidr.org/) - End-to-end web server for differential lipidomics preprocessing and enrichment analysis
- [LIPID MAPS](https://www.lipidmaps.org/) - Community-standard lipid classification system with databases, pathways, and bioinformatics
- [LipoCLEAN](https://github.com/stavis1/LipoCLEAN) - ML filter improving untargeted lipid ID confidence by rescoring with isotope and chromatographic metrics
- [LORA](https://lora.metabolomics.fgu.cas.cz/) - Web app for lipid over-representation analysis based on structural hierarchy
- [LPPtiger2](https://github.com/LMAI-TUD/lpptiger2) - Cross-platform software for predicting and identifying oxidized complex lipids
- [MS2Lipid](https://github.com/systemsomicslab/ms2lipid/tree/main/notebooks) - ML model predicting 97 lipid subclasses from MS/MS spectra with 97.4% accuracy
- [MS-RIDD](https://github.com/Metabolomics-CompMS/MS-RIDD) - C=C position-resolved untargeted lipidomics combining oxygen attachment dissociation with computational MS
- [Neurolipid Atlas](https://neurolipidatlas.com/) - Lipidomics data commons for neurodegenerative diseases using iPSC-derived cell data
- [RefLAS](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR002159) - Curated LC-MS lipidomics database based on NIST SRM 1950 human plasma
- [RPLC-IMS-CID-MS Lipid Database](https://tarheels.live/bakerlab/databases/) - Multidimensional database with RT, CCS, and m/z for 877 lipids across 24 classes
- [RTStaR](https://www.neurolipidomics.ca/rtstar/rtstar.html) - Retention time standardization and registration for nLC-nESI-MS/MS lipidomics
- [SimLipid](https://www.premierbiosoft.com/lipid/index.html) - Software for large-scale LC-MS untargeted lipidomics profiling with customizable library

---

## 🕸️ Metabolic Networks

- [DNEA](https://github.com/Karnovsky-Lab/DNEA/) - Differential network enrichment
- [FNICM](https://github.com/LiQi94/FNICM) - Framework for identifying disease-associated core metabolites via perturbed metabolic subnetworks
- [KGMN](https://github.com/ZhuMetLab/MetDNA2) - Multi-layer network for annotating unknown metabolites using metabolic reactions and MS/MS similarity
- [Kiphynet](http://pallab.cds.iisc.ac.in/kiphynet/) - Web app simulating metabolic transport and reactions in multi-scale microvascular models
- [M2R](https://github.com/e-weglarz-tomczak/m2r) - Python add-on to cobrapy for integrating gut microbiota metabolism into human metabolic models
- [MetaboliticsDB](https://github.com/alperdokay/metabolitics-frontend) - Database of metabolomics analyses with network-based flux analysis and AI disease association
- [MetNet](https://github.com/simeoni-biolab/MetNet) - Automated metabolic network reconstruction and comparison
- [MicroMap](https://dataverse.harvard.edu/dataverse/micromap) - Curated visualization of metabolism across 257,000+ microbial genome-scale metabolic reconstructions
- [MInfer](https://github.com/cellbiomaths/MInfer) - Derives metabolic interaction networks via Jacobian analysis
- [MINNO](https://lewisresearchgroup.github.io/MINNO/) - Web tool for refining metabolic networks using empirical metabolomics data in nonmodel organisms
- [NetAurHPD](https://github.com/TamirBar-Tov/NetAurHPD-Network-Auralization-Hyperlink-Prediction-Method) - Framework using graph auralization to predict metabolic pathways from correlation networks
- [Recon8D](https://github.com/sriram-lab/Recon8D) - ML-based metabolic regulome network from eight omics layers across ~1,000 cancer cell lines
- [Thermo-Flux](https://github.com/molecular-systems-biology/thermo-flux) - Semi-automated Python package adding thermodynamic constraints to genome-scale metabolic models

---

## 🗃️ Metadata

- [MatrixLM](https://github.com/senresearch/MatrixLM.jl) - Julia package for matrix linear models linking metabolite composition to sample traits
- [MetaXtract](https://github.com/Rappsilber-Laboratory/MetaXtract) - Tool extracting LC-MS acquisition metadata from Thermo Fisher raw files
- [PeakForest](https://github.com/peakforest/peakforest-webapp) - Open-source infrastructure for managing and sharing NMR and MS spectral data
- [SMetaS](https://github.com/metabolomics-us/metadatastandardizer) - Tool for automatically standardizing metabolomics sample metadata using ontology vocabularies

---

## 🔀 Multifunctional

- [DEIMoS](https://github.com/pnnl/deimos) - Data extraction for ion mobility
- [eMZed 3](https://emzed.ethz.ch) - Interactive analysis
- [FERMO](https://fermo.bioinformatics.nl/) - Feature-based molecular networking
- [iMAP](https://imap.metaboprofile.cloud/) - Web server for metabolomics data analysis with feature selection, network analysis, and visualization
- [maplet](https://github.com/krumsieklab/maplet) - R package for reproducible metabolomics statistical pipelines
- [Mass-Suite](https://github.com/XiminHu/mass-suite) - Python package for HRMS-based non-targeted analysis for water quality and environmental applications
- [matchms](https://github.com/matchms/matchms) - Python library for mass spectral data curation and similarity scoring
- [MargheRita](https://github.com/emosca-cnr/margheRita) - R package for LC-MS/MS SWATH metabolomics data analysis and metabolite identification
- [MetaboAnalyst 5.0](http://www.metaboanalyst.ca/) - Web platform for LC-HRMS processing, statistical analysis, and functional interpretation
- [MetaboAnalyst 6.0](https://www.metaboanalyst.ca) - Unified platform adding MS2 annotation, exposomics support, and dose-response analysis
- [MetaboAnalystR 4.0](http://www.metaboanalyst.ca/) - R package providing a unified LC-MS workflow from raw spectra to compound identification
- [MetaboLink](https://github.com/anitamnd/MetaboLink) - Web app for LC-MS metabolomics post-processing including correction, filtration, and normalization
- [MetDNA3](http://metdna.zhulab.cn/) - Two-layer interactive networking for metabolite annotation
- [MetMiner, MDAtoolkits](https://github.com/ShawnWx2019/MetMiner) - User-friendly pipeline for large-scale plant metabolomics data analysis
- [MetaProViz](https://saezlab.github.io/MetaProViz/) - R package for metabolomics functional analysis and visualization with curated prior knowledge
- [mpactR](https://www.mums2.org/mpactr/) - R package automating LC-MS/MS data pre-processing for microbiological samples
- [MS-DIAL 5](https://github.com/systemsomicslab/MsdialWorkbench) - Software for multimodal MS data mining with lipidome structural elucidation and imaging
- [MSOne](https://msone.claritybiosystems.com) - AI-powered web-based end-to-end LC-MS analysis suite
- [MZmine 3](https://github.com/mzmine/mzmine) - Scalable MS data analysis platform supporting LC/GC-MS, ion mobility, and MS imaging
- [NP3 MS Workflow](https://github.com/danielatrivella/NP3_MS_Workflow) - Open-source LC-MS/MS workflow for ranking bioactive natural products from complex mixtures
- [OpenMS 3](https://github.com/OpenMS/OpenMS) - Open-source framework of MS algorithms enabling reproducible large-scale data analysis
- [OpenMS WebApps](https://github.com/OpenMS/streamlit-template) - Framework for building user-friendly MS web applications using Streamlit and OpenMS
- [patRoon 2](https://github.com/rickhelmus/patRoon) - R-based non-target analysis platform with automated transformation product screening for HRMS
- [PMart](https://github.com/pmartR/PMart_ShinyApp) - Web tool for reproducible QC, statistical analysis, and visualization of multi-omics data
- [POMAShiny](https://github.com/nutrimetabolomics/POMAShiny) - Web-based tool for statistical analysis and visualization of metabolomics and proteomics data
- [Punc'data](https://github.com/WTVoe/puncdata) - JavaScript tool for molecular formula assignment and visualization of ultrahigh-resolution MS data
- [Rodin](https://rodin-meta.com/) - Python app for streamlined metabolomics preprocessing, statistics, and pathway analysis
- [SLAW](https://github.com/adelabriere/SLAW) - Scalable, self-optimizing workflow for large-scale untargeted LC-MS metabolomics
- [TidyMass](https://www.tidymass.org/) - R-based ecosystem for reproducible and traceable LC-MS untargeted metabolomics analysis
- [TidyMass2](https://github.com/tidymass/tidymass) - Enhanced LC-MS framework adding metabolite origin inference and functional module analysis
- [TraceMetrix](https://www.biosino.org/tracemetrix/) - Web platform providing interactive traceability across the full metabolomics pipeline
- [UmetaFlow](https://github.com/biosustain/snakemake_UmetaFlow) - Snakemake workflow for untargeted metabolomics with preprocessing, annotation, and networking
- [XCMS-METLIN](https://xcmsonline.scripps.edu/) - Integrates XCMS data processing with METLIN's 935,000+ experimental MS/MS spectra

---

## 🔗 Multiomics

- [AgeAnnoMO](https://relab.xidian.edu.cn/AgeAnnoMO/#/) - Multi-omics knowledgebase annotating aging-related genes, proteins, and metabolites across 50 species
- [BATMAN‑TCM 2.0](http://bionet.ncpsb.org.cn/batman-tcm/) - Database of ~2.3 million TCM ingredient–target protein interactions for pharmacology
- [BnIR](http://yanglab.hzau.edu.cn/BnIR) - Multi-omics database for rapeseed integrating genomics, transcriptomics, and metabolomics
- [haCCA](https://github.com/LittleLittleCloud/haCCA) - Workflow integrating spatial transcriptomics and MALDI-MSI metabolomics via correlated features
- [HoloFoodR](https://doi.org/10.18129/B9.bioc.HoloFoodR) - R/Bioconductor package for accessing and integrating holo-omics data from HoloFood
- [INTEGRATE](https://github.com/qLSLab/integrate) - Pipeline integrating metabolomics and transcriptomics to characterize metabolic regulation
- [iSODA](http://isoda.online/) - Web app for interactive single- and multi-omics data analysis with visualization
- [iTraNet](https://itranet.streamlit.app) - Web app for visualizing trans-omics networks spanning regulatory and metabolic layers
- [jMorp](https://jmorp.megabank.tohoku.ac.jp) - Japanese multi-omics reference panel providing metabolome, genome, and metagenome data
- [MetaNet](https://github.com/Asa12138/MetaNet) - High-performance R package for constructing and analyzing biological networks from multi-omics
- [MiMeNet](https://github.com/YDaiLab/MiMeNet) - Neural network for predicting metabolite abundances from microbiome composition
- [MODMS](https://modms.lzu.edu.cn) - Multi-omics database for alfalfa integrating genomes, transcriptomes, and metabolomics
- [MPOD](http://mpod.ict.cn) - Integrated multi-omics database and analysis platform for medicinal plants
- [MultiOmicsIntegrator](https://github.com/ASAGlab/MOI--An-integrated-solution-for-omics-analyses) - Nextflow pipeline for integrated analysis of RNA-seq, proteomics, and metabolomics
- [OmicsAnalyst](https://www.omicsanalyst.ca) - Web platform for multi-omics integration via correlation networks and dimension reduction
- [Omics Dashboard](https://biocyc.org/dashboard/dashboard-help.html) - Interactive tool for exploring metabolomics and multi-omics data by cellular system hierarchies
- [OmicsNet 2.0](https://www.omicsnet.ca) - Web tool for building multi-omics networks with 3D layouts and microbiome support
- [Paired Omics](https://github.com/iomega/paired-data-form) - Community platform documenting links between metabolome and genome data for natural products
- [SOmicsFusion](https://doi.org/10.5281/zenodo.7700528) - Toolbox for coregistering spatial metabolomics with biomedical imaging modalities like MRI
- [SVAtlas](https://www.svatlas.org/) - Database integrating single extracellular vesicle multi-omics data across 31 diseases
- [TurbOmics](https://proteomics.cnic.es/TurboPutative/TurbOmicsApp.html) - Web platform for multi-omics integration of metabolomics, proteomics, and transcriptomics
- [TurboPutative](https://proteomics.cnic.es/TurboPutative/) - Web server for reducing and prioritizing putative annotations in untargeted LC-MS

---

## 🔬 NMR

- [A-SIMA/A-MAP](https://poky.clas.ucdenver.edu) - Toolkit for NMR metabolite identification and multivariate analysis in the POKY suite
- [AQuA](https://pmc.ncbi.nlm.nih.gov/articles/instance/8253485/bin/ac0c04233_si_001.pdf) - Fast automated quantification of metabolites from 1D ¹H NMR spectra
- [ASICS](https://github.com/GaelleLefort/ASICS) - Joint automatic identification and quantification of metabolites across ¹H NMR spectral sets
- [Bucket Fuser](https://www.mdpi.com/article/10.3390/metabo12090812/s1) - Regularized regression for extracting metabolite signals from 1D ¹H NMR spectra
- [CASMDB](https://ccpn.ac.uk/software/analysismetabolomics/) - Integrated simulated metabolite reference database for 1D ¹H NMR metabolomics
- [COLMAR1d](https://spin.ccic.osu.edu/index.php/colmar1d) - Web server for automated quantitative 1D NMR metabolomics at arbitrary field strengths
- [COLMARppm](https://spin.ccic.osu.edu/index.php/colmar) - Web server for rapid prediction of ¹H and ¹³C NMR chemical shifts
- [COLMARvista](https://github.com/lidawei1975/colmarvista) - Open-source JavaScript tool for processing and visualizing 2D and pseudo-3D NMR spectra
- [CReSS](https://github.com/Qihoo360/CReSS) - Cross-modal retrieval matching ¹³C NMR spectra to molecular structures
- [DeepSAT](https://github.com/mwang87/DeepSAT) - Deep learning for molecular structure elucidation from NMR spectra
- [FlavorFormer](https://github.com/yfWang01/FlavorFormer) - Deep learning model for identifying compounds in flavor mixtures from ¹H NMR spectra
- [InRA](https://github.com/InRA-Software/InRA) - Automated interval-based tool for untargeted ¹H NMR multivariate fingerprinting
- [LAMAIS](https://github.com/Ariel-foliage/LAMAIS) - Automated qualitative analysis of 1D ¹H NMR mixtures
- [MADByTE](https://github.com/liningtonlab/madbyte) - NMR platform using 2D spectra to create chemical similarity networks for natural products
- [Magnetstein](https://github.com/ElsevierSoftwareX/SOFTX-D-25-00358) - Web application for quantitative NMR mixture analysis using Wasserstein distance
- [MagMet](http://magmet.ca/) - Web server for automated NMR processing and quantification of serum/plasma metabolites
- [MagMet-F](https://www.magmet.ca/) - Automated NMR identification and quantification of metabolites in human fecal extracts
- [mcfNMR](https://github.com/GeoMetabolomics-ICBM/mcfNMR) - Network flow method for compound identification in NMR mixture spectra without peak-picking
- [MetaboLabPy](https://github.com/ludwigc/metabolabpy) - Open-source Python package for 1D/2D NMR spectra processing
- [MetAssimulo 2.0](https://github.com/yanyan5420/MetAssimulo_2) - Web app for simulating realistic 1D and 2D metabolomic ¹H NMR spectra
- [MixONat](https://sourceforge.net/projects/mixonat/) - Software for ¹³C NMR-based dereplication of natural products in mixtures
- [MultiNMRFit](https://github.com/NMRTeamTBI/MultiNMRFit/) - Open-source tool for high-throughput fitting and parameter extraction from 1D NMR spectra
- [NAPROC-13](https://github.com/DIFACQUIM/naproc13_characterization) - Database of ¹³C NMR spectroscopic data for over 25,000 natural products
- [NMR2Struct](https://github.com/MarklandGroup/NMR2Struct) - ML framework for predicting molecular structures from 1D NMR spectra
- [NMRformer](https://github.com/zza1211/NMRformer) - Transformer-based deep learning for NMR peak assignment and metabolite identification
- [NMRFx](https://github.com/nanalysis/nmrfx) - Software for NMR data processing, visualization, peak assignment, and structure calculation
- [NMRhub](https://nmrhub.org/) - Integrated ecosystem spanning the complete NMR data lifecycle from acquisition to deposition
- [NMRInversions.jl](https://github.com/aris-mav/NMRInversions.jl) - Julia package for time-domain NMR data processing, inversions, and visualization
- [NMRium](https://www.nmrium.org/teaching) - Web-based platform for processing, interpreting, and teaching 1D and 2D NMR spectra
- [NMR molecular networking](https://github.com/enveda/NMR-Networking) - Framework adapting molecular networking to HSQC NMR spectra for structure annotation
- [NMRphasing](https://github.com/ajiangsfu/NMRphasing) - R package for automated NMR phase error correction using nonlinear shrinkage
- [NMRQNet](https://github.com/LiuzLab/NMRQNet) - Deep learning for automatic NMR identification and quantification of plasma metabolites
- [nmRanalysis](https://github.com/EMSL-Computing/nmRanalysis) - Web application for semi-automated NMR metabolite profiling with ML-driven identification
- [NP-MRD](https://np-mrd.org/) - Comprehensive database of NMR spectra and data for ~282,000 natural products
- [PRIMA-Panel](https://github.com/funkam/PRIMA) - R Shiny tool for assessing preanalytical variation impact on NMR metabolomics data
- [PROSPRE](https://prospre.ca/) - Deep learning tool for accurate ¹H NMR chemical shift prediction across solvents
- [Protomix](https://github.com/mzniber/Protomix) - Python package for preprocessing 1D ¹H-NMR metabolomics data with visualization
- [pSCNN](https://github.com/yuxuanliao/pSCNN) - Pseudo-Siamese CNN for compound identification in NMR mixture spectra
- [PyINETA](https://github.com/edisonomics/PyINETA) - Python platform integrating INADEQUATE and ¹³C-JRES NMR for metabolite annotation
- [ROIAL-NMR](https://github.com/Leo-Cheng-Lab/ROIAL-NMR) - Python tool for identifying metabolites from proton NMR regions-of-interest using HMDB
- [SAND](https://github.com/edisonomics/SAND) - Automated NMR spectral decomposition for metabolomics
- [SMolESY](https://github.com/pantakis/SMolESY-select) - Automated algorithm for metabolite assignment and quantification in blood ¹H NMR spectra
- [ukbnmr](https://github.com/sritchie73/ukbnmr/) - R package for quality control and technical variation removal in UK Biobank NMR data

---

## 🌱 Organism Specific

- [Aspergillus Metabolome Database](https://github.com/albertogilf/ceuMassMediator/tree/master/CMMAspergillusDB) - Database for m/z-based metabolite annotation in MS studies of Aspergillus
- [FoodAtlas](http://foodatlas.ai) - Knowledge graph of food-chemical relationships extracted from literature using LLMs
- [GMDP Database](https://lccl.shinyapps.io/GMDP/) - Web app for exploring metabolomics data from patients with inborn errors of metabolism
- [gutMGene v2.0](http://bio-computing.hrbmu.edu.cn/gutmgene) - Database of associations among gut microbes, their metabolites, and host genes
- [gutSMASH](https://gutsmash.bioinformatics.nl/) - Tool for systematic profiling of primary metabolic pathways in gut microbiome genomes
- [GutUDB](https://intestine.splicedb.net/) - Comprehensive multiomics database for analysis and visualization of intestinal diseases
- [HORDB](https://github.com/CPU-HORDB/HORDB) - Database of peptide hormones with activity, structure, and pharmaceutical information
- [IDBac](https://github.com/chasemc/IDBacApp) - MALDI-TOF MS platform for microbial strain prioritization based on phylogeny and metabolites
- [LiqBioer](http://www.medsysbio.org:8080/LiqBioer/) - Curated database of cancer biomarkers in body fluids from liquid biopsy studies
- [MDSi](http://sky.sxau.edu.cn/MDSi.htm) - Multi-omics database for foxtail millet integrating genomics, transcriptomics, and metabolomics
- [MetaDb](http://medmetadb.ynau.edu.cn) - Database of plant metabolites and their regulation, emphasizing medicinal plants
- [MetaboSeek](https://github.com/mjhelf/Metaboseek) - Open-source platform for untargeted comparative LC-MS metabolomics and feature annotation
- [MiMIR](https://github.com/DanieleBizzarri/MiMIR) - R Shiny app for projecting pre-trained risk models onto Nightingale NMR metabolomics data
- [mVOC 4.0](http://bioinformatics.charite.de/mvoc) - Database of ~3,500 microbial volatile organic compounds with MS matching
- [MyxoDB](https://figshare.com/collections/Constructing_a_Myxobacterial_Natural_Product_Database_to_Facilitate_NMR-Based_Metabolomics_Bioprospecting_of_Myxobacteria/6468946) - Database of myxobacterial natural products for NMR-based dereplication
- [Omics Untargeted Key Script](https://github.com/plyush1993/OUKS) - R-based scripts for comprehensive untargeted LC-MS metabolomics data processing
- [OrchidMD](https://www.orchidcomics.com) - Multi-omics database for orchids integrating five omics layers across 213 species
- [PaIRKAT](https://github.com/CharlieCarpenter/PaIRKAT) - Kernel regression method incorporating pathway topology for metabolomics association testing
- [PeanutOmics](https://cgm.sjtu.edu.cn/PeanutOmics/) - Multi-omics platform integrating transcriptomic, proteomic, and metabolomic data for peanut
- [PSC-db](http://pscdb.appsbio.utalca.cl/) - Database of plant secondary compounds with 3D structures and bioactivity data
- [SalivaDB](https://webs.iiitd.edu.in/raghava/salivadb/) - Comprehensive database of salivary biomarkers in humans for diagnostic research
- [SmilODB](http://www.isage.top:56789) - Multi-omics database for Salvia miltiorrhiza integrating genomic and metabolomic data
- [StreptomeDB 4.0](https://streptomedb.vm.uni-freiburg.de/streptomedb/) - Database of >8,500 streptomycete natural products with spectral and interaction data
- [The Molecular Human](https://littleswissriver.shinyapps.io/generate_network/) - Multiomics network integrating 18 platforms across biofluids to portray molecular phenotypes

---

## 🛤️ Pathway, Enrichment and Ontology Tools

- [ChemFOnt](https://www.chemfont.ca/) - Hierarchical ontology describing functions and actions of over 341,000 biologically important chemicals
- [EnrichMET](https://github.com/biodatalab/enrichmet) - Web tool for pathway enrichment analysis of metabolomics data with multiple databases
- [GINv2.0](https://github.com/BIGchix/GINv2.0) - Knowledge base integrating signaling and metabolic networks for pathway analysis
- [IDSL.GOA](https://goa.idsl.me/) - Gene ontology enrichment analysis extending biological interpretation beyond pathway-level metabolomics
- [iMSEA](https://github.com/BioNet-XMU/iMSEA) - Metabolite set enrichment analysis for deciphering drug interaction mechanisms via metabolomics
- [Lilikoi V2.0](https://github.com/lanagarmire/lilikoi2) - Personalized pathway-based method for disease classification using metabolomics data
- [massDatabase](https://massdatabase.tidymass.org/) - R package integrating compound and pathway databases for metabolite annotation and enrichment
- [MBROLE3](http://csbg.cnb.csic.es/mbrole3) - Web server for enrichment analysis of metabolite sets using annotations from multiple databases
- [metGWAS 1.0](https://github.com/saifurbd28/metGWAS-1.0) - Tool integrating GWAS summary statistics with metabolomics to identify disease-linked pathways
- [metLinkR](https://github.com/metLinkR/metLinkR) - R tool for cross-linking metabolite identifiers across spectral, structural, and pathway databases
- [MetChem](https://CRAN.R-project.org/package=MetChem) - R package for chemical-structure-based analysis and interpretation of metabolomic profiling data
- [MIMOSA2](http://www.borensteinlab.com/software_MIMOSA2.html) - Tool linking microbial community composition to metabolite variation in paired studies
- [MiNEApy](https://github.com/vpandey-om/mineapy) - Python tool for metabolic network enrichment analysis using elementary flux modes
- [MS2MP](https://github.com/ucasaccn/MS2MP) - Deep learning framework predicting KEGG pathways directly from MS/MS spectra
- [ORA](https://github.com/cwieder/metabolomics-ORA) - Best-practice guidelines and benchmarking for over-representation pathway analysis in metabolomics
- [PALS](https://pals.glasgowcompbio.org/app/) - Tool for scoring pathway-level activity changes from grouped metabolite intensities
- [PaintOmics 4](https://paintomics.org) - Web server for integrative multi-omics visualization on KEGG, Reactome, and MapMan pathways
- [PathBank 2.0](https://pathbank.org/) - Comprehensive pathway database providing detailed metabolite and protein pathway maps
- [SBGNview](https://github.com/datapplab/SBGNview) - R package for omics data visualization on pathway maps using standard SBGN format
- [secCellFie](https://github.com/LewisLabUCSD/secCellFie) - Extension of CellFie quantifying secretory pathway metabolic activity from transcriptomics
- [Viime-Path](https://github.com/girder/viime-path) - Tool for visualizing and integrating enrichment results of altered metabolic pathways
- [WebGestalt](https://www.webgestalt.org) - Web platform for integrated enrichment analysis across genomics, proteomics, and metabolomics
- [Xconnector](https://github.com/Proteomicslab57357/Xconnector) - Tool for connecting and bridging diverse metabolomics databases and their identifiers

---

## 🔎 Patterns

- [MCnebula](https://mcnebula.org/) - Multiple classes nebula
- [Metabokiller](https://pypi.org/project/Metabokiller/) - Ensemble classifier for carcinogen prediction from metabolite structures
- [MS2LDA 2.0](https://github.com/vdhooftcompmet/MS2LDA) - LDA for MS/MS
- [mzBucket](https://github.com/hildebrandtlab/mzBucket) - m/z bucketing

---

## ⚙️ Pre-processing

- [3D-MSNet](https://github.com/CSi-Studio/3D-MSNet) - 3D feature detection directly on lossless profile MS signals
- [AriumMS](https://github.com/AdrianHaun/AriumMS/) - Arium MS processing
- [asari](https://github.com/shuzhao-li-lab/asari) - Feature extraction
- [AVIR](https://github.com/HuanLab/AVIR) - Adduct and in-source verification
- [CorrDIA](https://github.com/zhangxiang1205/CorrDIA) - Correlation-based DIA
- [CycloBranch](https://ms.biomed.cas.cz/cyclobranch) - Cyclopeptide analysis
- [DaDIA](https://github.com/HuanLab/DaDIA) - Data-dependent acquisition
- [DBDIpy](https://github.com/leopold-weidner/DBDIpy) - DBDI processing
- [DecoID](https://github.com/pattilab/DecoID/releases) - Deconvolution ID
- [DeepMSProfiler](https://github.com/yjdeng9/DeepMSProfiler) - End-to-end deep learning on raw metabolomic MS signals
- [DIAMetAlyzer](https://github.com/OpenMS/DIAMetAlyzer) - DIA analysis
- [DNMS2Purifier](https://github.com/HuanLab/DNMS2Purifier) - MS/MS purification
- [DuReS](https://github.com/BiosystemEngineeringLab-IITB/dures) - Duplicate removal
- [Eclipse](https://github.com/broadinstitute/bmxp/tree/main/bmxp/eclipse) - Eclipse processing
- [EVA](https://github.com/HuanLab/EVA) - CNN-based metabolic feature recognition from extracted ion chromatogram peak shapes
- [Finnee 2024](https://github.com/glerny/Finnee2024) - Finnee update
- [G-Aligner](https://github.com/CSi-Studio/G-Aligner) - Global alignment
- [GCMSFormer](https://github.com/zxguocsu/GCMSFormer) - Transformer for GC-MS
- [HeuSMA](https://github.com/Lacterd/HeuSMA) - Heuristic spectral matching
- [IDSL.IPA](https://cran.r-project.org/package=IDSL.IPA) - Isotope pattern analysis
- [ImpLiMet](https://complimet.ca/shiny/implimet/) - Imputation for lipidomics
- [ISFrag](https://github.com/HuanLab/ISFrag) - In-source fragmentation
- [IsoFusion](https://github.com/xfcui/IsoFusion) - Isotope fusion
- [JPA](https://github.com/HuanLab/JPA) - Joint peak alignment
- [LAGF](https://github.com/zsspython/LAGF) - Lag filtering
- [LongitudinalProf-MSMS-Method](https://github.com/academicCQMD/LongitudinalProf-MSMS-Method) - Longitudinal profiling
- [maplet](https://github.com/krumsieklab/maplet) - R package for reproducible metabolomics statistical pipelines
- [MassDash](https://github.com/Roestlab/massdash) - Mass dashboard
- [MassLite](https://github.com/chemzzchem/MassLite/tree/main) - Lightweight MS
- [MESSES](https://github.com/MoseleyBioinformaticsLab/MESSES) - Metadata extraction
- [metabCombiner](https://bioconductor.org/packages/metabCombiner/) - Dataset combining
- [metabCombiner 2.0](https://metabcombiner.medicine.umich.edu/) - Updated combiner
- [Metabonaut](https://github.com/rformassspectrometry/Metabonaut) - Metabolomics explorer
- [Metaboprep](https://github.com/MRCIEU/metaboprep) - Metabolomics preparation
- [metaboprep v2](https://github.com/MRCIEU/metaboprep) - Updated preparation
- [Metanorm](https://github.com/UGent-LIMET/Metanorm) - Normalization
- [MOCCA](https://github.com/Bayer-Group/mocca) - Multivariate curve resolution
- [MPACT](https://github.com/BalunasLab/mpact) - Metabolomics processing
- [MS1FA](https://github.com/RuibingS/MS1FA_RShiny_dashboard) - MS1 feature analysis
- [MS2Planner](https://github.com/mohimanilab/MS2Planner) - MS/MS planning
- [MSMCE](https://github.com/WoodFY/MSMCE) - MS molecular connectivity
- [MS-TDF Software](https://github.com/CaixiaYuan/MS-TDF) - TIMS data filter
- [msFeaST](https://github.com/kevinmildau/msFeaST) - MS feature selection
- [MultiABLER](https://github.com/holab-hku/MultiABLER) - Multi-batch learning
- [mzLearn](https://github.com/ReviveMed/mzEmbed) - Parameter-free MS1 signal detection and alignment
- [mzRAPP](https://github.com/YasinEl/mzRAPP) - Reliability assessment of automated LC-MS data pre-processing
- [NeatMS](https://github.com/bihealth/NeatMS) - Neural network for MS
- [Nextflow4MS-DIAL](https://github.com/Nextflow4Metabolomics/nextflow4ms-dial) - Nextflow for MS-DIAL
- [NP-PRESS](https://github.com/MicroResearchLab/NP-PRESS) - Natural products processing
- [OpenNAU](https://github.com/zjuRong/openNAU) - Open NAU
- [PCPFM](https://github.com/shuzhao-li-lab/PythonCentricPipelineForMetabolomics) - Python pipeline
- [Peak Pair Pruner](https://github.com/QibinZhangLab/Peak_Pair_Pruner) - Peak pair analysis
- [PeakBot](https://github.com/christophuv/PeakBot) - Peak detection bot
- [PeakDecoder](https://github.com/EMSL-Computing/PeakDecoder) - Peak decoding
- [PeakDetective](https://github.com/pattilab/PeakDetective) - Peak investigation
- [PeakPerformance](https://peak-performance.readthedocs.io/en/latest/) - Peak performance
- [PFΔScreen](https://github.com/JonZwe/PFAScreen) - PFAS screening
- [QuanFormer](https://github.com/LinShuhaiLAB/QuanFormer) - Quantification transformer
- [RegFilter](https://github.com/computational-chemical-biology/regression_filter) - Regression filtering
- [rtmsEcho](https://github.com/nathanieltwarog/rtmsEcho) - RT-MS echo
- [Spectral Denoising](https://github.com/FangYuan717/spectral_denoising) - Denoising
- [TopNEXt](https://github.com/glasgowcompbio/vimms) - Top-N optimization
- [xcms](https://github.com/sneumann/xcms) - Peak detection and alignment

---

## ✅ Quality Control

- [ALISTER](https://itmp.shinyapps.io/alister/) - Pre-analytical quality guidance tool for assessing sample handling in lipidomics/metabolomics
- [CordBat](https://github.com/BorchLab/CordBat) - Batch effect correction for large-scale metabolomics spanning extended collection periods
- [Dbnorm](https://github.com/NBDZ/dbnorm) - R package for batch effect normalization in MS-based metabolomics across multiple batches
- [hRUV](https://github.com/SydneyBioX/hRUV) - Hierarchical method for removing unwanted variation in large-scale LC-MS metabolomics
- [iDIA-QC](https://github.com/guomics-lab/iDIA-QC) - Quality control tool using DIA-based consensus metrics to monitor LC-MS/MS performance
- [InjectionDesign](https://github.com/CSi-Studio/InjectionDesign) - Automated plate design for sample injection order with batch balancing and randomization
- [MAFFIN](https://github.com/HuanLab/MAFFIN) - Post-acquisition normalization removing sample amount variation in comparative metabolomics
- [marr](http://bioconductor.org/packages/release/bioc/html/marr.html) - Maximum rank reproducibility
- [MassQLab](https://github.com/JohnsonDylan/MassQLab) - Automated system suitability testing for mass spectrometry quality control
- [MatrixQCvis](https://bioconductor.org/packages/MatrixQCvis/) - Shiny app for interactive quality assessment and exploratory analysis of omics data
- [MetaMOPE](https://metamope.cmdm.tw) - Web tool for automated mobile phase optimization in HILIC LC-MS metabolomics
- [MetaPro](https://github.com/CSi-Studio/MetaPro) - Software for semi-targeted LC-MS with optimized quantification and spectral curation
- [MSNORM](https://pubs.acs.org/doi/suppl/10.1021/jasms.3c00295/suppl_file/js3c00295_si_001.txt) - R tool for evaluating and applying normalization methods to reduce systematic MS variability
- [Msquality](https://bioconductor.org/packages/MsQuality) - MS quality metrics
- [MultiBaC](https://github.com/ConesaLab/MultiBaC) - Batch effect correction for multi-omic datasets combining data from different batches
- [mzQuality](https://github.com/hankemeierlab/mzQuality) - Quality assurance tool for detecting experimental variation in MS-based metabolomics
- [OSCA Finder](https://github.com/dawnmengsjtu/OSCA-Finder) - Pipeline for osmolality-based calibration and deep learning-aided urine metabolomics
- [Paramounter](https://github.com/HuanLab/Paramounter) - Tool for optimizing LC-MS data processing parameters in untargeted metabolomics
- [PeakQC](https://github.com/pnnl/IonToolPack) - Automated quality control for monitoring peak integrity across large MS-based cohorts
- [QC4Metabolomics](https://github.com/stanstrup/QC4Metabolomics) - QC for metabolomics
- [QComics](https://github.com/ricoderks/QComics) - Standardized quality control protocol for ensuring reproducibility of metabolomics data
- [QuantyFey](https://github.com/CDLMarkus/QuantyFey) - QC-based signal drift correction in quantitative LC-MS analysis without internal standards
- [RALPS](https://github.com/zamboni-lab/RALPS) - Adversarial deep learning for inter-batch correction in untargeted MS metabolomics
- [Rapid QC-MS](https://github.com/czbiohub-sf/Rapid-QC-MS) - Rapid quality control
- [RawBeans](https://bitbucket.org/incpm/prot-qc/downloads/) - Tool for monitoring LC-MS system performance via raw data quality control
- [RawHummus](https://github.com/YonghuiDong/RawHummus) - Tool for inspecting raw LC-MS data quality including sensitivity drift and RT shifts
- [Shinyscreen](https://gitlab.com/uniluxembourg/lcsb/eci/shinyscreen) - R/Shiny app for exploration, visualization, and quality assessment of HR-MS raw data
- [WaveICA 2.0](https://github.com/dengkuistat/WaveICA_2.0) - Wavelet-based ICA method for batch effect removal in untargeted metabolomics

---

## ⏱️ RT (Retention Time)

- [ABCoRT](https://github.com/RiverCCC/ABCoRT) - Atom-bond co-learning framework for LC retention time prediction aiding metabolite identification
- [AsRTNet](https://github.com/piscookie/AS-RT) - Multimodal feature fusion model predicting retention times of arsenic compounds
- [CMM-RT](https://github.com/constantino-garcia/cmmrt) - ML regressors trained on 80,000 retention times for cross-method metabolite annotation
- [DeepGCN-RT](https://github.com/kangqiyue/DeepGCN-RT) - Deep graph convolutional network for accurate LC-MS retention time prediction
- [DeepRTAlign](https://github.com/PHOENIXcenter/deeprtalign) - Deep learning method for RT alignment in large-cohort LC-MS studies
- [GNN-RT](https://github.com/Qiong-Yang/GNN-RT) - Graph neural network predicting LC retention times directly from molecular structures
- [Graphormer-RT](https://github.com/HopkinsLaboratory/Graphormer-RT) - Graphormer model enabling RT prediction and comparison across different LC methods
- [QGeoGNN](https://github.com/woshixuhao/Retention-Time-Prediction-for-Chromatographic-Enantioseparation/tree/main) - ML framework predicting chiral molecule retention times to facilitate enantioseparation
- [ReTimeML](https://mikeallwright23-retime-app-lipid3-021zpv.streamlit.app/) - ML tool automating ceramide and sphingomyelin identification using LC-MS/MS retention times
- [retention_time_gnn](https://github.com/seokhokang/retention_time_gnn) - Transfer learning GNN predicting retention times for target chromatographic systems
- [RI-based-CPSC](https://github.com/WHU-Fenglab/RI-based-CPSC) - Retention index-based peak shift correction for peak alignment in untargeted metabolomics
- [ROASMI](https://github.com/FangYuan717/ROASMI) - Model for reliable retention order prediction to aid small molecule identification
- [RT-Ensemble Pred](https://gitlab.com/mikic93/rt-ensemble-pred) - Ensemble method predicting retention times across different chromatographic methods
- [RT-Pred](https://rtpred.ca/) - ML tool predicting HPLC retention times from molecular structures
- [RT-Transformer](https://github.com/01dadada/RT-Transformer) - Transformer model for scalable retention time prediction across diverse chromatographic conditions

---

## 🦠 Single Cell Metabolomics

- [MetaPhenotype](https://github.com/songyuan93/MetaPhenotype) - Meta-learning framework for single-cell MS phenotype classification with limited samples
- [MMEASE](https://idrblab.org/mmease/) - Comprehensive workflow for single-cell metabolomics data analysis
- [scFUMES](https://github.com/ChengF-Lab/scFUMES) - Algorithm integrating single-cell transcriptomics to identify cell-type-specific metabolic regulators
- [SCMeTA](https://github.com/ChenLab/SCMeTA) - Modular Python library for standardized single-cell metabolomics data processing

---

## 🗺️ Spatial Metabolomics

- [scSpaMet](https://github.com/coskunlab/ScSpaMet) - Framework for joint protein-metabolite profiling of single cells in tissues
- [SMART](https://github.com/bioinfo-ibms-pumc/SMART) - Open-source platform for precise molecular formula assignment in MSI
- [SMAnalyst](https://github.com/mzlab-research/SManalyst) - Integrated web platform for QC, preprocessing, identification, and analysis of spatial metabolomics
- [SMQVP](https://github.com/mzlab-research/SMQVP) - GUI software for spatial metabolomics data quality visualization, control, and preprocessing
- [SpaMTP](https://github.com/GenomicsMachineLearning/SpaMTP) - End-to-end analysis software linking spatial metabolomics with other spatial-omics modalities

---

## 🧩 Specialized

- [Amanida](https://github.com/mariallr/amanida) - R package for metabolomics meta-analysis combining p-values and fold-changes weighted by study size
- [arcMS](https://github.com/leesulab/arcMS) - R package for collecting MSE data from Waters UNIFI and storing in Apache Parquet format
- [ATLASx](https://lcsb-databases.epfl.ch/Atlas2) - Biochemical reaction atlas linking known metabolism to predicted enzymatic reactions
- [BAGO](https://github.com/huaxuyu/bago) - Bayesian optimization for autonomous LC gradient optimization using limited data
- [BioPKS Pipeline, RetroTide](https://github.com/JBEI/BioPKS-Pipeline) - Automated retrobiosynthesis tool integrating polyketide synthases for pathway design
- [BioTransformer 3.0](https://biotransformer.ca) - Biotransformation prediction
- [BitterMasS](https://github.com/Niv-Lab/BitterPredict1) - ML model predicting bitter compounds from tandem mass spectrometry data
- [BreathXplorer](https://github.com/HuanLab/breathXplorer) - Python package for processing real-time HRMS-based exhaled breath metabolomics data
- [CCWeights](https://github.com/YonghuiDong/CCWeights) - R package for automated calibration curve weighting factor selection for LC-MS quantification
- [ChemEcho](https://github.com/biorack/chemecho) - Tool translating fragmentation pattern knowledge into human-readable rules via MassQL queries
- [ChemSpectra](https://github.com/ComPlat/chem-spectra-app) - Web-based software for visualizing and analyzing IR, MS, and NMR spectroscopic data
- [COMMIT](https://github.com/pwendering/COMMIT/tree/v1.0.0) - Gap-filling approach for microbial community metabolic models considering metabolite permeability
- [COVRECON](https://bitbucket.org/mosys-univie/covrecon) - ML framework identifying metabolite biomarkers of fitness from metabolomics cohort data
- [DarkNPS](https://github.com/skinniderlab/DarkNPS) - Deep learning tool for identifying emerging new psychoactive substances from MS data
- [DeepMet](https://deepmet.org/) - Large language model approach for interpreting and discovering novel mammalian metabolites
- [Easy-Amanida](https://github.com/mariallr/easy-amanida) - Shiny app simplifying metabolomics meta-analysis using reported p-values and fold-changes
- [ECDFormer](https://github.com/HowardLi1984/ECDFormer) - Transformer model for accurate electronic circular dichroism spectrum prediction
- [FORUM](https://github.com/eMetaboHUB/Forum-DiseasesChem) - Knowledge graph linking chemicals and biomedical concepts for metabolic signature interpretation
- [GlyKAn AZ](https://github.com/ZaylaSchaeffer/GlyKAn-AZ-application) - Open-source GUI for automated glycan structure identification from LC-MS/MS data
- [homologueDiscoverer](https://github.com/kevinmildau/homologueDiscoverer) - R package for detecting homologous polymer series in untargeted metabolomics data
- [IonFlow](https://github.com/wanchanglin/ionflow) - Galaxy/R tool for ionomics data analysis using multivariate statistics and network analysis
- [M2S](https://github.com/rjdossan/M2S) - Tool for aligning and matching untargeted LC-MS features across separate datasets
- [Maldi Transformer](https://github.com/gdewael/maldi-nn) - Transformer neural network architecture for MALDI-TOF MS data analysis
- [MassQL](https://github.com/mwang87/MassQueryLanguage) - Query language for directly querying raw mass spectrometry data with user-defined filters
- [MEISTER](https://github.com/richardxie1119/MEISTER) - Integrative MS framework enabling brain-wide 3D biochemical mapping with single-cell resolution
- [MetaboliteChat](https://github.com/13501274828/MetaboliteChat) - Multimodal LLM integrating molecular-graph and image reasoning for metabolite analysis
- [MetaboLM](https://github.com/Qiu-Shizheng/MetaboLM) - Transformer language model pre-trained on plasma metabolomics for chronic disease prediction
- [Metaboverse](https://github.com/Metaboverse/Metaboverse) - Tool for contextualized exploration of metabolic network data and hypothesis generation
- [MetCohort](https://github.com/JunYang2021/MetCohort) - Large-scale LC-HRMS metabolomics raw data alignment, feature detection, and quantification
- [MetDIT](https://github.com/Li-OmicsLab/MetDIT) - Deep CNN framework transforming 1D metabolomics data into 2D images for clinical analysis
- [mGWAS-Explorer 2.0](https://www.mgwas.ca/) - Platform for causal analysis between metabolites and phenotypes in metabolomics GWAS
- [MINE 2.0](https://minedatabase.ci.northwestern.edu) - In silico metabolic network expansion database for annotating unknown metabolomics peaks
- [MODAPro](https://github.com/zhaoxiaoqi0714/MODAPro) - Deep learning framework integrating variational graph autoencoders for multi-omics integration
- [MolSpectLLM](https://github.com/Eurekashen/MolSpectLLM) - Foundation model integrating SMILES, experimental spectra, and 3D structure for molecular analysis
- [MSConnect](https://github.com/RTKlab-BYU/MSConnect) - Integrated data management platform enabling traceable multi-software MS analysis workflows
- [MSident](https://github.com/cplqam/GUI) - GUI tool automating metabolite identification from ROIMCR chemometric analysis
- [MS-RT](https://github.com/XianghuWang-287/Metabolomics_Clustering_Benchmark) - Benchmark for evaluating MS/MS spectral clustering algorithms in metabolomics
- [MS2toImg](https://github.com/jsehhs/MS2toImg-CNN) - CNN approach converting MS/MS spectra to images for identification
- [Multipass CCS Refiner](https://ericgier.shinyapps.io/Multipass-CCS-Refiner/) - Shiny app for calibrating multipass CCS measurements from cyclic ion mobility
- [NP Analyst](https://www.npanalyst.org/) - Open platform mapping bioassay activity to MS features for natural product discovery
- [NPFimg](http://github.com/poomcj/NPFimg) - Image processing and ML method identifying chemo/biomarker features in chromatography-MS
- [Pickaxe](https://github.com/tyo-nu/MINE-Database) - Flexible biochemical reaction prediction tool for generating enzymatic reaction networks
- [ptairMS](https://bioconductor.org/packages/release/bioc/html/ptairMS.html) - Bioconductor suite for processing PTR-TOF-MS exhaled breath data in cohort studies
- [SMITER](https://github.com/LeidelLab/SMITER) - Python tool simulating LC-MS/MS runs for any biomolecule using formula-based calculations
- [SPIFFED](https://github.com/bio-it-station/SPIFFED) - Computational method for analyzing co-fractionation MS data to construct protein interaction networks
- [Umatrix](https://cran.r-project.org/package=Umatrix) - Generative algorithm using self-organizing maps to augment small biomedical sample sizes
- [VIMMS 2.0](https://github.com/glasgowcompbio/vimms/) - Framework for developing and testing MS fragmentation strategies in silico and on instruments

---

## 📕 Spectral Library

- [CIeaD](https://www.moleculardetective.org/Links.html) - Open-access plant metabolite spectral library with complementary CID and EAD fragmentation spectra
- [DNA Adduct Portal](https://sites.google.com/umn.edu/dnaadductportal) - Curated database and portal for DNA adduct mass spectrometric detection and characterization
- [EMBL-MCF 2.0](https://www.embl.org/groups/metabolomics/instrumentation-and-software/#MCF-library) - Updated LC-MS metabolite spectral library with authenticated reference standards
- [GNPS Drug Library](https://github.com/ninahaoqizhao/Manuscript_GNPS_Drug_Library) - MS/MS reference library of drugs and metabolites for screening drug exposure
- [HighResNPS](https://highresnps.com/) - HRMS spectral database for suspect screening of new psychoactive substances
- [In silico infrared spectral library of human metabolites](https://zenodo.org/records/7706021) - IR library
- [MassBank](https://massbank.jp/) - Mass spectral database
- [metScribeR](https://github.com/ncats/metScribeR) - R/Shiny package accelerating creation of in-house retention time metabolite standard libraries
- [mFAM](https://github.com/MassBank/MassBank-data/commit/a91b1ca4841aea536545f7c1d452c1f80d225e84) - Approach for evaluating MS features at the metabolite family level for functional annotation
- [MIADB](https://gnps.ucsd.edu) - Public MS/MS spectral database of monoterpene indole alkaloids on the GNPS platform
- [MS2extrcat](https://www.cooperstonelab.com/MS2extract/) - R package facilitating automated creation of MS/MS reference spectral libraries
- [PASL](https://gnps.ucsd.edu/ProteoSAFe/gnpslibrary.jsp?library=PYRROLIZIDINE-ALKALOID-SPECTRAL-LIBRARY) - High-resolution Orbitrap MS/MS spectral library of pyrrolizidine alkaloids
- [RMB-mix](https://github.com/rmassbank/rmbmix) - Tool generating spectral libraries from complex mixtures to expand exposomics coverage
- [UCL-MetIsoLib](https://github.com/kserafimov10/UCLMetIsoLib) - HILIC-based isomer-resolved HRMS/MS metabolite library with 334 annotated metabolites
- [WFSR Food Safety Mass Spectral Library](https://www.wur.nl/en/show/food-safety-mass-spectral-library.htm) - Food safety library

---

## 📈 Statistical

- [cluster-CV](https://github.com/cluster-cv) - Integrative approach combining univariate, multivariate, and pathway analysis for metabolomics
- [DisCo P-ad](https://github.com/KechrisLab/DisCoPad) - Multiple testing correction method for controlling false discoveries in metabolomics
- [GetFeatistics](https://github.com/FrigerioGianfranco/GetFeatistics) - R package for calibration-based quantification and statistical analysis of metabolomics data
- [imputomics](https://github.com/BioGenies/imputomics) - R framework benchmarking and applying 52 missing value imputation algorithms for metabolomics
- [MAI](https://bioconductor.org/packages/release/bioc/html/MAI.html) - Multiple annotation imputation
- [MAMSI](https://pypi.org/project/mamsi/) - Workflow integrating multi-assay untargeted LC-MS metabolomics for cross-assay biomarker discovery
- [MetaboLINK/ PCA-GLASSO](https://github.com/jlichtarge/pcaGLASSO) - PCA combined with Graphical Lasso for dimensionality reduction of large metabolomics datasets
- [MetaHD](https://bookdown.org/a2delivera/MetaHD/) - Multivariate meta-analysis model handling correlations and missing values in high-dimensional metabolomics
- [MeTEor](https://github.com/scibiome/meteor) - Interactive Shiny app for explorative analysis and visualization of longitudinal metabolomics data
- [MultiClassMetabo](http://idrblab.cn/multiclassmetabo/) - Platform for multiclass metabolomics classification and biomarker identification
- [NOREVA](https://idrblab.org/noreva/) - Protocol for evaluating and selecting optimal metabolomics data processing workflows
- [omicsMIC](https://github.com/WQLin8/omicsMIC) - Comprehensive tool comparing missing value imputation methods for MS-based omics data
- [PLSKO](https://github.com/guannan-yang/PLSKO) - Assumption-free knockoff generator for FDR-controlled differential expression in omics data
- [PredCMB](https://www.sysbiolab.org/predcmb) - Tool predicting changes in individual metabolite levels from shotgun metagenome data
- [rMisbeta](https://cran.r-project.org/package=rMisbeta) - R package for robust imputation of missing values in the presence of outliers
- [SMART 2.0](https://github.com/YuJenL/SMART) - Integrated tool for targeted/untargeted metabolomics with normalization, QC, and pathway analysis
- [SpectraX](https://github.com/weantony/SpectraX) - MATLAB app for automated analysis of complex direct/ambient mass spectrometry spectra
- [wscaling.R](https://github.com/nishithkumarpaul/robustScaling/blob/main/wscaling.R) - R script providing outlier-robust weighted scaling methods for metabolomics preprocessing

---

## 🎯 Targeted

- [automRm](https://gitlab.gwdg.de/joerg.buescher/automrm) - R package for fully automatic MRM LC-MS data preprocessing using ML peak detection
- [CLAW-MRM](https://github.com/chopralab/CLAW) - Automated platform for lipid annotation, statistical analysis, and MRM data parsing
- [MRMPro](http://mrmpro.csibio.com/) - Open-source tool for efficient MRM peak review and correction in large cohort studies
- [MRMQuant](https://github.com/kslynn128171/MRMQuant) - Automated tool for accurate targeted metabolomic quantitation from MRM LC-MS/MS
- [Norm ISWSVR](https://github.com/Dingxian666/Norm-ISWSVR) - Normalization method correcting batch effects and analytical drift in targeted metabolomics
- [RHermes](https://github.com/RogerGinBer/RHermes) - Formula-oriented, peak-detection-free method optimizing MS2 acquisition from raw LC/MS1 data
- [SCALiR](https://lewisresearchgroup-ms-conc-streamlit-app-w9e64f.streamlit.app/) - Tool for automated metabolite concentration calculation from LC-MS standard curves
- [TARDIS](https://github.com/UGent-LIMET/TARDIS) - Open-source R package for automated preprocessing of targeted LC-MS metabolomics/lipidomics data

---

## 📊 Visualization

- [ClusterApp](http://ccbl-apps.fcfrp.usp.br/ClusterApp) - Web app for interactive clustering and dimensionality reduction of metabolomics data
- [EDaViS](https://github.com/LRuehrmund/RoMBAT) - Script for visualizing complex time-resolved volatile metabolomics data
- [GraphBio](https://github.com/databio2022/GraphBio) - Shiny web app providing 15 visualization methods for omics data without programming
- [lcmsWorld](https://github.com/PGB-LIV/lcmsWorld) - 3D visualization tool for quality-control inspection of raw LC-MS data
- [MODE](https://github.com/pmartR/MODE_ShinyApp) - Interactive Shiny app for detailed exploration of multidimensional omics data
- [MS-VIS](https://github.com/hkevgill/MS-VIS) - R Shiny tool for visualizing mass spectra and analyzing mass lists
- [PathwayNexus](https://www.cls.uni-konstanz.de/software/pathway-nexus) - Tool for mapping metabolomics data across multiple conditions onto metabolic pathways
- [RAIChU](https://github.com/BTheDragonMaster/RAIChU) - Tool for predicting and visualizing NRPS and PKS biosynthetic pathway products
- [RDD metabolomics platform](https://gnps-rdd.streamlit.app/) - Platform contextualizing unannotated MS/MS spectra by comparison to reference datasets
- [SpecXplore](https://github.com/kevinmildau/specXplore) - Exploratory tool for manual annotation and visualization of mass spectral data

---

## 📄 License

This project is licensed under the terms in the [LICENSE](LICENSE) file.
