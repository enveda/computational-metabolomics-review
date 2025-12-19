# Trends in Computational Metabolomics (2021–2025)
[![Awesome](https://cdn.jsdelivr.net/gh/sindresorhus/awesome@d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)
A curated list of computational tools, databases, and software for metabolomics research published between 2021-2025.

## Citation

> Domingo-Fernández, D., Healy, D., Kind, T., Allen, A., Colluru, V., and Misra, B. (2025). Trends in computational metabolomics in the last five years (2021–2025). chemRxiv.

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
- [Mass Spectrum Transformer](https://github.com/chensaian/TransG-Net) - Transformer-based mass spectrum analysis
- [MWFormer](https://github.com/zhanghailiangcsu/MWFormer) - Molecular weight prediction
- [Spec2Class](https://huggingface.co/VickiPol/binary_models) - Spectrum to class prediction

### 🔍 De Novo Generation

- [DiffMS](https://github.com/coleygroup/DiffMS) - Diffusion model for MS
- [MASSISTANT](https://github.com/BartaLazar/MASSISTANT/) - Mass spectrometry assistant
- [MS-BART](https://github.com/OpenDFM/MS-BART) - BART model for mass spectrometry
- [MS2Mol](https://doi.org/10.26434/chemrxiv-2023-vsmpx-v3) - MS/MS to molecule prediction
- [MSGo](http://github.com/aaronma2020/MSGO) - Mass spectrometry structure generation
- [MSNovelist](https://github.com/meowcat/MSNovelist) - De novo structure elucidation
- [Mass2SMILES](https://github.com/volvox292/mass2smiles) - Mass to SMILES conversion
- [MassGenie](https://github.com/neilswainston/FragGenie) - Fragment generation
- [OMG](https://github.com/HassounLab/OMGG) - Optimal molecule generation
- [SEISMiQ](https://github.com/Boehringer-Ingelheim/seismiq) - Structure elucidation
- [Spec2Mol](https://github.com/KavrakiLab/Spec2Mol) - Spectrum to molecule
- [TeFT](https://github.com/thumingo/TeFT) - Transformer for fragmentation

### 🧠 Learned Spectrum Representations

- [ChemEmbed](https://github.com/massspecdl/ChemEmbed) - Chemical embeddings
- [CLERMS](https://github.com/HaldamirS/CLERMS) - Contrastive learning for MS
- [DeepMASS](https://github.com/hcji/DeepMASS2_Data_Processing) - Deep learning for mass spectrometry
- [DreaMS](https://github.com/pluskal-lab/DreaMS) - Deep learning representations
- [LSM1-MS2](https://github.com/matterworksbio/LSM1-MS2) - Language-spectrum model
- [MS2DeepScore](https://github.com/matchms/ms2deepscore) - Deep learning similarity scoring
- [MS2DeepScore 2.0](https://github.com/matchms/MS2DeepScore) - Updated deep learning similarity
- [MSBERT](https://github.com/zhanghailiangcsu/MSBERT) - BERT for mass spectrometry
- [NaFM](https://github.com/TomAIDD/NaFM-Official) - Neural attention foundation model
- [Spec2Vec](https://github.com/iomega/spec2vec) - Word2Vec for spectra
- [SpecEmbedding](https://huggingface.co/spaces/xp113280/SpecEmbedding) - Spectrum embedding

### 🏎️ Molecular Formula Prediction

- [CRB-FCC](https://github.com/iconSS/FCC/releases/tag/v1) - Fragmental chain characterization
- [FIDDLE](https://github.com/JosieHong/FIDDLE) - Formula identification
- [FSA](https://github.com/kslynn128171/FSA) - Formula subset analysis
- [IDSL.UFA](https://cran.r-project.org/package=IDSL.UFA) - United formula annotation
- [MIST-CF](https://github.com/samgoldman97/mist-cf) - Chemical formula prediction
- [RASSP](https://github.com/thejonaslab/rassp-public) - Rapid annotation

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

- [3D-MPEA](https://github.com/ZibianFan/3D-MPEA) - 3D molecular property embedding
- [BAM](https://github.com/HassounLab/BAM) - Biotransformation-based annotation method
- [ChemWalker](https://github.com/computational-chemical-biology/ChemWalker) - Annotation propagation via random walks
- [ConCISE](http://github.com/zquinlan/concise) - Consensus classifications
- [E-SGMN](https://sourceforge.net/projects/e-sgmn/) - Enhanced structure-guided molecular networking
- [IIMN](https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/) - Ion identity molecular networking
- [IMN4NPD](https://github.com/mwang87/NP-Classifier) - Integrated molecular networking for natural products
- [MCN](https://github.com/Alexander0/molecular_communities) - Molecular community networking
- [MMSA](https://multiplealignment.gnps2.org/setscreation) - Multiple mass spectral alignment
- [ModiFinder](https://github.com/Wang-Bioinformatics-Lab/ModiFinder_base) - Modification site localization
- [MolNotator](https://github.com/ZzakB/MolNotator) - Mass spectral feature to molecule workflow
- [MS-Net](https://zenodo.org/records/17669288) - MS networking
- [SGMNS](https://github.com/DLUT-datas/SGMNS) - Structure-guided molecular network strategy
- [SIMILE](https://github.com/biorack/simile) - Spectral alignment with statistical significance

### Others

- [AnnoMe](https://github.com/chrboku/AnnoMe) - MS/MS spectra classification
- [AnnoSM](https://github.com/YangHuaLab/AnnoSM) - Substituent mode annotation
- [BioTransformer4.0](https://github.com/Wishartlab-openscience/Biotransformer) - Biotransformation prediction
- [IDSL.CSA](https://cran.r-project.org/package=IDSL.CSA) - Composite spectra analysis
- [Inventa](https://luigiquiros.github.io/inventa/) - Structural novelty discovery
- [ipaPy2](https://github.com/francescodc87/ipaPy2) - Integrated probabilistic annotation
- [MADGEN](https://github.com/HassounLab/MADGEN) - Metabolite annotation
- [MCheM](https://github.com/sirius-ms/sirius) - Multiplexed chemical metabolomics
- [MetaboAnnotatoR](https://github.com/gggraca/MetaboAnnotatoR) - All-ion fragmentation annotation
- [MMST](https://github.com/mpriessner/MultiModalSpectralTransformer) - Multimodal spectral transformer
- [MS2DECIDE](https://github.com/MejriY/MS2DECIDE) - Decision theory for annotation
- [mWISE](https://github.com/b2slab/mWISE) - Metabolite identification
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
- [DeepCDM](https://github.com/ADNLab-SCU/DeepCDMs) - Deep collision dissociation model
- [ESP](https://github.com/HassounLab/ESP) - Enumerated substructure prediction
- [FIORA](https://github.com/BAMeScience/fiora) - Fragmentation prediction
- [fragnnet](https://github.com/FraGNNet/fragnnet) - Fragmentation neural network
- [GrAFF-MS](https://www.envedabio.com/posts/a-new-scalable-ml-model-for-accurate-prediction-of-small-molecule-mass-spectra) - Graph attention for fragmentation
- [HDSE-MS](https://github.com/lzjforyou/HDSE-MS) - High-dimensional spectral encoding
- [ICEBERG](https://github.com/coleygroup/ms-pred) - In silico chemical biology
- [MARASON](https://github.com/coleygroup/ms-pred) - Mass spectrum prediction
- [MassFormer](https://github.com/Roestlab/massformer/) - Transformer for mass spectra
- [MassKG](https://xomics.com.cn/masskg) - Mass spectrometry knowledge graph
- [MS2Compound](https://github.com/beherasan/MS2Compound) - Compound identification
- [ms-pred](https://github.com/coleygroup/ms-pred) - Mass spectrum prediction
- [NPS-MS](https://nps-ms.ca/) - Novel psychoactive substances MS
- [PPGB-MS2](https://github.com/zhengfj1994/PPGB_MS2) - Pathway-guided prediction
- [QC-GN2oMS2](https://github.com/PNNL-m-q/qcgnoms) - Quantum chemistry MS prediction
- [SingleFrag](https://github.com/MaribelPR/SingleFrag) - Single fragmentation prediction

### 🧰 Workflow Tools

- [ENTAiLS Toolkit](https://whitehead-heather.github.io/ENTAiLSToolkit/) - Annotation toolkit
- [MetaboCoreUtils, MetaboAnnotation and CompoundDb](https://github.com/jorainer/MetaboAnnotationTutorials) - R packages for annotation
- [MetaFetcheR](https://github.com/komorowskilab/MetaFetcheR/) - Metabolite fetching
- [MS2Query](https://github.com/iomega/ms2query) - Query-based annotation
- [mscompiler](https://github.com/QizhiSu/mspcompiler) - MS compiler
- [mspepsearchr](https://github.com/AndreySamokhin/mssearchr) - Peptide search
- [PS2MS](https://github.com/jhhung/PS2MS) - Protein to MS
- [RapidMass](https://github.com/Katherine00689/RapidMass) - Rapid mass analysis
- [Scannotation](https://github.com/scannotation/Scannotation_software) - Scan annotation
- [ShinyMetID](https://github.com/jjs3098/ShinyMetID) - Shiny app for metabolite ID
- [SiMD](https://si-simd.com/) - Similarity-based metabolite discovery

---

## 🧪 Benchmark/Dataset

- [MassSpecGym](https://github.com/pluskal-lab/MassSpecGym) - Mass spectrometry benchmark
- [MetaBench](https://github.com/metabench/metabench) - Metabolomics benchmarking

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

- [AMDB](https://amdb.online/) - Antimicrobial metabolite database
- [BinDiscover database](https://bindiscover.metabolomics.us/) - Binary discovery database
- [CCDB](https://github.com/idslme/CCDB/tree/main) - Compound class database
- [COCONUT](https://coconut.naturalproducts.net/) - Natural products database
- [DNA adduct database](https://gitlab.nexs-metabolomics/projects/dna_adductomics_database) - DNA adducts
- [EnzyMine](http://www.rxnfinder.org/enzymine/) - Enzyme mining
- [foodMASST](https://masst.ucsd.edu/foodmasst) - Food mass spectrometry
- [FragHub](https://github.com/eMetaboHUB/FragHub) - Fragment database hub
- [GNPS Dashboard](https://dashboard.gnps2.org/) - GNPS data visualization
- [HMDB 5.0](https://hmdb.ca) - Human metabolome database
- [HREI-MSDB](https://github.com/mtshn/gchrmsexplain) - High-resolution EI-MS database
- [Human Hair Atlas](https://metabolomics.cloud/hair/) - Hair metabolome atlas
- [LEAFBot](https://ccms-ucsd.github.io/GNPSDocumentation/gnpslibraries/) - Leaf metabolomics
- [MarkerDB 2.0](https://markerdb.ca/) - Biomarker database
- [MassBase](http://webs2.kazusa.or.jp/massbase/) - Mass spectrometry database
- [MassSpecBlocks](https://github.com/privrja/MassSpecBlocks) - Building blocks for MS
- [MCID database](http://www.mycompoundid.org/) - Compound identification
- [MedMeta](https://cbcb.cdutcm.edu.cn/MedMeta) - Medical metabolomics
- [MetaboLights](https://www.ebi.ac.uk/metabolights) - Metabolomics data repository
- [Metabolome atlas of the aging mouse brain](https://mouse.atlas.metabolomics.us/) - Mouse brain metabolome
- [MetalinksDB](https://metalinks.omnipathdb.org/) - Metabolite links database
- [MetaNetX 2025 update](https://www.metanetx.org/) - Metabolic network database
- [MetHoS](https://methos.cebitec.uni-bielefeld.de/welcomeview) - Metabolomics hosting
- [METLIN-CCS](https://metlin.scripps.edu/) - METLIN with CCS values
- [microbeMASST](https://github.com/robinschmid/microbe_masst) - Microbial mass spectrometry
- [MiMeDB](https://mimedb.org/) - Microbiome metabolome database
- [MiMeDB 2.0](https://mimedb.org/) - Updated microbiome metabolome database
- [MSCAT](https://mscat.metabolomicsworkbench.org/) - MS catalog
- [MSnLib](https://github.com/corinnabrungs/msn_tree_library) - MSn library
- [Natural Products Atlas 2.0](https://www.npatlas.org/) - Natural products database
- [NGlycDB](https://metaspace2020.org/project/velickovic_Nglyc_2021?tab=datasets) - N-glycan database
- [NIST23](https://chemdata.nist.gov/dokuwiki/doku.php?id=chemdata:nist23-msms) - NIST mass spectral library
- [NMRlipids Databank](https://nmrlipids.github.io/) - NMR lipids database
- [NPMine](https://github.com/computational-chemical-biology/npmine) - Natural products mining
- [PharmMet DB](https://github.com/ClinicalBiomarkersLaborabory/PharmMet) - Pharmacometabolomics
- [plantMASST](https://github.com/helenamrusso/plantmasst) - Plant mass spectrometry
- [Pyrfume](https://pyrfume.org/) - Olfactory metabolomics
- [RepoRT](https://github.com/michaelwitting/RepoRT) - Retention time repository
- [Spectraverse](https://github.com/skinniderlab/spectraverse-analysis) - Spectral universe
- [TOMATOMET](http://metabolites.in/tomato‐fruits/) - Tomato metabolome

---

## 💊 Drug Discovery

- [BitBIRCH-Lean](https://github.com/mqcomplab/bblean) - Chemical space exploration
- [ChromaQuant](https://hplcfdu.shinyapps.io/ChromaQuant/) - Chromatographic quantification
- [DMetFinder](https://github.com/KeWHs/dmetdata) - Drug metabolite finder
- [Limelight](https://limelight-ms.org) - MS data sharing
- [SynFrag](https://synfrag.simm.ac.cn/) - Synthetic fragment analysis

---

## 🌱 Exposomics

- [CMDN](https://github.com/LeaveMeNotTonight/CMDN) - Chemical mixture detection network
- [EISA-EXPOSOME](https://github.com/Lab-XUE/EISA-EXPOSOME) - Exposome analysis
- [ExposomeX](http://www.exposomex.cn/) - Exposomics platform
- [FeatureHunter](https://msomics.abrc.sinica.edu.tw/FeatureHunter/) - Feature detection
- [FluoroMatch IM](https://fluoromatch.com) - Fluorinated compound matching
- [HalogenFinder](https://halogenfinder.com) - Halogenated compound detection
- [HExpMetDB](https://github.com/FangLabNTU/HExpMetDB) - Human exposome database
- [MDRB](https://github.com/933ZhangDD/MDRB) - Metabolite-disease relationship
- [MetabFlow](https://bddg.hznu.edu.cn/metabflow/) - Metabolomics workflow
- [MSThunder](https://github.com/LQZ0123/MSThunder.git) - MS data processing

---

## 🌀 FT-ICR MS

- [MetaboDirect](https://github.com/Coayala/MetaboDirect) - Direct analysis
- [MoleTrans](https://github.com/JibaoLiu/MoleTrans) - Molecular transformation

---

## 🗂️ Formats

- [Aird](https://github.com/CSi-Studio/AirdPro) - Aird format tools
- [AlphaTims](https://github.com/MannLabs/alphatims) - TIMS data handling
- [CloMet](https://github.com/rmallol/clomet) - Cloud metabolomics
- [Dear-OMG](https://github.com/jianweishuai/Dear-OMG) - OMG format tools
- [mspack](https://github.com/fhanau/mspack) - MS data packing
- [mwtab Python Library for RESTful Access](https://github.com/MoseleyBioinformaticsLab/mwtab) - mwTab access
- [MZA](https://github.com/PNNL-m-q/mza) - MZ archive format
- [mzapy](https://github.com/PNNL-m-q/mzapy) - Python mza tools
- [mzPeak](https://github.com/mobiusklein/mzpeak_prototyping) - Peak handling
- [pyOpenMS-viz](https://github.com/OpenMS/pyopenms_viz) - OpenMS visualization
- [rawrr R](https://github.com/fgcz/rawrr) - Raw file reading in R
- [SpectriPy](https://github.com/rformassspectrometry/SpectriPy) - Spectral data in Python
- [TIMSCONVERT](https://github.com/gtluu/timsconvert) - TIMS conversion

---

## 🔥 GC-MS

- [ADAP-KDB](https://adap.cloud/) - ADAP knowledge database
- [CRISP](https://github.com/vivekmathema/GCxGC-CRISP) - GCxGC data processing
- [eROI](https://github.com/Gaozhuang-big/eROI) - Enhanced ROI detection
- [gc-ims-tools](https://github.com/Charisma-Mannheim/gc-ims-tools) - GC-IMS tools
- [GcDUO](https://github.com/mariallr/GcDuo) - GC dual detection
- [GCMS-ID](https://gcms-id.ca) - GC-MS identification
- [isoSCAN](https://github.com/jcapelladesto/isoSCAN) - Isotope scanning
- [MACE](http://www.oc.tu-bs.de/schulz/html/MACE.html) - Mass spectral analysis
- [MetaboPAC](https://github.com/gtStyLab/MetaboPAC) - Metabolome pathway analysis
- [mzrtsim](https://github.com/yufree/mzrtsim) - Simulation tools
- [RIpred](https://ripred.ca/) - Retention index prediction
- [SERDA](https://github.com/slfan2013/SERDA) - Spectral entropy analysis
- [SpecTUS](https://github.com/hejjack/SpecTUS) - Spectral tools
- [uafR](github.org/castratton/uafR) - UAF R package
- [UniqPy](https://github.com/DNKonanov/uni_cli) - Unique peak detection

---

## 🔋 Holistic/Standalone Tools

- [FluoroMatch 2.0](http://innovativeomics.com/software/fluoromatch-flow-covers-entire-pfas-workflow/) - PFAS workflow
- [GNPS2](https://gnps2.org) - Global natural products social networking
- [MAVEN2](https://github.com/eugenemel/maven/releases/tag/2.11.8) - Metabolomics analysis
- [MAW](https://github.com/zmahnoor14/MAW) - Metabolomics analysis workflow
- [MassCube](https://github.com/huaxuyu/masscube) - Mass spectrometry cube
- [MeRgeION](https://github.com/daniellyz/meRgeION2) - Ion merging
- [MetaboReport](https://metaboreport.com/) - Metabolomics reporting
- [MetaboScape 2025](https://www.bruker.com/en/products-and-solutions/mass-spectrometry/ms-software/metaboscape.html) - Bruker metabolomics
- [Metabox 2.0](http://metsysbio.com/metabox) - Metabolomics toolbox
- [MetEx](https://mo.princeton.edu/MetEx/) - Metabolite extraction
- [PlantMetSuite](https://github.com/YaonanLab/PlantMetSuite) - Plant metabolomics suite
- [Workflow4Metabolomics](https://workflow4metabolomics.org/) - Galaxy workflows

---

## 🌈 IR

- [Chemprop-IR](https://github.com/gfm-collab/chemprop-IR) - IR property prediction
- [Graphormer-IR](https://github.com/HopkinsLaboratory/Graphormer-IR) - IR prediction with Graphormer
- [Spectro](https://github.com/ChemAI-Lab/spectro) - IR/NMR annotation

---

## 🖼️ Imaging MS

- [13C-SpaceM](https://github.com/Buglakova/13C-SpaceM) - 13C spatial metabolomics
- [Aird-MSI](https://github.com/CSi-Studio/mzmine3/tree/Aird2.0) - Aird for MSI
- [Cardinal v3](https://github.com/kuwisdelu/Cardinal) - MSI analysis in R
- [DeepION](https://github.com/BioNet-XMU/DeepION) - Deep learning for ion images
- [deepPseudoMSI](https://www.deeppseudomsi.org/) - Deep learning pseudo MSI
- [DIMPLE](https://github.com/dickinsonlab/DIMPLE-code) - DIMPLE imaging
- [HT SpaceM](https://github.com/delafior/HT-SpaceM) - High-throughput SpaceM
- [i2i](https://github.com/LanekoffLab/i2i) - Image to image
- [imzML Writer](https://github.com/VIU-Metabolomics/imzML_Writer) - imzML writing
- [LipidQMap](https://github.com/swinnenteam/LipidQMap) - Lipid quantitative mapping
- [M2aia](https://github.com/m2aia/m2aia) - MSI analysis platform
- [macroMS](https://github.com/macroMS/macroMS_web) - Macro mass spectrometry
- [mass2adduct](https://github.com/kbseah/mass2adduct) - Adduct detection
- [massNet](https://github.com/wabdelmoula/massNet) - Mass networking for MSI
- [MassVision](https://SlicerMassVision.readthedocs.io) - 3D MSI visualization
- [Met-ID](https://github.com/pbjarterot/Met-ID) - Metabolite identification for MSI
- [METASPACE-ML](https://github.com/metaspace2020/metaspace) - METASPACE machine learning
- [MetaVision3D](https://metavision3d.rc.ufl.edu) - 3D metabolomics vision
- [MSI-Explorer](https://github.com/MMV-Lab/MSI-Explorer) - MSI exploration
- [MSIannotator](https://github.com/Yingzhu96/MSIannotator) - MSI annotation
- [MSIGen](https://github.com/LabLaskin/MSIGen) - MSI data generation
- [MSight](https://github.com/laurenfields/MSIght) - MSI insight
- [MSIpixel](https://github.com/gmartano/MSIpixel) - Pixel-level MSI
- [msiFlow](https://github.com/Immunodynamics-Engel-Lab/msiflow) - MSI workflow
- [MSroi](http://www.mcrals.info/) - MSI ROI analysis
- [Multi‑MSI Processor](https://github.com/multimsiproc/MMP) - Multi-MSI processing
- [Pew2](https://github.com/djdt/pewpew) - MSI visualization
- [rMSIannotation](https://github.com/prafols/rMSIproc) - R MSI annotation
- [S2IsoMEr](https://github.com/alexandrovteam/S2IsoMEr) - Isotope analysis for MSI
- [SagMSI](https://github.com/BioNet-XMU/SagMSI) - Sagittal MSI
- [SmartGate](https://github.com/zhanglabtools/SmartGate) - Smart gating for MSI
- [SONAR-MSI](https://github.com/GrantBellic/SONAR-MSI) - SONAR for MSI
- [SpaceM](https://github.com/alexandrovteam/SpaceM) - Spatial metabolomics
- [SpatialMETA](https://github.com/WanluLiuLab/SpatialMETA) - Spatial meta-analysis
- [SSN](https://github.com/LanekoffLab/i2i) - Spectral similarity network
- [subMALDI](https://github.com/wesleyburr/subMaldi) - Sub-MALDI analysis

---

## 🌬️ Ion Mobility MS

- [AllCCS2](http://allccs.zhulab.cn/) - CCS prediction database
- [AutoCCS](https://github.com/PNNL-Comp-Mass-Spec/AutoCCS) - Automated CCS calculation
- [AutonoMS](https://autonoms.readthedocs.io) - Autonomous MS
- [CCS Predictor 2.0](https://github.com/facundof2016/CCSP2.0) - CCS prediction
- [CCSfind](https://github.com/sangeeta97/ccs_find) - CCS finding
- [GCIMS](https://github.com/sipss/GCIMS) - GC-IMS tools
- [HyperCCS](https://github.com/NeoNexusX/HyperCCS) - Hyperparameter CCS
- [Met4DX](https://met4dx.zhulab.cn/) - 4D metabolomics
- [MobiLipid](https://github.com/FelinaHildebrand/MobiLipid) - Lipid mobility
- [MOCCal](https://github.com/HinesLab/MOCCal) - Mobility calibration
- [Mol2CCS](https://github.com/enveda/ccs-prediction) - Molecule to CCS
- [OpenTIMS, TimsPy, and TimsR](www.github.com/michalsta/opentims) - TIMS data tools
- [PACCS](https://github.com/yuxuanliao/PACCS) - CCS prediction
- [PNNL PreProcessor](https://github.com/PNNL-Comp-Mass-Spec/PNNL-PreProcessor/releases/tag/v6.0) - PNNL preprocessing
- [POMICS](https://www.pomics.org/) - Proteomics and metabolomics
- [PubChemLite plus CCS](https://pubchemlite.lcsb.uni.lu) - PubChemLite with CCS
- [SigmaCCS](https://github.com/zmzhang/SigmaCCS) - Sigma CCS
- [Snakemake CCS](https://github.com/DasSusanta/snakemake_ccs) - Snakemake for CCS
- [SSN-CCSPIF](https://github.com/Anqi-Guo0105/Trendline-based-IM-MS-Feature-Filtering-Software-V1.0) - Trendline-based filtering

---

## 🧮 Isotopic

- [Aerith](https://github.com/thepanlab/Aerith) - Isotope analysis
- [CIL-ExPMRM](http://www.exposomemrm.com/) - Chemical isotope labeling
- [FAMetA](https://www.fameta.es) - Fatty acid metabolomics
- [IsoPairFinder](https://github.com/DoddLab/IsoPairFinder) - Isotope pair finding
- [isopair](https://github.com/kilgain/MassSpec/tree/main) - Isotope pairing
- [Isoreader](https://github.com/isoverse/isoreader) - Isotope ratio reading
- [IsoSolve](https://pypi.org/project/isosolve/) - Isotope solving
- [khipu](https://github.com/shuzhao-li-lab/khipu) - Isotope pattern analysis
- [Khipu-web](https://metabolomics.cloud/khipu) - Web-based khipu
- [SIMPEL](https://github.com/SIMPELmetabolism/SIMPEL) - Stable isotope metabolomics

---

## 🔭 Large Scale

- [LargeMetabo](https://github.com/LargeMetabo/LargeMetabo) - Large-scale metabolomics
- [peakPantheR](https://bioconductor.org/packages/peakPantheR/) - Peak integration at scale

---

## 🧬 Lipidomics

- [ADViSELipidomics](https://github.com/ShinyFabio/ADViSELipidomics) - Advanced lipidomics
- [BATL](complimet.ca/batl/) - Batch tools for lipidomics
- [DBLiPro](http://lipid.cloudna.cn/home) - Lipid database profiler
- [DIATAGe](https://github.com/Velenosi-Lab/DIATAGeR) - DIA tagging
- [Doxlipid](https://figshare.com/articles/journal_contribution/_em_strong_In_silico_strong_em_strong_-predicted_dynamic_oxlipidomics_MS_MS_library_High-throughput_discovery_and_characterization_of_unknown_oxidized_lipids_strong_/23700177/4) - Oxidized lipids library
- [LINEX2](https://gitlab.lrz.de/lipitum-projects/linex) - Lipid network explorer
- [LipiDetective](https://github.com/LipiTUM/lipidetective) - Lipid detection
- [LipiDex 2](https://github.com/coongroup/LipiDex-2) - Lipid identification
- [LipiDisease](http://cbdm-01.zdv.uni-mainz.de:3838/piyusmor/LipiDisease/) - Lipid-disease associations
- [Lipid Spectrum Generator](https://github.com/98104781/LSG) - Lipid spectrum generation
- [Lipid Wizard](https://github.com/RaoboXu/Lipidwizard) - Lipid analysis wizard
- [Lipid4Danalyzer](http://lipid4danalyzer.zhulab.cn/) - 4D lipid analysis
- [LipidA-IDER](https://github.com/Systems-Biology-Of-Lipid-Metabolism-Lab/LipidA-IDER) - Lipid A identification
- [LipidCruncher](https://lipidcruncher.org/) - Lipid data crunching
- [Lipidepifind](https://github.com/Xinguang-Liu/Lipidepifind) - Lipid epitope finder
- [LipidFinder 2.0](https://github.com/ODonnell-Lipidomics/LipidFinder) - Lipid finding
- [LipidIN](https://github.com/LinShuhaiLAB/LipidIN) - Lipid identification network
- [LipidOne 2.0](https://lipidone.eu) - Lipid analysis platform
- [LipidOracle](https://hub.docker.com/r/zambonilab/lipidoracle) - Lipid prediction
- [LipidOz](https://github.com/PNNL-m-q/lipidoz) - Ozone-induced lipid analysis
- [LipidQuant 1.0](https://holcapek.upce.cz/lipidquant.php) - Lipid quantification
- [LipidQuant 2.1](https://zenodo.org/record/10364585) - Updated lipid quantification
- [LipidSig](chenglab.cmu.edu.tw/lipidsig) - Lipid signatures
- [LipidSig 2.0](https://lipidsig.bioinfomics.org/) - Updated lipid signatures
- [LipidSigR](https://github.com/BioinfOMICS/LipidSigR) - R package for LipidSig
- [LipidSpace](https://lifs-tools.org/tools/lipidspace.html) - Lipid space visualization
- [LipidSuite](https://suite.lipidr.org/) - Lipid analysis suite
- [LIPID MAPS](https://www.lipidmaps.org/) - Lipid maps database
- [LipoCLEAN](https://github.com/stavis1/LipoCLEAN) - Lipidomics cleaning
- [LORA](https://lora.metabolomics.fgu.cas.cz/) - Lipid ontology reasoning
- [LPPtiger2](https://github.com/LMAI-TUD/lpptiger2) - Lipid pattern prediction
- [MS2Lipid](https://github.com/systemsomicslab/ms2lipid/tree/main/notebooks) - MS/MS lipid identification
- [MS-RIDD](https://github.com/Metabolomics-CompMS/MS-RIDD) - Lipid double bond analysis
- [Neurolipid Atlas](https://neurolipidatlas.com/) - Neural lipid atlas
- [RefLAS](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR002159) - Reference lipid analysis
- [RPLC-IMS-CID-MS Lipid Database](https://tarheels.live/bakerlab/databases/) - Lipid database
- [RTStaR](www.neurolipidomics.ca/rtstar/rtstar.html) - RT standardization
- [SimLipid](https://www.premierbiosoft.com/lipid/index.html) - Lipid simulation

---

## 🕸️ Metabolic Networks

- [DNEA](https://github.com/Karnovsky-Lab/DNEA/) - Differential network enrichment
- [FNICM](https://github.com/LiQi94/FNICM) - Flux network integration
- [KGMN](https://github.com/ZhuMetLab/MetDNA2) - Knowledge graph metabolic network
- [Kiphynet](http://pallab.cds.iisc.ac.in/kiphynet/) - Kinetic-phylogenetic network
- [M2R](https://github.com/e-weglarz-tomczak/m2r) - Metabolite to reaction
- [MetaboliticsDB](https://github.com/alperdokay/metabolitics-frontend) - Metabolite database
- [MetNet](https://github.com/simeoni-biolab/MetNet) - Metabolic network
- [MicroMap](https://dataverse.harvard.edu/dataverse/micromap) - Microbiome mapping
- [MInfer](https://github.com/cellbiomaths/MInfer) - Metabolite inference
- [MINNO](https://lewisresearchgroup.github.io/MINNO/) - Minimal network notation
- [NetAurHPD](https://github.com/TamirBar-Tov/NetAurHPD-Network-Auralization-Hyperlink-Prediction-Method) - Network auralization
- [Recon8D](https://github.com/sriram-lab/Recon8D) - Reconstruction 8D
- [Thermo-Flux](https://github.com/molecular-systems-biology/thermo-flux) - Thermodynamic flux

---

## 🗃️ Metadata

- [MatrixLM](https://github.com/senresearch/MatrixLM.jl) - Matrix linear models
- [MetaXtract](https://github.com/Rappsilber-Laboratory/MetaXtract) - Metadata extraction
- [PeakForest](https://github.com/peakforest/peakforest-webapp) - Peak forest database
- [SMetaS](https://github.com/metabolomics-us/metadatastandardizer) - Metadata standardization

---

## 🔀 Multifunctional

- [BUDDY](https://github.com/HuanLab/BUDDY) - Bottom-up deconvolution
- [DEIMoS](https://github.com/pnnl/deimos) - Data extraction for ion mobility
- [eMZed 3](https://emzed.ethz.ch) - Interactive analysis
- [FERMO](https://fermo.bioinformatics.nl/) - Feature-based molecular networking
- [iMAP](https://imap.metaboprofile.cloud/) - Integrated metabolomics analysis
- [maplet](https://github.com/krumsieklab/maplet) - Metabolomics analysis
- [Mass-Suite](https://github.com/XiminHu/mass-suite) - Mass spectrometry suite
- [matchms](https://github.com/matchms/matchms) - MS matching
- [MargheRita](https://github.com/emosca-cnr/margheRita) - R metabolomics
- [MetaboAnalyst 5.0](http://www.metaboanalyst.ca/) - Metabolomics analysis
- [MetaboAnalyst 6.0](https://www.metaboanalyst.ca) - Updated MetaboAnalyst
- [MetaboAnalystR 4.0](http://www.metaboanalyst.ca/) - R package for MetaboAnalyst
- [MetaboLink](https://github.com/anitamnd/MetaboLink) - Metabolomics linking
- [MetDNA3](http://metdna.zhulab.cn/) - Metabolite annotation
- [MetMiner, MDAtoolkits](https://github.com/ShawnWx2019/MetMiner) - Metabolite mining
- [MetaProViz](https://saezlab.github.io/MetaProViz/) - Metabolomics visualization
- [mpactR](https://www.mums2.org/mpactr/) - Metabolomics impact
- [MS-DIAL 5](https://github.com/systemsomicslab/MsdialWorkbench) - MS data analysis
- [MSOne](https://msone.claritybiosystems.com) - MS1-based analysis
- [MZmine 3](https://github.com/mzmine/mzmine) - Mass spectrometry analysis
- [NP3 MS Workflow](https://github.com/danielatrivella/NP3_MS_Workflow) - Natural products workflow
- [OpenMS 3](https://github.com/OpenMS/OpenMS) - Open-source MS analysis
- [OpenMS WebApps](https://github.com/OpenMS/streamlit-template) - Web applications
- [patRoon 2](https://github.com/rickhelmus/patRoon) - Pattern recognition
- [PMart](https://github.com/pmartR/PMart_ShinyApp) - Proteomics/metabolomics
- [POMAShiny](https://github.com/nutrimetabolomics/POMAShiny) - POMA Shiny app
- [Punc'data](https://github.com/WTVoe/puncdata) - Punctured data analysis
- [Rodin](https://rodin-meta.com/) - Metabolomics platform
- [SLAW](https://github.com/adelabriere/SLAW) - Sample-level analysis
- [TidyMass](https://www.tidymass.org/) - Tidy mass spectrometry
- [TidyMass2](https://github.com/tidymass/tidymass) - Updated TidyMass
- [TraceMetrix](https://www.biosino.org/tracemetrix/) - Trace analysis
- [UmetaFlow](https://github.com/biosustain/snakemake_UmetaFlow) - Untargeted metabolomics
- [XCMS-METLIN](https://xcmsonline.scripps.edu/) - XCMS online

---

## 🔗 Multiomics

- [AgeAnnoMO](https://relab.xidian.edu.cn/AgeAnnoMO/#/) - Age annotation multiomics
- [BATMAN‑TCM 2.0](http://bionet.ncpsb.org.cn/batman-tcm/) - TCM database
- [BnIR](http://yanglab.hzau.edu.cn/BnIR) - Brassica napus integration
- [haCCA](https://github.com/LittleLittleCloud/haCCA) - Hierarchical CCA
- [HoloFoodR](https://doi.org/10.18129/B9.bioc.HoloFoodR) - Hologenome food
- [INTEGRATE](https://github.com/qLSLab/integrate) - Multiomics integration
- [iSODA](http://isoda.online/) - Integrated systems omics
- [iTraNet](https://itranet.streamlit.app) - Integrative network
- [jMorp](https://jmorp.megabank.tohoku.ac.jp) - Japanese multiomics
- [MetaNet](https://github.com/Asa12138/MetaNet) - Meta network analysis
- [MiMeNet](https://github.com/YDaiLab/MiMeNet) - Microbiome-metabolome network
- [MODMS](https://modms.lzu.edu.cn) - Multiomics data management
- [MPOD](http://mpod.ict.cn) - Multiomics platform
- [MultiOmicsIntegrator](https://github.com/ASAGlab/MOI--An-integrated-solution-for-omics-analyses) - Omics integration
- [OmicsAnalyst](https://www.omicsanalyst.ca) - Omics analysis
- [Omics Dashboard](https://biocyc.org/dashboard/dashboard-help.html) - Omics visualization
- [OmicsNet 2.0](https://www.omicsnet.ca) - Omics networking
- [Paired Omics](https://github.com/iomega/paired-data-form) - Paired omics data
- [SOmicsFusion](10.5281/zenodo.7700528) - Omics fusion
- [SVAtlas](https://www.svatlas.org/) - Structural variant atlas
- [TurbOmics](https://proteomics.cnic.es/TurboPutative/TurbOmicsApp.html) - Turbocharged omics
- [TurboPutative](https://proteomics.cnic.es/TurboPutative/) - Turbo putative annotation

---

## 🔬 NMR

- [A-SIMA/A-MAP](https://poky.clas.ucdenver.edu) - NMR analysis
- [AQuA](https://pmc.ncbi.nlm.nih.gov/articles/instance/8253485/bin/ac0c04233_si_001.pdf) - Automated quantification
- [ASICS](https://github.com/GaelleLefort/ASICS) - Automatic spectrum identification
- [Bucket Fuser](https://www.mdpi.com/article/10.3390/metabo12090812/s1) - Bucket fusion
- [CASMDB](https://ccpn.ac.uk/software/analysismetabolomics/) - CAS metabolomics database
- [COLMAR1d](https://spin.ccic.osu.edu/index.php/colmar1d) - 1D COLMAR
- [COLMARppm](https://spin.ccic.osu.edu/index.php/colmar) - COLMAR ppm
- [COLMARvista](https://github.com/lidawei1975/colmarvista) - COLMAR visualization
- [CReSS](https://github.com/Qihoo360/CReSS) - Cross-referenced spectroscopy
- [DeepSAT](https://github.com/mwang87/DeepSAT) - Deep learning for SAT
- [FlavorFormer](https://github.com/yfWang01/FlavorFormer) - Flavor analysis
- [InRA](https://github.com/InRA-Software/InRA) - Integrated NMR analysis
- [LAMAIS](https://github.com/Ariel-foliage/LAMAIS) - Language model for MS
- [MADByTE](https://github.com/liningtonlab/madbyte) - Bioactivity typing
- [Magnetstein](https://github.com/ElsevierSoftwareX/SOFTX-D-25-00358) - NMR processing
- [MagMet](http://magmet.ca/) - Magnetic metabolomics
- [MagMet-F](https://www.magmet.ca/) - MagMet fluorescence
- [mcfNMR](https://github.com/GeoMetabolomics-ICBM/mcfNMR) - MCF NMR
- [MetaboLabPy](https://github.com/ludwigc/metabolabpy) - Python metabolomics
- [MetAssimulo 2.0](https://github.com/yanyan5420/MetAssimulo_2) - Metabolite simulation
- [MixONat](https://sourceforge.net/projects/mixonat/) - Mixture analysis
- [MultiNMRFit](https://github.com/NMRTeamTBI/MultiNMRFit/) - Multi-NMR fitting
- [NAPROC-13](https://github.com/DIFACQUIM/naproc13_characterization) - Natural products 13C
- [NMR2Struct](https://github.com/MarklandGroup/NMR2Struct) - NMR to structure
- [NMRformer](https://github.com/zza1211/NMRformer) - Transformer for NMR
- [NMRFx](https://github.com/nanalysis/nmrfx) - NMR processing
- [NMRhub](https://nmrhub.org/) - NMR hub
- [NMRInversions.jl](https://github.com/aris-mav/NMRInversions.jl) - NMR inversions in Julia
- [NMRium](https://www.nmrium.org/teaching) - NMR in browser
- [NMR molecular networking](https://github.com/enveda/NMR-Networking) - NMR networking
- [NMRphasing](https://github.com/ajiangsfu/NMRphasing) - NMR phasing
- [NMRQNet](https://github.com/LiuzLab/NMRQNet) - NMR quantification network
- [nmRanalysis](https://github.com/EMSL-Computing/nmRanalysis) - R NMR analysis
- [NP-MRD](https://np-mrd.org/) - Natural products MRD
- [PRIMA-Panel](https://github.com/funkam/PRIMA) - PRIMA panel
- [PROSPRE](https://prospre.ca/) - Prospective NMR
- [Protomix](https://github.com/mzniber/Protomix) - Protomix analysis
- [pSCNN](https://github.com/yuxuanliao/pSCNN) - Parallel spectral CNN
- [PyINETA](https://github.com/edisonomics/PyINETA) - Python INETA
- [ROIAL-NMR](https://github.com/Leo-Cheng-Lab/ROIAL-NMR) - ROIAL for NMR
- [SAND](https://github.com/edisonomics/SAND) - Spectral analysis
- [SMolESY](https://github.com/pantakis/SMolESY-select) - Small molecule ESY
- [ukbnmr](https://github.com/sritchie73/ukbnmr/) - UK Biobank NMR

---

## 🌱 Organism Specific

- [Aspergillus Metabolome Database](https://github.com/albertogilf/ceuMassMediator/tree/master/CMMAspergillusDB) - Aspergillus metabolome
- [FoodAtlas](http://foodatlas.ai) - Food metabolomics
- [GMDP Database](https://lccl.shinyapps.io/GMDP/) - Gut microbiome database
- [gutMGene v2.0](http://bio-computing.hrbmu.edu.cn/gutmgene) - Gut microbiome genes
- [gutSMASH](https://gutsmash.bioinformatics.nl/) - Gut secondary metabolites
- [GutUDB](https://intestine.splicedb.net/) - Gut database
- [HORDB](https://github.com/CPU-HORDB/HORDB) - Herbal odor database
- [IDBac](https://github.com/chasemc/IDBacApp) - Bacterial identification
- [LiqBioer](http://www.medsysbio.org:8080/LiqBioer/) - Liquid biopsy
- [MDSi](http://sky.sxau.edu.cn/MDSi.htm) - Maize database
- [MetaDb](http://medmetadb.ynau.edu.cn) - Medicinal metabolome database
- [MetaboSeek](https://github.com/mjhelf/Metaboseek) - Metabolomics seeking
- [MiMIR](https://github.com/DanieleBizzarri/MiMIR) - Microbiome integration
- [mVOC 4.0](http://bioinformatics.charite.de/mvoc) - Microbial VOCs
- [MyxoDB](https://figshare.com/collections/Constructing_a_Myxobacterial_Natural_Product_Database_to_Facilitate_NMR-Based_Metabolomics_Bioprospecting_of_Myxobacteria/6468946) - Myxobacteria database
- [Omics Untargeted Key Script](https://github.com/plyush1993/OUKS) - OUKS
- [OrchidMD](www.orchidcomics.com) - Orchid metabolomics
- [PaIRKAT](https://github.com/CharlieCarpenter/PaIRKAT) - Pathway integration
- [PeanutOmics](https://cgm.sjtu.edu.cn/PeanutOmics/) - Peanut omics
- [PSC-db](http://pscdb.appsbio.utalca.cl/) - Plant secondary compounds
- [SalivaDB](https://webs.iiitd.edu.in/raghava/salivadb/) - Saliva database
- [SmilODB](http://www.isage.top:56789) - Smile odor database
- [StreptomeDB 4.0](https://streptomedb.vm.uni-freiburg.de/streptomedb/) - Streptomyces database
- [The Molecular Human](https://littleswissriver.shinyapps.io/generate_network/) - Human molecular network

---

## 🛤️ Pathway, Enrichment and Ontology Tools

- [ChemFOnt](https://www.chemfont.ca/) - Chemical ontology
- [EnrichMET](https://github.com/biodatalab/enrichmet) - Enrichment for metabolomics
- [GINv2.0](https://github.com/BIGchix/GINv2.0) - Gene interaction network
- [IDSL.GOA](https://goa.idsl.me/) - Gene ontology annotation
- [iMSEA](https://github.com/BioNet-XMU/iMSEA) - Integrative MSEA
- [Lilikoi V2.0](https://github.com/lanagarmire/lilikoi2) - Pathway analysis
- [massDatabase](https://massdatabase.tidymass.org/) - Mass database
- [MBROLE3](http://csbg.cnb.csic.es/mbrole3) - Metabolite biological role
- [metGWAS 1.0](https://github.com/saifurbd28/metGWAS-1.0) - Metabolite GWAS
- [metLinkR](https://github.com/metLinkR/metLinkR) - Metabolite linking
- [MetChem](https://CRAN.R-project.org/package=MetChem) - Metabolite chemistry
- [MIMOSA2](http://www.borensteinlab.com/software_MIMOSA2.html) - Microbiome metabolome
- [MiNEApy](https://github.com/vpandey-om/mineapy) - Minimal network enrichment
- [MS2MP](https://github.com/ucasaccn/MS2MP) - MS to metabolic pathway
- [ORA](https://github.com/cwieder/metabolomics-ORA) - Over-representation analysis
- [PALS](https://pals.glasgowcompbio.org/app/) - Pathway activity levels
- [PaintOmics 4](https://paintomics.org) - Pathway visualization
- [PathBank 2.0](https://pathbank.org/) - Pathway database
- [SBGNview](https://github.com/datapplab/SBGNview) - SBGN visualization
- [secCellFie](https://github.com/LewisLabUCSD/secCellFie) - Cell-type metabolic functions
- [Viime-Path](https://github.com/girder/viime-path) - Pathway analysis
- [WebGestalt](https://www.webgestalt.org) - Web-based enrichment
- [Xconnector](https://github.com/Proteomicslab57357/Xconnector) - Cross-omics connector

---

## 🔎 Patterns

- [MCnebula](https://mcnebula.org/) - Multiple classes nebula
- [Metabokiller](https://pypi.org/project/Metabokiller/) - Pattern killing
- [MS2LDA 2.0](https://github.com/vdhooftcompmet/MS2LDA) - LDA for MS/MS
- [mzBucket](https://github.com/hildebrandtlab/mzBucket) - m/z bucketing

---

## ⚙️ Pre-processing

- [3D-MSNet](https://github.com/CSi-Studio/3D-MSNet) - 3D MS network
- [AriumMS](https://github.com/AdrianHaun/AriumMS/) - Arium MS processing
- [asari](https://github.com/shuzhao-li-lab/asari) - Feature extraction
- [AVIR](https://github.com/HuanLab/AVIR) - Adduct and in-source verification
- [CorrDIA](https://github.com/zhangxiang1205/CorrDIA) - Correlation-based DIA
- [CycloBranch](https://ms.biomed.cas.cz/cyclobranch) - Cyclopeptide analysis
- [DaDIA](https://github.com/HuanLab/DaDIA) - Data-dependent acquisition
- [DBDIpy](https://github.com/leopold-weidner/DBDIpy) - DBDI processing
- [DecoID](https://github.com/pattilab/DecoID/releases) - Deconvolution ID
- [DeepMSProfiler](https://github.com/yjdeng9/DeepMSProfiler) - Deep learning profiler
- [DIAMetAlyzer](https://github.com/OpenMS/DIAMetAlyzer) - DIA analysis
- [DNMS2Purifier](https://github.com/HuanLab/DNMS2Purifier) - MS/MS purification
- [DuReS](https://github.com/BiosystemEngineeringLab-IITB/dures) - Duplicate removal
- [Eclipse](https://github.com/broadinstitute/bmxp/tree/main/bmxp/eclipse) - Eclipse processing
- [EVA](https://github.com/HuanLab/EVA) - Enhanced visualization
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
- [maplet](https://github.com/krumsieklab/maplet) - Metabolomics mapping
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
- [mzLearn](https://github.com/ReviveMed/mzEmbed) - m/z learning
- [mzRAPP](https://github.com/YasinEl/mzRAPP) - Rapid processing
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

- [ALISTER](https://itmp.shinyapps.io/alister/) - Quality assessment
- [CordBat](https://github.com/BorchLab/CordBat) - Batch correction
- [Dbnorm](https://github.com/NBDZ/dbnorm) - Database normalization
- [hRUV](https://github.com/SydneyBioX/hRUV) - Hierarchical RUV
- [iDIA-QC](https://github.com/guomics-lab/iDIA-QC) - DIA quality control
- [InjectionDesign](https://github.com/CSi-Studio/InjectionDesign) - Injection optimization
- [MAFFIN](https://github.com/HuanLab/MAFFIN) - Mass feature filtering
- [marr](http://bioconductor.org/packages/release/bioc/html/marr.html) - Maximum rank reproducibility
- [MassQLab](https://github.com/JohnsonDylan/MassQLab) - Mass quality lab
- [MatrixQCvis](https://bioconductor.org/packages/MatrixQCvis/) - Matrix QC visualization
- [MetaMOPE](https://metamope.cmdm.tw) - Metabolomics MOPE
- [MetaPro](https://github.com/CSi-Studio/MetaPro) - Meta processing
- [MSNORM](https://pubs.acs.org/doi/suppl/10.1021/jasms.3c00295/suppl_file/js3c00295_si_001.txt) - MS normalization
- [Msquality](https://bioconductor.org/packages/MsQuality) - MS quality metrics
- [MultiBaC](https://github.com/ConesaLab/MultiBaC) - Multiple batch correction
- [mzQuality](https://github.com/hankemeierlab/mzQuality) - m/z quality
- [OSCA Finder](https://github.com/dawnmengsjtu/OSCA-Finder) - Outlier finding
- [Paramounter](https://github.com/HuanLab/Paramounter) - Parameter optimization
- [PeakQC](https://github.com/pnnl/IonToolPack) - Peak quality control
- [QC4Metabolomics](https://github.com/stanstrup/QC4Metabolomics) - QC for metabolomics
- [QComics](https://github.com/ricoderks/QComics) - QC for omics
- [QuantyFey](https://github.com/CDLMarkus/QuantyFey) - Quantification quality
- [RALPS](https://github.com/zamboni-lab/RALPS) - Regularized adversarial learning
- [Rapid QC-MS](https://github.com/czbiohub-sf/Rapid-QC-MS) - Rapid quality control
- [RawBeans](https://bitbucket.org/incpm/prot-qc/downloads/) - Raw data QC
- [RawHummus](https://github.com/YonghuiDong/RawHummus) - Raw data visualization
- [Shinyscreen](https://gitlab.com/uniluxembourg/lcsb/eci/shinyscreen) - Shiny screening
- [WaveICA 2.0](https://github.com/dengkuistat/WaveICA_2.0) - Wavelet ICA

---

## ⏱️ RT (Retention Time)

- [ABCoRT](https://github.com/RiverCCC/ABCoRT) - Antibody RT prediction
- [AsRTNet](https://github.com/piscookie/AS-RT) - Attention-based RT
- [CMM-RT](https://github.com/constantino-garcia/cmmrt) - CMM for RT
- [DeepGCN-RT](https://github.com/kangqiyue/DeepGCN-RT) - Deep GCN for RT
- [DeepRTAlign](https://github.com/PHOENIXcenter/deeprtalign) - Deep RT alignment
- [GNN-RT](https://github.com/Qiong-Yang/GNN-RT) - GNN for RT prediction
- [Graphormer-RT](https://github.com/HopkinsLaboratory/Graphormer-RT) - Graphormer for RT
- [QGeoGNN](https://github.com/woshixuhao/Retention-Time-Prediction-for-Chromatographic-Enantioseparation/tree/main) - Geometric GNN for RT
- [ReTimeML](https://mikeallwright23-retime-app-lipid3-021zpv.streamlit.app/) - RT machine learning
- [retention_time_gnn](https://github.com/seokhokang/retention_time_gnn) - RT GNN
- [RI-based-CPSC](https://github.com/WHU-Fenglab/RI-based-CPSC) - RI-based compound prediction
- [ROASMI](https://github.com/FangYuan717/ROASMI) - RT optimization
- [RT-Ensemble Pred](https://gitlab.com/mikic93/rt-ensemble-pred) - Ensemble RT prediction
- [RT-Pred](https://rtpred.ca/) - RT prediction
- [RT-Transformer](https://github.com/01dadada/RT-Transformer) - Transformer for RT

---

## 🦠 Single Cell Metabolomics

- [MetaPhenotype](https://github.com/songyuan93/MetaPhenotype) - Metabolic phenotyping
- [MMEASE](https://idrblab.org/mmease/) - Single cell metabolomics
- [scFUMES](https://github.com/ChengF-Lab/scFUMES) - Single cell functional metabolomics
- [SCMeTA](https://github.com/ChenLab/SCMeTA) - Single cell meta-analysis

---

## 🗺️ Spatial Metabolomics

- [scSpaMet](https://github.com/coskunlab/ScSpaMet) - Spatial metabolomics
- [SMART](https://github.com/bioinfo-ibms-pumc/SMART) - Spatial metabolomics analysis
- [SMAnalyst](https://github.com/mzlab-research/SManalyst) - Spatial metabolomics analyst
- [SMQVP](https://github.com/mzlab-research/SMQVP) - Spatial quality visualization
- [SpaMTP](https://github.com/GenomicsMachineLearning/SpaMTP) - Spatial metabolomics platform

---

## 🧩 Specialized

- [Amanida](https://github.com/mariallr/amanida) - Amanita metabolomics
- [arcMS](https://github.com/leesulab/arcMS) - Archaea MS
- [ATLASx](https://lcsb-databases.epfl.ch/Atlas2) - Comprehensive atlas
- [BAGO](https://github.com/huaxuyu/bago) - Bioactivity-guided optimization
- [BioPKS Pipeline, RetroTide](https://github.com/JBEI/BioPKS-Pipeline) - Polyketide analysis
- [BioTransformer 3.0](https://biotransformer.ca) - Biotransformation prediction
- [BitterMasS](https://github.com/Niv-Lab/BitterPredict1) - Bitter compound prediction
- [BreathXplorer](https://github.com/HuanLab/breathXplorer) - Breath analysis
- [CCWeights](https://github.com/YonghuiDong/CCWeights) - CCS weighting
- [ChemEcho](https://github.com/biorack/chemecho) - Chemical echo
- [ChemSpectra](https://github.com/ComPlat/chem-spectra-app) - Chemical spectra app
- [COMMIT](https://github.com/pwendering/COMMIT/tree/v1.0.0) - Community integration
- [COVRECON](https://bitbucket.org/mosys-univie/covrecon) - COVID reconstruction
- [DarkNPS](https://github.com/skinniderlab/DarkNPS) - Dark chemical space
- [DeepMet](https://deepmet.org/) - Deep learning metabolomics
- [Easy-Amanida](https://github.com/mariallr/easy-amanida) - Easy Amanida
- [ECDFormer](https://github.com/HowardLi1984/ECDFormer) - ECD transformer
- [FORUM](https://github.com/eMetaboHUB/Forum-DiseasesChem) - Disease-chemical forum
- [GlyKAn AZ](https://github.com/ZaylaSchaeffer/GlyKAn-AZ-application) - Glycan analysis
- [homologueDiscoverer](https://github.com/kevinmildau/homologueDiscoverer) - Homologue discovery
- [IonFlow](https://github.com/wanchanglin/ionflow) - Ion flow analysis
- [M2S](https://github.com/rjdossan/M2S) - Metabolite to structure
- [Maldi Transformer](https://github.com/gdewael/maldi-nn) - MALDI analysis
- [MassQL](https://github.com/mwang87/MassQueryLanguage) - Mass query language
- [MEISTER](https://github.com/richardxie1119/MEISTER) - Metabolite identification
- [MetaboliteChat](https://github.com/13501274828/MetaboliteChat) - LLM for metabolites
- [MetaboLM](https://github.com/Qiu-Shizheng/MetaboLM) - Language model for metabolomics
- [Metaboverse](https://github.com/Metaboverse/Metaboverse) - Metabolic visualization
- [MetCohort](https://github.com/JunYang2021/MetCohort) - Cohort metabolomics
- [MetDIT](https://github.com/Li-OmicsLab/MetDIT) - Metabolite data integration
- [mGWAS-Explorer 2.0](https://www.mgwas.ca/) - Metabolite GWAS explorer
- [MINE 2.0](https://minedatabase.ci.northwestern.edu) - Metabolomics MINE
- [MODAPro](https://github.com/zhaoxiaoqi0714/MODAPro) - MODA processing
- [MolSpectLLM](https://github.com/Eurekashen/MolSpectLLM) - LLM for spectroscopy
- [MSConnect](https://github.com/RTKlab-BYU/MSConnect) - MS connectivity
- [MSident](https://github.com/cplqam/GUI) - MS identification
- [MS-RT](https://github.com/XianghuWang-287/Metabolomics_Clustering_Benchmark) - MS-RT clustering
- [MS2toImg](https://github.com/jsehhs/MS2toImg-CNN) - MS/MS to image
- [Multipass CCS Refiner](https://ericgier.shinyapps.io/Multipass-CCS-Refiner/) - CCS refinement
- [NP Analyst](https://www.npanalyst.org/) - Natural products analysis
- [NPFimg](http://github.com/poomcj/NPFimg) - Natural products fingerprint
- [Pickaxe](https://github.com/tyo-nu/MINE-Database) - Reaction rule mining
- [ptairMS](https://bioconductor.org/packages/release/bioc/html/ptairMS.html) - PTR-MS analysis
- [SMITER](https://github.com/LeidelLab/SMITER) - Structure elucidation
- [SPIFFED](https://github.com/bio-it-station/SPIFFED) - Spectral processing
- [Umatrix](https://cran.r-project.org/package=Umatrix) - Unified matrix
- [VIMMS 2.0](https://github.com/glasgowcompbio/vimms/) - Virtual metabolomics MS

---

## 📕 Spectral Library

- [CIeaD](https://www.moleculardetective.org/Links.html) - Chemical investigation
- [DNA Adduct Portal](https://sites.google.com/umn.edu/dnaadductportal) - DNA adducts
- [EMBL-MCF 2.0](https://www.embl.org/groups/metabolomics/instrumentation-and-software/#MCF-library) - EMBL library
- [GNPS Drug Library](https://github.com/ninahaoqizhao/Manuscript_GNPS_Drug_Library) - Drug spectra
- [HighResNPS](https://highresnps.com/) - High-resolution NPS
- [In silico infrared spectral library of human metabolites](https://zenodo.org/records/7706021) - IR library
- [MassBank](https://massbank.jp/) - Mass spectral database
- [metScribeR](https://github.com/ncats/metScribeR) - Metabolite scribing
- [mFAM](https://github.com/MassBank/MassBank-data/commit/a91b1ca4841aea536545f7c1d452c1f80d225e84) - Mass family
- [MIADB](https://gnps.ucsd.edu) - MIADB database
- [MS2extrcat](https://www.cooperstonelab.com/MS2extract/) - MS/MS extraction
- [PASL](https://gnps.ucsd.edu/ProteoSAFe/gnpslibrary.jsp?library=PYRROLIZIDINE-ALKALOID-SPECTRAL-LIBRARY) - Pyrrolizidine alkaloids
- [RMB-mix](https://github.com/rmassbank/rmbmix) - RMassBank mixing
- [UCL-MetIsoLib](https://github.com/kserafimov10/UCLMetIsoLib) - UCL isotope library
- [WFSR Food Safety Mass Spectral Library](https://www.wur.nl/en/show/food-safety-mass-spectral-library.htm) - Food safety library

---

## 📈 Statistical

- [cluster-CV](https://github.com/cluster-cv) - Cluster cross-validation
- [DisCo P-ad](https://github.com/KechrisLab/DisCoPad) - Discovery co-analysis
- [GetFeatistics](https://github.com/FrigerioGianfranco/GetFeatistics) - Feature statistics
- [imputomics](github.com/BioGenies/imputomics) - Imputation methods
- [MAI](https://bioconductor.org/packages/release/bioc/html/MAI.html) - Multiple annotation imputation
- [MAMSI](https://pypi.org/project/mamsi/) - Multi-method analysis
- [MetaboLINK/ PCA-GLASSO](https://github.com/jlichtarge/pcaGLASSO) - PCA GLASSO
- [MetaHD](https://bookdown.org/a2delivera/MetaHD/) - High-dimensional metabolomics
- [MeTEor](https://github.com/scibiome/meteor) - Metabolite testing
- [MultiClassMetabo](http://idrblab.cn/multiclassmetabo/) - Multi-class classification
- [NOREVA](https://idrblab.org/noreva/) - Normalization evaluation
- [omicsMIC](https://github.com/WQLin8/omicsMIC) - Omics integration
- [PLSKO](https://github.com/guannan-yang/PLSKO) - PLS knockoff
- [PredCMB](https://www.sysbiolab.org/predcmb) - Prediction for CMB
- [rMisbeta](https://cran.r-project.org/package=rMisbeta) - R Misbeta
- [SMART 2.0](https://github.com/YuJenL/SMART) - Statistical metabolomics
- [SpectraX](https://github.com/weantony/SpectraX) - Spectral analysis
- [wscaling.R](https://github.com/nishithkumarpaul/robustScaling/blob/main/wscaling.R) - Weighted scaling

---

## 🎯 Targeted

- [automRm](https://gitlab.gwdg.de/joerg.buescher/automrm) - Automated MRM
- [CLAW-MRM](https://github.com/chopralab/CLAW) - CLAW for MRM
- [MRMPro](http://mrmpro.csibio.com/) - MRM processing
- [MRMQuant](https://github.com/kslynn128171/MRMQuant) - MRM quantification
- [Norm ISWSVR](https://github.com/Dingxian666/Norm-ISWSVR) - IS normalization
- [RHermes](https://github.com/RogerGinBer/RHermes) - Targeted screening
- [SCALiR](https://lewisresearchgroup-ms-conc-streamlit-app-w9e64f.streamlit.app/) - Scaling in R
- [TARDIS](https://github.com/UGent-LIMET/TARDIS) - Targeted analysis

---

## 📊 Visualization

- [ClusterApp](http://ccbl-apps.fcfrp.usp.br/ClusterApp) - Clustering application
- [EDaViS](https://github.com/LRuehrmund/RoMBAT) - Enhanced data visualization
- [GraphBio](https://github.com/databio2022/GraphBio) - Graph visualization
- [lcmsWorld](https://github.com/PGB-LIV/lcmsWorld) - LC-MS world
- [MODE](https://github.com/pmartR/MODE_ShinyApp) - Metabolomics data explorer
- [MS-VIS](https://github.com/hkevgill/MS-VIS) - MS visualization
- [PathwayNexus](www.cls.uni-konstanz.de/software/pathway-nexus) - Pathway visualization
- [RAIChU](https://github.com/BTheDragonMaster/RAIChU) - Chemical visualization
- [RDD metabolomics platform](https://gnps-rdd.streamlit.app/) - RDD visualization
- [SpecXplore](https://github.com/kevinmildau/specXplore) - Spectral exploration

---

## 📄 License

This project is licensed under the terms in the [LICENSE](LICENSE) file.
