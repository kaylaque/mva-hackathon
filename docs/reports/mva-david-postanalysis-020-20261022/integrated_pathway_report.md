# Integrated pathway report

Research-only post-analysis of saved g:Profiler and DAVID exports. Descriptive annotation is not
efficacy evidence, not a diagnostic finding, and not independent replication.

## Definitions used

- **Pathway-supported**: a provider reported the term below its own adjusted-q threshold; each provider is judged on its own statistics.
- **DAVID sensitivity threshold**: Benjamini <= 0.2; g:Profiler BH_Q remains at 0.05.
- **Nonspecific for an MVA mechanism**: the term is not itself a spindle-assembly-checkpoint, chromosome-segregation or cohesin process, so it is context rather than mechanism.
- **Modification annotation**: membership in a post-translational-modification term names members, not measured activity.
- **not_comparable**: the other provider produced no row for this term, so no cross-provider statement is made.
- **Jaccard**: gene overlap between two complete nonempty contributor sets; unavailable sets report nothing rather than an overlap of zero.

## Verified provenance identities

- run_key: `2acba21c73a0c6236cd2cf701f3097601632e3faff27906e21250379ddc3bf11`
- import_version: `5a3cf729a7471ef7e91ece6837821fe346b31aa8e448ce21ed3705d5e5e6cbc3`
- baseline projection SHA-256: `{'p0_gene_qc.csv': '4b977f06add4b2c4be16fd63f5c7b9cabb82fb95aafa5bd4ff9abb503b212bcb', 'p1_enrichment_results.csv': 'faee884061e8d76b11c85d7060d88f4725f575d6a87755db770cba1bbae6a8e6'}`
- analysis threshold: 0.2
- DAVID categories: OMIM_DISEASE, UP_KW_DISEASE, COG_ONTOLOGY, UP_KW_BIOLOGICAL_PROCESS, UP_KW_CELLULAR_COMPONENT, UP_KW_MOLECULAR_FUNCTION, UP_KW_PTM, UP_SEQ_FEATURE, GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, UP_KW_LIGAND, KEGG_PATHWAY, INTERPRO, PIR_SUPERFAMILY, SMART, UP_KW_DOMAIN
- distinct nonempty clusters: 226
- cluster annotation rows: 1812

## Native provider results

### DAVID-significant rows

- GO:0005515: protein binding (native Benjamini 0.2; g:Profiler not_comparable)
- GO:0005829: cytosol (native Benjamini 0.0776; g:Profiler not_comparable)
- nan: Acetylation (native Benjamini 3.94e-08; g:Profiler not_comparable)
- nan: Phosphoprotein (native Benjamini 0.00612; g:Profiler not_comparable)
- nan: Host-virus interaction (native Benjamini 0.114; g:Profiler not_comparable)
- nan: Cytoskeleton (native Benjamini 0.0599; g:Profiler not_comparable)
- nan: Mitochondrion (native Benjamini 0.0623; g:Profiler not_comparable)
- nan: Nucleus (native Benjamini 0.0805; g:Profiler not_comparable)
- nan: Cytoplasm (native Benjamini 0.0805; g:Profiler not_comparable)
- nan: Endosome (native Benjamini 0.0805; g:Profiler not_comparable)

### g:Profiler-significant unique terms

- none: no g:Profiler row met its own adjusted-q threshold

### Unmapped significant annotations requiring review

- WIKIPATHWAYS EV-356d21414978e80464655eaa: not_significant / not_comparable
- UNMAPPED EV-f71fc83c355f3793c61d9cc1: not_comparable / significant
- UNMAPPED EV-0d85fa2a802eea461ef3f3ce: not_comparable / significant
- UNMAPPED EV-9a26e07c29af22dfbe9e04f2: not_comparable / significant
- UNMAPPED EV-23f8865dfab0979acd9ca0b5: not_comparable / significant
- UNMAPPED EV-7f327d0210c848d1c741423b: not_comparable / significant
- UNMAPPED EV-3a5566e96b69a97197b02752: not_comparable / not_significant
- UNMAPPED EV-bb854e625b3ae237ef194534: not_comparable / significant
- UNMAPPED EV-5c00676f0a3a21cb45623dea: not_comparable / significant
- UNMAPPED EV-9e78ca1bc58cce572f0a5524: not_comparable / not_significant
- UNMAPPED EV-0a0b83123a57aa8af58d41e6: not_comparable / not_significant
- UNMAPPED EV-48d16d7123fb1f329858bb43: not_comparable / significant
- UNMAPPED EV-f0056eccde1f7496f3eeb756: not_comparable / not_significant
- UNMAPPED EV-6d4be3f80fb1cea496869531: not_comparable / not_significant
- UNMAPPED EV-ff1c2c621595e00ff3b2b0e2: not_comparable / not_significant
- UNMAPPED EV-d6f494b07a5cf769305993a0: not_comparable / not_significant
- UNMAPPED EV-e2d7e52de6a7a5ee573d1bb9: not_comparable / not_significant
- UNMAPPED EV-0f856a70348a92fcb131cb26: not_comparable / not_significant
- UNMAPPED EV-4316a83a1116b1fc1b5efdfa: not_comparable / not_significant
- UNMAPPED EV-c91bf5d3f3177b6f56b862c8: not_comparable / not_significant
- UNMAPPED EV-590a35df25f7e2800897a406: not_comparable / not_significant
- UNMAPPED EV-0c7c00ec9ccca7c45ad6df23: not_comparable / not_significant
- UNMAPPED EV-8f1d815403a0229c30df112e: not_comparable / not_significant
- UNMAPPED EV-264f5b302b1db7bf701e5ca4: not_comparable / not_significant
- UNMAPPED EV-cc0975224a2fbe1840cb5cb9: not_comparable / not_significant
- UNMAPPED EV-539d10d7bea953b9353d5c37: not_comparable / not_significant
- UNMAPPED EV-75a64a0f5b59d163cc703896: not_comparable / not_significant
- UNMAPPED EV-9dade5dad0eeee89d91700b6: not_comparable / not_significant
- UNMAPPED EV-0539c2df0a182e288aaab6e0: not_comparable / not_significant
- UNMAPPED EV-80a1472aae4066bba21cff91: not_comparable / not_significant
- UNMAPPED EV-872232d1bb9a8a91b007dd65: not_comparable / not_significant
- UNMAPPED EV-912cea68df8783adb9990e4e: not_comparable / not_significant
- UNMAPPED EV-b62d6748e2b20aa295cf1403: not_comparable / not_significant
- UNMAPPED EV-4f45bb0287b1c9e310ea3214: not_comparable / not_significant
- UNMAPPED EV-31065c04b0738ada16a0a6c8: not_comparable / not_significant
- UNMAPPED EV-bc1431aa988d847f52ccecff: not_comparable / not_significant
- UNMAPPED EV-9aa6f0a0b89afb935d525010: not_comparable / not_significant
- UNMAPPED EV-3c32a94f94872203bc8c99dc: not_comparable / not_significant
- UNMAPPED EV-bb6b121965ceb6a4af0e9456: not_comparable / not_significant
- UNMAPPED EV-1458ab212640af0770e7591e: not_comparable / not_significant
- UNMAPPED EV-0e33b0142f6ad06a422debb2: not_comparable / not_significant
- UNMAPPED EV-0adf35ab23ef63866aac953d: not_comparable / not_significant
- UNMAPPED EV-e432b3e5e79ad6c32d330ff0: not_comparable / not_significant
- UNMAPPED EV-b71f045bbaed91cde46e51be: not_comparable / not_significant
- UNMAPPED EV-a39bd0f5daeb43f9e4f0d153: not_comparable / not_significant
- UNMAPPED EV-e1432be58498417ea1cfc6a0: not_comparable / not_significant
- UNMAPPED EV-8b5b6cf1d3bccc74f7c8e557: not_comparable / not_significant
- UNMAPPED EV-0d3af674617307d1d2042717: not_comparable / not_significant
- UNMAPPED EV-fbc22d5c53390a0a81beb442: not_comparable / not_significant
- UNMAPPED EV-7f58e3e87958f1acb2c6dc3b: not_comparable / not_significant
- UNMAPPED EV-52d5c48eb1a99017da10b3ca: not_comparable / not_significant
- UNMAPPED EV-1c2640fda869f9ad297ad671: not_comparable / not_significant
- UNMAPPED EV-272c16667196fbb67cb329b2: not_comparable / not_significant
- UNMAPPED EV-a8f6c544a4d0cc1d74ed713c: not_comparable / not_significant
- UNMAPPED EV-e8297882c3bfe0b113ff5d32: not_comparable / not_significant
- UNMAPPED EV-ba5b2732f428da7cc4b6f217: not_comparable / not_significant
- UNMAPPED EV-4e9e8d2b5afb5a5116fe5d20: not_comparable / not_significant
- UNMAPPED EV-687f3bfcb1ac8908df5ec4e5: not_comparable / not_significant
- UNMAPPED EV-eae2e17cb17671d925b2b253: not_comparable / not_significant
- UNMAPPED EV-3b1bc1d563a36cea1e0415bc: not_comparable / not_significant
- UNMAPPED EV-5ffed98886952aee7cf2d336: not_comparable / not_significant
- UNMAPPED EV-5760e1cfe23308658bdb3cad: not_comparable / not_significant
- UNMAPPED EV-78b96561873c053758d4f29e: not_comparable / not_significant
- UNMAPPED EV-06d4947e76e837e809cef1cd: not_comparable / not_significant
- UNMAPPED EV-740d07f032dbdeac51b094a6: not_comparable / not_significant
- UNMAPPED EV-c9cdb88c79ecd75ae12a0a62: not_comparable / not_significant
- UNMAPPED EV-947505e345c37105e2afcc4c: not_comparable / not_significant
- UNMAPPED EV-22de8e658dfdeef8a6a175cd: not_comparable / not_significant
- UNMAPPED EV-40d2aadd924b569312ea1e61: not_comparable / not_significant
- UNMAPPED EV-e060db5ebc36a81cd5f394f1: not_comparable / not_significant
- UNMAPPED EV-76715901fe9129693467257d: not_comparable / not_significant
- UNMAPPED EV-08c04e01fb40f1ede45f4f9d: not_comparable / not_significant
- UNMAPPED EV-d119701c986808727bc79d5a: not_comparable / not_significant
- UNMAPPED EV-e98a1cbb0dc0ac3560111e46: not_comparable / not_significant
- UNMAPPED EV-e56786cc25093fccb39f37fe: not_comparable / not_significant
- UNMAPPED EV-810cc6ca1d028fbe53500320: not_comparable / not_significant
- UNMAPPED EV-03291c8a76390aa2f0a41959: not_comparable / not_significant
- UNMAPPED EV-e77ceff016f8978bbd1c38ea: not_comparable / not_significant
- UNMAPPED EV-a07bcde48a1539b3a3ac63b0: not_comparable / not_significant
- UNMAPPED EV-59358f88bf41ca42582e83b9: not_comparable / not_significant
- UNMAPPED EV-84a7277e96e28fb2bbe1f10d: not_comparable / not_significant
- UNMAPPED EV-f717f28ae37a9789caff9941: not_comparable / not_significant
- UNMAPPED EV-f34452e224416412d9ea0be0: not_comparable / not_significant
- UNMAPPED EV-c71c1b3c7eb1770a80b59fba: not_comparable / not_significant
- UNMAPPED EV-3aa3f569cf531781a8c661f7: not_comparable / not_significant
- UNMAPPED EV-66a85adbf0ac3bf04cd20d23: not_comparable / not_significant
- UNMAPPED EV-70695020d68f3311ab33fbe3: not_comparable / not_significant
- UNMAPPED EV-5e972c18e3d453b554f7bddb: not_comparable / not_significant
- UNMAPPED EV-6f9fa39ccd4a8d2e518733f6: not_comparable / not_significant
- UNMAPPED EV-0a6f778e7bd7c1a9d28408a3: not_comparable / not_significant
- UNMAPPED EV-195544bd5a3ea4e8125c3be0: not_comparable / not_significant
- UNMAPPED EV-3a07cd684c8f9fac2321a2f5: not_comparable / not_significant
- UNMAPPED EV-f648c7a8f6e64d482b462737: not_comparable / not_significant
- UNMAPPED EV-88f69d1efc4e789e90de6cf9: not_comparable / not_significant
- UNMAPPED EV-30a902069458c8fae62a2cc8: not_comparable / not_significant
- UNMAPPED EV-8fe9a41770c0989877dfd3d8: not_comparable / not_significant
- UNMAPPED EV-36287874fb0ddbf67d028151: not_comparable / not_significant
- UNMAPPED EV-f4237f3384dfedd5e095df58: not_comparable / not_significant
- UNMAPPED EV-1cfcd28d11ae5d14d38cf598: not_comparable / not_significant
- UNMAPPED EV-aed11a745f558f9a6e935e76: not_comparable / not_significant
- UNMAPPED EV-77b162aa29f0b02cbda41fdc: not_comparable / not_significant
- UNMAPPED EV-e8064ca1ef4a1a2f3cfe706b: not_comparable / not_significant
- UNMAPPED EV-8e85cacd8ce4f6b56752ef0a: not_comparable / not_significant
- UNMAPPED EV-5f8bc95c458130b00ee17b2f: not_comparable / not_significant
- UNMAPPED EV-917ab6906a2c007301ce00b6: not_comparable / not_significant
- UNMAPPED EV-64ba71b501d58537ed6a6dcb: not_comparable / not_significant
- UNMAPPED EV-452793ad44fd5c9a9dbf54bd: not_comparable / not_significant
- UNMAPPED EV-275c9c2e89891fb76c1f1860: not_comparable / not_significant
- UNMAPPED EV-99ce2e94ba6c6e472e3e5b3e: not_comparable / not_significant
- UNMAPPED EV-8420574f1dadc95f49f73370: not_comparable / not_significant
- UNMAPPED EV-b96855faa057422e0db33ceb: not_comparable / not_significant

## Per-provider contributor coverage

| Provider | complete | partial | unavailable | absent |
| --- | --- | --- | --- | --- |
| g:Profiler | 0 | 0 | 16455 | 244 |
| DAVID | 318 | 0 | 0 | 16381 |

- combined complete (shared term, two complete nonempty sets): 0

## Top clusters

| Cluster | Score | Rows | Mapped terms | Shared | Context genes | Context source |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.77 | 3 | 2 | 0 | BUB1B, MAD1L1, MAD2L2, NCAPD2, SMC4 | cluster_observed |
| 2 | 2.71 | 3 | 2 | 0 | BUB1B, MAD1L1, MAD2L2, NCAPD2, SMC4 | cluster_observed |
| 3 | 2.34 | 5 | 2 | 0 | none | cluster_observed |
| 4 | 1.92 | 3 | 3 | 0 | none | cluster_observed |
| 5 | 1.49 | 6 | 0 | 0 | none | cluster_observed |
| 6 | 1.4 | 6 | 1 | 0 | none | cluster_observed |
| 7 | 1.33 | 9 | 1 | 0 | none | cluster_observed |
| 8 | 1.27 | 5 | 0 | 0 | none | cluster_observed |
| 9 | 1.24 | 3 | 2 | 0 | none | cluster_observed |
| 10 | 1.11 | 5 | 1 | 0 | none | cluster_observed |
| 11 | 1.07 | 5 | 1 | 0 | none | cluster_observed |
| 12 | 1.07 | 3 | 2 | 0 | none | cluster_observed |
| 13 | 0.99 | 5 | 0 | 0 | none | cluster_observed |
| 14 | 0.98 | 3 | 3 | 1 | none | cluster_observed |
| 15 | 0.96 | 3 | 3 | 0 | none | cluster_observed |
| 16 | 0.96 | 3 | 0 | 0 | none | cluster_observed |
| 17 | 0.95 | 9 | 6 | 1 | none | cluster_observed |
| 18 | 0.93 | 22 | 14 | 1 | none | cluster_observed |
| 19 | 0.88 | 5 | 5 | 1 | none | cluster_observed |
| 20 | 0.85 | 6 | 5 | 3 | none | cluster_observed |

## Observed gene membership

- DAVID (1268 genes, 10717 term-gene rows): A1BG, A2M, AAMDC, ABCC1, ACADM, ACADS, ACADVL, ACBD6, ACCS, ACE, ACLY, ACOT13, ACOX2, ACP5, ACSL6, ACSM1, ACSM2B, ACSS3, ACTG2, ACTL7B, ACTN1, ACTN2, ADA, ADAMTS1, ADAMTS12, ADAMTS17, ADAMTS7, ADAMTSL1, ADARB1, ADCK5, ADCY8, ADGRA3, ADH1A, ADH4, ADK, AFAP1L1, AFG3L2, AGAP1, AGAP6, AGMO, AGXT, AGXT2, AIFM1, AIP, AJUBA, AK9, AKR1C2, AKR1C3, AKR1E2, ALDH16A1, ALDH1B1, ALDH3A1, ALDH3A2, ALDH3B2, ALDH7A1, ALG1, ALG11, ALG12, ALG14, ALG1L2, ALK, ALOX5AP, AMACR, ANKFY1, ANKLE2, ANKRD2, ANKRD23, ANKRD27, ANKRD36C, ANKRD44, ANTXR1, ANXA11, AP1AR, AP1S1, AP1S3, AP4S1, APOA1, APOA4, APOBEC3D, APOC1, APOL3, ARAP1, ARAP2, ARFRP1, ARG1, ARHGAP26, ARHGAP9, ARHGEF17, ARHGEF19, ARHGEF37, ARID1B, ARID5B, ARL9, ARMC2, ARMC6, ARMC9, ARMS2, ARPC3, ASAH1, ASGR1, ASIC3, ATE1, ATG3, ATL1, ATM, ATOX1, ATP1B3, ATP5MC1, ATP5MJ, ATP6V0A4, ATP6V1B1, ATP6V1B2, ATXN2L, AVPI1, BABAM2, BAIAP3, BANK1, BAZ1A, BBS4, BBS9, BBX, BCAR1, BCAS2, BCAS3, BCKDHA, BCKDHB, BCL2L13, BCR, BEGAIN, BEND7, BICDL2, BID, BIN1, BLVRA, BMP2K, BMS1, BNIP3, BPHL, BPI, BRAT1, BRD2, BTN2A1, BUB1B, C11ORF52, C12ORF60, C14ORF180, C15ORF40, C1GALT1, C1QL1, C1QL2, C22ORF23, C2CD2, C2CD3, C6ORF132, C9ORF78, CA1, CACNB1, CALD1, CALM2, CAMK1G, CAMK2A, CAMK2B, CANX, CAPN12, CARD14, CARMIL1, CASKIN2, CASP9, CCDC116, CCDC124, CCDC137, CCDC149, CCDC169, CCDC170, CCDC198, CCDC24, CCDC32, CCDC66, CCDC68, CCDC69, CCDC74A, CCHCR1, CCL24, CCN4, CCN6, CCNB1IP1, CCNL2, CCR5, CCT3, CD200R1L, CD27, CD300H, CD302, CD34, CD3E, CD55, CD82, CDH23, CDK1, CDKN1A, CELA1, CENPM, CEP131, CEP170, CEP76, CEP89, CERS4, CES1, CETP, CFAP263, CFAP44, CFL1, CFLAR, CFP, CGA, CGB7, CGREF1, CHCHD10, CHIA, CHMP4A, CHPT1, CHST8, CHTF18, CHUK, CIAO1, CIDEA, CIMIP2A, CIROP, CITED2, CKM, CLDN5, CLDN7, CLEC4C, CLMN, CLN5, CLYBL, CMBL, CMC1, CNEP1R1, CNOT1, CNOT10, CNOT9, CNTN2, COA5, COG3, COL11A2, COL18A1, COL4A2, COL6A2, COL6A3, COL6A5, COLEC12, COPRS, COX10, COX16, CPB2, CPED1, CPN2, CPSF6, CRACR2B, CRCP, CRNKL1, CRYAA, CRYAB, CSGALNACT1, CSPP1, CSRP3, CTDSP1, CTDSPL2, CTNNB1, CTSC, CTSL, CTSZ, CTU2, CYB5A, CYB5B, CYBA, CYC1, CYFIP2, CYLC2, CYP11B1, CYP21A2, CYP27A1, CYP2A6, CYP2A7, CYP2C19, CYP2C8, CYP2D6, CYP2J2, D2HGDH, DARS1, DBNL, DCAF10, DCAF15, DCD, DCHS2, DDIT4L, DDOST, DDX11, DDX24, DDX28, DDX47, DDX49, DDX60L, DEFB121, DENND2B, DGCR2, DGCR6, DGKE, DGKH, DHCR24, DHPS, DIDO1, DIP2A, DLG2, DNAAF4, DNAAF5, DNAH10, DNAH11, DNAH14, DNAJB6, DNAJC24, DNAJC7, DOCK2, DPF3, DPP7, DPPA2, DPPA4, DR1, DRAM2, DTX1, DTYMK, DUSP1, DYM, DYNC1LI1, DYNC2I2, DYNLT1, DZANK1, DZIP3, ECSIT, EEF1A2, EEF1E1, EFCAB13, EIF2S3, EIF4E, EIF4ENIF1, EIPR1, ELAC2, ELOB, ELP1, ELP3, ELP5, EMP2, ENAH, ENSG00000291072, ENSG00000291130, EPB41L4A, EPRS1, EPS15L1, ERCC6L2, ERGIC1, ERICH1, ETV2, ETV4, EVA1A, EVI5, EVPLL, EWSR1, EXD3, EXOC3, EXOC3L2, EXOC3L4, EXOC6B, EXTL3, F11, F13B, FAAH, FABP3, FAM111B, FAM131C, FAM221A, FAM222B, FAM228A, FAM241B, FAM3B, FAM3C, FAN1, FARP1, FARS2, FARSA, FAT4, FAXDC2, FBF1, FBXL4, FBXO42, FCGR2A, FCGR3B, FEZ1, FEZ2, FGA, FGD5, FGFRL1, FHL1, FHL5, FIZ1, FLG, FLT4, FMO3, FN3KRP, FNTA, FOCAD, FOXD4, FOXN4, FOXO1, FOXRED1, FRA10AC1, FRG2C, FSTL1, FTCD, FTO, FTSJ3, FXYD5, G3BP1, GABRB2, GABRG3, GALNT3, GALNT8, GALNTL5, GAPVD1, GARIN1A, GBA3, GBF1, GCNT3, GDI2, GDPD5, GEM, GEMIN2, GEMIN7, GEN1, GFM1, GFUS, GINS1, GIPC1, GIPC3, GKAP1, GLB1, GLDC, GLRB, GLT6D1, GLT8D1, GLYCTK, GMCL1, GMEB2, GNAI1, GNB1, GOLGA6L2, GOLGA6L6, GOLGA8S, GOLM1, GOSR1, GOSR2, GOT2, GPAM, GPAT4, GPATCH4, GPX7, GREB1, GRIN3B, GRM3, GRM7, GRTP1, GSDMB, GSE1, GTF2A1, GTF2IRD1, GTPBP10, GUCA1C, GUCY1B1, GYPB, HADH, HAS1, HAUS6, HDHD5, HDLBP, HEATR1, HHAT, HIBCH, HIF1A, HINT1, HLA-A, HLA-B, HLA-C, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5, HMBOX1, HMGB1, HNF1A, HNRNPA1L2, HNRNPAB, HOXD1, HPGD, HPN, HPRT1, HRNR, HSCB, HSD17B2, HSPA8, HYCC2, HYDIN, ICE2, IDH2, IDI2, IDO2, IFIT3, IFITM2, IFT25, IGBP1, IGDCC3, IGHMBP2, IGSF10, IL32, IL34, ILRUN, INF2, ING5, INO80, INPP5D, INPP5E, INTS4, IPO11, IQGAP1, IRAG1, IRAK1BP1, IRF2, IRF8, ISCU, ITGB1, ITGB2, ITGB3, ITIH1, ITPKB, ITSN1, JAG2, JPH2, JSRP1, KANSL1, KCNK17, KDM4C, KDM5B, KIAA0586, KIF11, KIF12, KIF16B, KIF1B, KISS1, KISS1R, KIT, KIZ, KLHL23, KLHL3, KLHL4, KMO, KMT2B, KMT2C, KRI1, KRT24, KRT31, KRTAP1-1, KRTAP1-5, KRTAP10-6, KRTAP13-2, KRTAP19-6, KRTAP4-8, KRTAP9-1, L2HGDH, L3MBTL4, LAMA3, LAMA5, LAMC1, LAPTM4B, LAPTM5, LARGE2, LCMT1, LCN10, LDB3, LFNG, LGALS8, LGI1, LHPP, LILRA1, LILRB3, LILRB4, LINC03126, LINS1, LMBRD1, LMBRD2, LOXL1, LPIN1, LRP4, LRPAP1, LRRC75A, LRRIQ1, LRRK2, LRSAM1, LRWD1, LSR, LTBP1, LTK, LYST, LZTR1, MAD1L1, MAD2L2, MAGI2, MAK16, MAML3, MAN2A1, MAP1S, MAP2K3, MAP3K21, MAP6, MAPK10, MAPRE2, MARCHF10, MARCHF7, MARK3, MATN2, MC3R, MCC, MCF2L, MCL1, MDM2, MDM4, MDP1, MED28, MEIG1, MEP1B, MFAP2, MFHAS1, MFSD11, MFSD12, MGST3, MIB2, MICA, MICU3, MIER3, MIPOL1, MITF, MLEC, MLYCD, MMAB, MMD2, MPHOSPH6, MPPE1, MPRIP, MROH2A, MRPL2, MRPL24, MRPL27, MRPL32, MRPS2, MRPS22, MRPS27, MRPS9, MRTFA, MS4A12, MS4A14, MS4A6A, MS4A7, MSH5, MSRB2, MST1L, MT1X, MTARC1, MTCL1, MTERF1, MTHFR, MTMR2, MTNAP1, MUC3A, MUC7, MVD, MXRA8, MYBPC1, MYCT1, MYH10, MYH9, MYL12B, MYL2, MYL7, MYO15A, MYO15B, MYO5A, MYO7A, MYOC, MYOF, MYRF, MYT1L, NAA16, NAALAD2, NAMPT, NANS, NAV3, NAXD, NBN, NBR1, NCAPD2, NCAPH2, NCL, NDE1, NDUFAF4, NDUFAF6, NDUFB1, NDUFB5, NDUFB9, NDUFS6, NECTIN2, NF1, NFKBIL1, NFX1, NFYC, NHERF2, NINJ1, NIPA2, NIT2, NKAPD1, NKRF, NLRC5, NLRP1, NOM1, NPDC1, NPHP1, NPIPB2, NPSR1, NPW, NPY, NREP, NRG4, NRXN1, NSL1, NT5C, NT5C2, NT5C3B, NTN1, NTRK1, NUCB1, NUCB2, NUDT1, NUDT11, NUDT14, NUP188, NUP35, NUP93, NUP98, OAS1, OAS2, OLFM1, OPRM1, OPTN, OR5K2, OR6C1, P2RX3, PAFAH1B1, PAK4, PALM, PARD3, PARP1, PARP10, PARP14, PAX4, PBLD, PBOV1, PBRM1, PBXIP1, PCCB, PCGF2, PCK1, PCSK4, PCTP, PDCD4, PDCD6, PDIA5, PDIA6, PDP1, PDS5A, PEAK3, PEMT, PER3, PERM1, PEX11G, PFDN5, PGLYRP2, PGPEP1L, PHC1, PHF3, PHKG1, PIGQ, PIK3C3, PINK1, PIP5K1A, PIP5KL1, PITPNM2, PITRM1, PJA2, PKD1, PKD1L2, PKHD1, PKIB, PLA2G2C, PLAUR, PLBD1, PLCG2, PLEK, PLEKHH2, PLGRKT, PML, PMS1, PNLIPRP2, PNN, PNPLA3, POLD2, POLDIP2, POLR2D, POLR3F, POLR3H, POLR3K, POT1, POU2F1, POU6F2, PPEF2, PPFIBP1, PPIL6, PPP1CA, PPP1R12B, PPP1R9B, PPP3R1, PPP4R3B, PRDM16, PRDM5, PRICKLE1, PRKAG2, PRKAR1A, PRKRA, PRMT6, PRNP, PRPF18, PRPH, PRPSAP1, PRSS1, PRSS22, PRSS48, PSG1, PSG11, PSG2, PSG5, PSG8, PSMA4, PSMD12, PSMD14, PSMD3, PSMG3, PSMG4, PSORS1C1, PTAR1, PTBP3, PTCD3, PTCHD3, PTGES3, PTGES3L, PTGR2, PTGR3, PTP4A3, PTPA, PTPMT1, PTPN2, PTPN4, PTPRJ, PXMP2, PYROXD1, PYROXD2, R3HCC1L, RAB3IP, RAB40C, RAB4B, RACK1, RAD17, RAD23B, RAF1, RAG1, RBFOX1, RBIS, RBKS, RBM48, RBM6, RBM7, RC3H2, RCN2, RD3, RDH13, RELN, RER1, RETSAT, REV3L, RGL4, RHBG, RHNO1, RHOBTB3, RIMS1, RIN2, RIPOR2, RMDN2, RNASE9, RNF11, RNF13, RNF130, RNF212, RNF40, RNF44, RP1, RP1L1, RPF1, RPF2, RPH3AL, RPL13, RPP30, RPS16, RPS17, RPS21, RPS29, RPS6KA1, RRAGC, RUNX2, S100B, S100P, SAA1, SCAF4, SCAPER, SCARB1, SCD5, SCIN, SCLY, SCN2A, SCN5A, SCRN3, SDHD, SDR39U1, SDR42E1, SEC23B, SEC63, SELENBP1, SEMA3C, SEMA6C, SEPTIN2, SERAC1, SERP1, SERPINA1, SERPINA2, SERPINB6, SET, SETBP1, SF3B2, SFMBT1, SFSWAP, SFTPA1, SFTPA2, SGF29, SH3YL1, SHLD1, SIGLEC10, SIGLEC12, SIM2, SIRT3, SIX1, SKP2, SLAIN1, SLC11A2, SLC16A12, SLC19A1, SLC1A5, SLC22A1, SLC22A24, SLC23A1, SLC25A10, SLC25A5, SLC29A1, SLC2A8, SLC30A4, SLC34A1, SLC35B1, SLC37A1, SLC37A4, SLC38A10, SLC38A9, SLC3A1, SLC3A2, SLC41A3, SLC6A12, SLC8A1, SLC9B2, SMARCD2, SMC1B, SMC4, SMIM26, SMPDL3B, SNRNP35, SNRNP40, SNRPB, SNRPE, SNX13, SNX15, SNX29, SNX7, SOCS5, SON, SORBS3, SP100, SPATA22, SPATA31F3, SPATA7, SPDEF, SPEN, SPG7, SPIDR, SPINT1, SPP1, SPTA1, SPTAN1, SRA1, SRL, SSBP4, ST3GAL3, STAM2, STAP2, STARD3NL, STAU1, STAU2, STC2, STIMATE, STK17B, STK32B, STMN1, STOX1, STXBP4, STYK1, SUCLA2, SUMO4, SUN1, SUPT3H, SVIL, SYNE2, SYNE3, SYT1, TACC2, TAF9, TAGLN3, TANC1, TAOK1, TAPBPL, TAS2R19, TBC1D14, TBC1D31, TBC1D7, TC2N, TCEAL6, TDG, TECR, TEX26, TFR2, TGOLN2, THAP6, THUMPD3, TIE1, TIGD6, TLCD3A, TLE2, TLR8, TM4SF4, TMC6, TMC8, TMCC2, TMEM106B, TMEM107, TMEM120A, TMEM165, TMEM176A, TMEM176B, TMEM182, TMEM186, TMEM216, TMEM25, TMEM254, TMEM255B, TMEM259, TMEM276, TMEM63A, TMEM67, TMEM87A, TMEM8B, TMEM91, TMPRSS2, TMPRSS3, TMPRSS4, TNC, TNFRSF10B, TNFRSF10C, TNIK, TNIP1, TNK1, TNNT2, TOM1L2, TOP1MT, TOR1AIP1, TOR3A, TOX4, TP53, TP53I11, TP53I13, TPD52L2, TPGS2, TPK1, TPM2, TPSB2, TPT1, TPTE2, TRDMT1, TRIM4, TRIM5, TRIM51G, TRIP4, TRMT61A, TRPC4, TRPV3, TSEN54, TSFM, TSPAN16, TSSK6, TTBK2, TTC12, TTC27, TTC7A, TTLL12, TTLL3, TTLL5, TTN, TTR, TTYH1, TUBA8, TUBB1, TUBD1, TUBGCP2, TXN2, TYR, TYSND1, U2AF1L4, UBAC2, UBE2D3, UBE2D4, UBE2F, UBE2I, UBE2QL1, UBE2Z, UBR4, UCHL5, UEVLD, UFM1, UGP2, UGT2B7, UIMC1, UMAD1, UNC13B, UNC50, UPB1, UPF3A, UPP1, UQCC2, URGCP, USP29, USP45, USP48, VAMP2, VARS2, VAT1, VDR, VEZT, VGLL3, VHL, VIPAS39, VLDLR, VOPP1, VPS13B, VPS35L, VPS41, VPS45, VPS51, VWDE, WAS, WBP11, WDHD1, WDR18, WDR37, WDR41, WDR7, WDR73, WFS1, WIF1, WNK1, WWOX, WWP2, WWTR1, XCR1, YBX3, ZBTB10, ZDHHC1, ZDHHC11B, ZFPM1, ZFR2, ZFYVE26, ZIM3, ZKSCAN3, ZMAT1, ZMPSTE24, ZMYM5, ZNF117, ZNF133, ZNF146, ZNF148, ZNF197, ZNF254, ZNF3, ZNF304, ZNF329, ZNF331, ZNF33A, ZNF417, ZNF44, ZNF480, ZNF519, ZNF544, ZNF552, ZNF566, ZNF567, ZNF573, ZNF586, ZNF594, ZNF598, ZNF599, ZNF625, ZNF626, ZNF66, ZNF677, ZNF681, ZNF7, ZNF714, ZNF717, ZNF76, ZNF765, ZNF778, ZNF8, ZNF80, ZNF805, ZNF84, ZNF880, ZNF98, ZPBP, ZRSR2, ZSCAN2, ZSCAN9

## Open review rows: 113

- david_only GO GO:0005515: mapped on a stable identifier; the two providers differ in coverage, effective background, thresholds or mappings
- david_only GO GO:0005829: mapped on a stable identifier; the two providers differ in coverage, effective background, thresholds or mappings
- unmapped_term WIKIPATHWAYS nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists
- unmapped_term UNMAPPED nan: unmapped: the provider term carries no supported stable identifier, so no counterpart state exists

## Output links

- [integrated_pathway_evidence.csv](./integrated_pathway_evidence.csv)
- [gene_pathway_support.csv](./gene_pathway_support.csv)
- [david_cluster_summary.csv](./david_cluster_summary.csv)
- [provider_discordance_review.csv](./provider_discordance_review.csv)
- [integrated_pathway_report.md](./integrated_pathway_report.md)
- [postanalysis_provenance.json](./postanalysis_provenance.json)

## Limitations

- Provider statistics remain separate; no pooled hypothesis test was calculated.
- Authoritative public definitions were not retrieved in this offline run, so term interpretation stays provisional.
- File hashes for every generated artifact are recorded in `postanalysis_provenance.json`, which is excluded from its own hash map by convention.
- Null Jaccard denotes unavailable complete contributor sets, not disjoint gene membership.
- Gene-set membership is observed annotation context, not pathogenicity evidence.
- A nonspecific term for an MVA mechanism is context; modification terms annotate membership, not measured activity.
