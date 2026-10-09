#!/usr/bin/env Rscript
# ==============================================================================
# Module scores of ASCL1 target gene sets in Ascl1+ vs Ascl1- hepatocytes
# ==============================================================================
#
# Description:
#   Pre-geriatric hepatocytes are split into Ascl1+ (normalized Ascl1 > 0) and
#   Ascl1- cells. For each gene set, a Seurat module score is compared between
#   Ascl1+ and Ascl1- cells within each sex:
#     - animal level: mean score per mouse in Ascl1+ and Ascl1- cells (mice
#       with >= 10 cells in both groups), paired Wilcoxon and Cohen's dz
#     - within-animal label permutation: Ascl1 status shuffled inside each
#       mouse 10,000 times, P = (r + 1) / (n + 1)
#
#   Figure 7, panels i-l (box + jitter; stars from the permutation P)
#     i  reactome_female_promoter
#     j  reactome_P2G
#     k  P2G_C4
#     l  P2G_C2
#
# Input:
#   - rna_wnn.h5ad (from 04_multiome_integration/02_wnn_integration.py)
#     with cell_type, age, sex and sample
#
# Output (in OUTPUT_DIR, per gene set <tag>):
#   - Ascl1_module_score_<tag>_stats.csv        animal-level statistics
#   - Ascl1_module_score_<tag>_permutation.csv  label permutation test
#   - Ascl1_module_score_<tag>_box_jitter.pdf   box + jitter (figure panel)
#   - Ascl1_module_score_<tag>_ALLPANELS.pdf    violin, box, box + jitter,
#                                               density, paired mice, null
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(schard)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
H5AD_PATH  <- "rna_wnn.h5ad"
OUTPUT_DIR <- "module_score"
N_PERM     <- 10000

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# GENE SETS
# ==============================================================================
GENE_SETS <- list(

  # panel i: genes in Reactome pathways enriched among female promoter peaks
  reactome_female_promoter = c(
    "Gldc", "Eno1", "Ndst1", "As3mt", "Cyp2d6", "Nampt", "Psmd3", "Akt1",
    "Mlycd", "Mmut", "Prkaca", "Gamt", "Stard4", "Acot7", "Stard5", "Nup210",
    "Phka1", "Prkca", "Armc8", "Upb1", "Rufy1", "Slc5a6", "Rbp4", "Plbd1",
    "Cyp2e1", "Tmem126b", "Hpd", "Uba52", "Acot4", "Mtmr3", "Seh1l", "Shmt2",
    "Shmt1", "Iqgap1", "Adh7", "Hacd3", "Mtmr4", "Pla2g6", "Sdsl", "Gns",
    "Adh5", "Mtmr6", "Hsd11b1", "Adh4", "Ldha", "Cyb5r3", "Nup85", "Agmo",
    "Prkar2a", "Ttpa", "St3gal4", "Fdft1", "Lhpp", "Slc19a2", "St3gal3",
    "Crebbp", "Gstm1", "Chka", "Cyp4f2", "Stard10", "Sqle", "Bhmt", "Tcn2",
    "Adi1", "Pah", "Ogdh", "Mms19", "Gpi", "Mtmr10", "Qprt", "Adk", "Prkag2",
    "Pygl", "Csad", "Serpina6", "Med15", "Impa1", "Tk1", "Elovl1", "Sphk2",
    "Coq10b", "Coq10a", "Pgm2l1", "Ncor2", "Itpkc", "Plcb3", "Bdh2", "Med24",
    "Mtf1", "Agps", "Pgam5", "Plin2", "Angptl4", "Ppara", "Tkt", "Mgll",
    "Ppard", "Rapgef4", "Asah1", "Samd8", "Sgms1", "Idua", "Gstp1", "Odc1",
    "Adcy6", "Gpat4", "Ldlr", "Bche", "Pklr", "Aadat", "Gss", "Mocos", "Bbox1",
    "Fmo1", "G0s2", "Gnmt", "Mkln1", "Lyrm4", "Etnk2", "Lpcat3", "Ciapin1",
    "Crat", "Aldh1l1", "Rora", "Fabp12", "Fads2", "Iyd", "Urod", "Mpc2",
    "Slco2b1", "Osbp", "Nudt19", "Slc25a44", "Dgat2", "Dgat1", "Gpt2", "Dio1",
    "Ipmk", "Pipox", "Tpst2", "Slc27a2", "Pfkfb2", "Slc26a1", "Abhd3", "Insig2",
    "Itpr1", "Itpr2", "Mri1", "Ak4", "Acat2", "Tkfc", "Esd", "Them4", "Hacl1",
    "Nadk", "Esrra", "Osbpl9", "Osbpl8", "Slc37a1", "Got1", "Uros", "Got2",
    "Dhcr24", "Cyp4f11", "Cyp4f12", "Uckl1", "Mlxipl", "Alb", "Aco1", "Dhcr7",
    "Osbpl1a", "Alas2", "Alas1", "Hsp90ab1", "Etfa", "Alad", "Nadk2", "Gm2a",
    "Sin3a", "Acp5", "Upp2", "Idh3a", "Cers6", "Gpx1", "Apoa4", "Apoa5",
    "Coasy", "Inmt", "Cyp2a7", "Dusp23", "Cyp2a6", "Slc25a19", "Slc25a10",
    "Agxt", "Aspa", "Gapdh", "Cers2", "Slc25a13", "St6galnac6", "Ndufb9",
    "Acss3", "Carnmt1", "Lrp1", "Ppm1l", "Rpel1", "Gls2", "Rpe", "Ppm1k",
    "Cyp3a4", "Bphl", "Cyp3a5", "Cyp3a7", "Ndufa7", "Mmab", "Mid1ip1", "Apoc3",
    "Ddhd1", "Ass1", "Gstz1", "Pc", "Sucla2", "Uroc1", "Vapb", "Apoc2", "Cth",
    "Nags", "Faah", "Pank3", "Acsm3", "Acsm1", "Eci2", "Nudt5", "Hs6st1",
    "Rps14", "Pnp", "Ppat", "Pdk4", "Sephs2", "Pdk2", "Rps13", "Pdk1", "Mccc2",
    "Cpt1a", "Acsl1", "Phykpl", "Entpd5", "Sacm1l", "Mccc1", "Acsl5", "Slc6a12",
    "Taldo1", "Acsl4", "Cmbl", "Pisd", "Cyp27a1", "Cox7a2l", "Por", "Ppp1r3c",
    "Suclg2", "Ces2", "Ces1", "Pla1a", "Ptges", "Gcdh", "Sdc4", "Mgst1",
    "Pik3r1", "Pgs1", "Inpp5a", "Ppcs", "Inpp5e", "Scap", "Apoe", "Arsg",
    "Decr2", "Pck1", "Arsb", "Surf1", "Parp16", "Aimp1", "Ext1", "Parp10",
    "Pm20d1", "Cyp26a1", "Hal", "Pemt", "Kyat1", "Acy3", "Acy1", "Hibadh",
    "Slc44a1", "Gbe1", "Tecr", "Gys2", "AcadL", "Hyal2", "Nmrk1", "Tnfaip8l1",
    "Acadm", "Dlat", "Aass", "Vkorc1", "Ncoa2", "Srebf1", "Manba", "Bhmt2",
    "Gaa", "Acad9", "Ampd2", "Aprt", "Nme1", "Ahcy", "Adh1c", "Adh1b", "Adh1a",
    "Haao", "Akr1d1", "Hsd17b4", "Cyp7a1", "Cyp17a1", "Nt5e", "Rxra", "Cyp2a13",
    "Hsd17b2", "Xdh", "Pnpla7", "Nfyb", "Pnpla8", "Gcgr", "Plekha5", "Fah",
    "Gck", "Etnppl", "Ppp1ca", "Aldh4a1", "Nudt9", "Gsta5", "Pdp2", "Gpam",
    "Pnpla2", "Acaa2", "Inppl1", "Cpox", "Msmo1", "Comt", "Fpgs", "Acaa1",
    "Sbf1", "Ca14", "Phyh", "Miga2", "Hgd", "Acot11", "Acot12", "Mat1a",
    "Aldh3a2", "Pla2g15", "Tbl1xr1", "Gpd2", "Gpd1", "Cmpk1", "Idi1", "Nnmt",
    "Cacna1d", "Glyat", "Higd1a", "Pld1", "Acaca", "Pld2", "Ubc", "Prodh",
    "Mpst", "Gckr", "Ephx2", "Ephx1", "Crls1", "Azin1", "Qdpr", "Grhpr",
    "Sepsecs", "Amacr", "Gnb1", "Cox20", "Galk1", "B4galt1", "Chd9", "Cyp4a22",
    "Slc2a2", "Dbi", "Pctp", "Aldh2", "Ca5a", "Naprt", "Aox1", "Man2c1",
    "Bckdha", "Hs3st3b1", "Abcc2", "Cyp3a43", "Cyp4a11", "Gpcpd1", "Prodh2",
    "Ndufs8", "Acox2", "Acox1", "Chdh", "Ehhadh", "Ndufs5", "Aldh1a1", "Hagh",
    "B4galt7", "Fbp1", "B4galt5", "Vac14", "Ddc", "Hpgd", "Pon2", "Hmgcr",
    "Nt5c2", "Agpat2", "Agpat3", "Ddo", "Rmnd5b", "Gna14", "Hmgcl", "Rmnd5a",
    "Cbs", "Mbtps1", "Hsd3b2", "Hsd3b1", "Suox", "Ttc19", "Pnkd", "Sardh",
    "Chpt1", "Thrsp", "Scly", "Lpin1", "Lpin2", "Ran"
  ),

  # panel j: genes in Reactome pathways enriched among female peaks overlapping hepatocyte P2G links
  reactome_P2G = c(
    "Dmgdh", "Slc23a2", "Cerk", "Gldc", "Pank1", "Eci2", "Ahr", "Nudt4",
    "Fads2", "Scp2", "Nampt", "Nudt19", "Pip4k2b", "Slc25a42", "Fads1",
    "Sult2a1", "Cpt1a", "Sds", "Hgd", "Acsl1", "Acot11", "Dio1", "Acsl3",
    "Mat1a", "Aldh3a2", "Cyp27a1", "Nup93", "Sult2b1", "Por", "Gpd2",
    "Mapkapk2", "Acot2", "Acot1", "Plbd1", "Trib3", "Cyp2e1", "Slc27a2",
    "Acot4", "Ces2", "Stx1a", "Ces1", "Pla1a", "Nnmt", "Slc22a5", "Sdc4",
    "Insig2", "Sar1b", "Insig1", "Ak2", "Gstt2", "Ak4", "Pla2g6", "Pld1",
    "Sdsl", "Aldh3b2", "Amdhd1", "Agmo", "Ubc", "St3gal5", "Asl", "St3gal6",
    "Arsg", "Decr2", "Pck1", "Prodh", "Gstt2b", "St3gal1", "Lhpp", "Chkb",
    "Slc37a1", "Slc10a1", "Nr1d1", "Azin1", "Acer1", "Cyp26a1", "Sqle", "Hal",
    "Uck1", "Cps1", "Pah", "Kyat3", "Sdc1", "Serinc3", "Gpi", "Alas1", "Tat",
    "Slc2a2", "Pik3cb", "Serpina6", "Gphn", "Gys2", "Aox1", "Upp2", "Cd36",
    "Sult4a1", "Aass", "Pla2g12a", "Abcc2", "Arg1", "Elovl5", "Acad9",
    "Cyp3a43", "Elovl6", "Adra2c", "Cyp7b1", "Gpcpd1", "Coq10b", "Prodh2",
    "Inmt", "Cyp2a7", "Cyp2a6", "Acox1", "Aldh1a1", "Aspg", "Pitpnm2", "Mvd",
    "Aldob", "Agxt", "Ppara", "Aspa", "B4galt5", "Ppard", "Rapgef4", "Slc25a13",
    "Acss3", "Ahcy", "Oat", "Adh1c", "Sgms1", "Adh1b", "Gls2", "Adh1a", "Pdhb",
    "Cyp3a4", "Cyp3a5", "Cyp7a1", "Adcy6", "Papss2", "Cyp3a7", "Gna14", "Nt5e",
    "Aldh1b1", "Cyp2a13", "Gpat4", "Ucp2", "Gpat3", "Cyp2f1", "Cyb5a", "Fahd1",
    "Sec24a", "Pklr", "Plekha1", "Nr1h4", "G0s2", "Cyp8b1", "Ass1", "Etnppl",
    "Gnmt", "Gsta5", "Uroc1", "Etnk2", "Cth", "Lpcat3", "Pnpla3", "Sardh",
    "Pnpla5", "Thrsp", "Lpin2"
  ),

  # panel k: hepatocyte P2G cluster C4
  P2G_C4 = c(
    "Ascl1", "Cux2", "4933404O12Rik", "Cyp3a41a", "Efnb2", "St3gal6", "Cyp3a44",
    "Cyp3a16", "Scn5a", "Lhpp", "Abcc4", "Zc2hc1a", "Sult2a1", "Slc1a2",
    "Pnpla3", "Aldh3a2", "Sult2a3", "Fmo3", "Sult2a2", "Sult2a4", "Cyp2a4",
    "Nipal1", "Gm13872", "Rap1gap2", "Lgr5", "Nrp2", "Cyp2e1", "Tmem64", "Glul",
    "Rcan2", "Cap2", "Mme", "Acot3", "Cyp2b9", "Oat", "Adcy1", "Slc22a27",
    "Prex2", "Colgalt2", "Adra2c", "Cyp2b13", "Slc13a3", "Tspan7", "Cyp2b10",
    "Tbx3", "Bmp7", "Acot4", "Cyp2c38", "Bnc2", "Tgfbr3", "A1bg", "Smad6",
    "Nt5e", "Pdzrn3", "Slc22a26", "Gulo", "Cyp2c40", "Setbp1", "Rmdn2", "Hamp2",
    "Myo9a", "Nkd1", "Slc24a3", "Cyp2a5", "Rorc", "Npr2", "Nxpe2", "Aacs",
    "Tox", "Sult3a1", "Tsix", "Pnpla5", "Ces1g", "Cyp2a22", "Cyp2c69",
    "1810008I18Rik", "Id4", "Cnst", "Vldlr", "Cd36", "Ucp2", "Lifr", "Srgap3",
    "Cyp2c39", "Greb1l", "Ahr", "Cyp2c37", "Pltp", "Slc16a10", "Cyp2c29",
    "Kif26b", "Col23a1", "Prodh", "Xist", "Slc22a3", "Rgs16", "Plbd1", "Pparg",
    "Serpina6", "Slco1b2", "Adh1", "Insrr", "Per3", "Ccnd1", "Airn", "Gna14",
    "Agtpbp1", "Insig1", "Cavin2", "Papss2", "Eci3", "Cyp3a25", "Clec2d",
    "Prom1", "Pik3cb", "Slc6a9", "Sorcs3", "A930009A15Rik", "Cyp7a1", "Usp2",
    "Rapgef5", "Tnfrsf19", "Sc5d", "Tdrd12", "Peg3", "Dbp", "Ginm1", "Sgms1",
    "Ces1b", "Ski", "Rnf43", "Pced1b", "Inhbb", "Cidea", "Esrrg", "Fmo2",
    "Fam126a", "Sult2a5", "Cyp1a2", "A930033H14Rik", "Acaa1b", "Sorcs2",
    "Endod1", "Vmn1r184", "Ces1d", "Ppard", "Adgrg2", "Stk33", "Sybu", "Mup19",
    "Clstn3", "Mroh9", "Gas6", "Lrtm1", "Sult4a1", "Il7", "Nxpe3", "Sytl5",
    "Teddm2", "Nxpe4", "Snrk", "Angptl8", "Gpx3", "Tsku", "Abcb1a", "Mpc2",
    "Acss3", "Pag1", "Susd4", "Synpo", "Mvd", "Gba2", "Tomm40l", "Slc16a11",
    "Sall1", "Parvb", "Fads2", "Spo11", "Cyp3a59", "Socs2", "Dclk3", "Stambpl1",
    "Adrb2", "Igf2r", "Inmt", "Tlr5", "Gsta3", "Slc22a2", "1810058I24Rik",
    "Notum", "Smox", "Rnase4", "Rfx4", "Slc40a1", "Fam102a", "Uck1", "Hmgcr",
    "Pgm3", "Fads1", "Pax8", "Slc16a5", "Zfp292", "Fut8", "Pex11a", "Tmem189",
    "Rdh16", "Hao2", "Abhd6", "Trim24", "Ube2e2", "Pon1", "Slc22a1",
    "2410131K14Rik", "Nlrp9c", "Cyp17a1", "Me1", "Atp10a", "Mapre2", "Dpf3",
    "Akr1c20", "Slc7a15", "Nhlrc2", "Cyp2c68", "0610005C13Rik", "Axin2", "Sp5",
    "Hacl1", "Fmo1", "B430212C06Rik", "Ccdc157", "Snai3", "Fads3", "Hspa12a",
    "Serpina16", "Crybg2", "Phldb2", "Slc23a2", "Pcp4l1", "Fitm1", "Nr1i3",
    "Nat8f1", "1700110K17Rik", "Srebf1", "Per2", "Id3", "Sord", "Gstt3",
    "Rab43", "Slc16a13", "Hip1r", "Fkbp4", "Gm16063", "Slc22a29", "Ces2c",
    "Fasn", "Ntrk2", "Cerk", "Gbp3", "Slc16a7", "Igfbp5", "Slc47a1", "Sbk1",
    "Ngef", "Usp10", "Rnase2a", "Kif1b", "Adcy10", "Dio1", "Nr1d1", "Car3",
    "Prlr", "Prodh2", "Abhd2", "Cyp4a14", "Klc1", "Hamp", "Aatk", "Akr1c6",
    "Khdrbs3", "Prex1", "Sh2d4a", "C030013G03Rik", "Dhdh", "Lingo4", "F7",
    "Dop1b", "Ces1e", "Dazap1", "Acot11", "Nrn1", "C630043F03Rik", "Mup15",
    "E330021D16Rik", "Hsf2", "Cdk6", "Stx1a", "Tbcel", "Rora", "Pitpnm2",
    "Sugp1", "Igfbp1", "Cyp2c50", "Cited2", "Tent5a", "Paqr9", "Zfp532",
    "Ces1c", "Ces1f", "Chic1", "Hs3st3a1", "Acer1", "Adra1b", "Mertk", "Dcaf6",
    "Tmem131l", "Smco4", "Miga2", "Ugt1a10", "Dnase2a", "1700007F19Rik",
    "Slc25a13", "Dio3os", "Cblb", "Ang", "Prelid2", "Morc3", "Id1", "Tnik",
    "Serpinb1a", "Gk5", "Aldh1a7", "Gypc", "Dct", "Wipf3", "Fgf1", "Nxph1",
    "Slc1a4", "S1pr5", "Acly", "4930523C07Rik", "Cipc", "H2-Q1", "Cyp2g1",
    "Esr1", "Blvrb", "Epm2aip1", "Rgn", "Bbs7", "Marveld1", "Pde1a", "Nampt",
    "Cyp27a1", "C730002L08Rik", "N4bp2l1", "Acot1", "Rgp1", "9330175M20Rik",
    "Tapt1", "Cyb5b", "Ankrd46", "Gm1110", "4931406C07Rik", "Cyp2c54",
    "Cyp4a10", "Zfp773", "Ube2d2a", "Bcl2l11", "Gstt2", "Ccdc85c", "Efr3b",
    "Dsg2", "Ccdc24", "Pip4k2b", "Mras", "Mcc", "Pcx", "Slc31a1", "Cdyl2",
    "Them4", "Pfkfb4", "Rnf186", "Ugt3a1", "Lipt2", "Dleu2", "Vwf", "Aoc1",
    "Slc1a5", "Irs2", "Abcd2", "Avpr1a", "A1cf", "Hif1a", "Kyat3", "Tjp1",
    "Klhl32", "Ar", "AI661453", "Ttc39b", "Grk3", "Por", "9430038I01Rik",
    "Vtcn1", "Aldh3b3", "Aqp8", "Tango2", "Aspscr1", "Gimap8", "Hsd3b7",
    "Ppp1r3b", "Fam49a", "Creb3l2", "Tmem79", "Sgsm1", "Cyp2c67",
    "C230024C17Rik", "Nedd9", "Lncppara", "Arrdc2", "Tns4", "Tnfsf15", "Zfp423",
    "Themis", "Lect2", "Gm19522", "Neu2", "Orm3", "Gse1", "2810013P06Rik",
    "Hmgcs1", "2310001H17Rik", "Micall2", "Tmem97", "Bcl2l1", "Enpp1",
    "Dclre1a", "Rhbg", "Pla2g15", "Tert", "Pycrl", "Helz", "Sult2a7", "S1pr1",
    "Mllt3", "Bag3", "Ttc23", "Pygo1", "Gbp7", "Gimap9", "Fam174a", "Fdps",
    "Camkk2", "Ccr6", "Trmt9b", "Gramd1b", "Rgs12", "Adamts17", "Ism2",
    "Gm2788", "Gimap1", "Spry4", "Stx1b", "Tiparp", "Mcf2l", "Aldh1a1", "Bmper",
    "Tmem25", "Gpam", "Cbx5", "Med12l", "Cobll1", "Tmem184b", "Rarb", "Gm20735",
    "Cldn2", "Bcat2", "Urgcp", "Ddx6", "Sardh", "G0s2", "Sult3a2", "Nav2",
    "Nuak2", "Zfp395", "Dlat", "Bbox1", "Irx1", "Sft2d3", "Gm5122", "Kcnk5",
    "2210417A02Rik", "Dntt", "Gata6", "Trim28", "1010001N08Rik", "Acsl3",
    "Cxxc5", "Ctse", "Afp", "Jup", "Pemt", "Rap2a", "Msl2", "Kank1", "Acad9",
    "Acsl5", "Olfr45", "Mast3", "Myo19", "P2ry14", "Ppara", "Dap", "Gls2",
    "Elovl6", "Hpd", "Ccdc69", "Egln3", "Tfpi2", "Sgk3", "2300009A05Rik",
    "Pknox2", "Slc16a6", "Tmem14a", "Elmo1", "Lipg", "Slc36a1", "Pxdc1",
    "Tmprss4", "Myo10", "Sned1", "Adcy6", "D930048N14Rik", "Plin2", "Srebf2",
    "Cyp2j11", "Adamts2", "Atp2b2", "Hnrnpc", "Cdc42ep4", "Acnat2", "Olfr541",
    "Csad", "Rnf166", "BC049352", "Slbp", "Tcf24", "Fahd1", "D630039A03Rik",
    "Retreg3", "Gstt1", "Angpt1", "Evc2", "Aox1", "Dusp6", "Acacb", "Lmna",
    "Hsd17b4", "Unc79", "Dym", "Znrf3", "Celsr1", "Me3", "Nr2f2", "Mycn",
    "Wdr93", "Sipa1l2", "Kank2", "Alg11", "Pfkfb3", "Aadac", "Asap2",
    "5730420D15Rik", "Fam107b", "Hes6", "Mindy4", "Cers6", "Rnf24", "Ppp2r5a",
    "Adam10", "Adnp", "Mfsd6", "Fam117a", "Rcor3", "LTO1", "Smarca2", "Vnn1",
    "Slc25a21", "Plpp3", "Ugt1a6a", "Noct", "Gm9949", "E2f8", "Dhrs4", "Notch1",
    "Tceanc2", "Arrdc4", "Nceh1", "9030616G12Rik", "Ptpn9", "4632404H12Rik",
    "D630033O11Rik", "Dnah11", "Megf11", "Gm11437", "Apol9a", "Ftx", "Plgrkt",
    "Gata4", "Fmo4", "Nptx1", "Jpx", "Dsc2", "Ccdc62", "Rsph4a", "Gas1",
    "Chpt1", "Abcd4", "Cep152", "D1Ertd622e", "Paip2", "Tmc6", "Cmbl", "Mia2",
    "Dlec1", "Ctdspl", "Atl3", "Grid1", "Cga", "Mylk"
  ),

  # panel l: hepatocyte P2G cluster C2
  P2G_C2 = c(
    "Cyp2d9", "Scara5", "Susd4", "Fgfr1", "1700042O10Rik", "C6", "Serpina1e",
    "Serpina9", "Lama3", "Slc8a1", "Cyp2f2", "Ppp1r9a", "Dpy19l1", "Cdh1",
    "Hsd3b5", "Nudt7", "C4a", "Cyp7b1", "Ttc39c", "Ncam2", "Dpy19l3", "Gm16551",
    "Znrf2", "Cyp4a12a", "Cyp4a12b", "Selenbp2", "Serpina12", "Ddc", "Cabyr",
    "Aadat", "Ces4a", "Zbtb7c", "Phlda1", "Serpina3k", "Vwc2l", "Aspg", "Scp2",
    "Rassf3", "Pknox2", "Zyg11a", "Gstp1", "Hao1", "Serpina11", "Zfp445",
    "Tcf7l1", "Gm12718", "Nat8", "Cryl1", "Egfr", "Cyp21a1", "Gadd45g",
    "Clec2d", "Slco1a1", "Tmem19", "Slc25a30", "Olfm3", "Clec2h", "Etfbkmt",
    "Cmah", "Gulp1", "Reck", "Phf20l1", "Adgrf1", "Snx29", "C8a", "Ociad2",
    "Pard3b", "Tmx4", "Capn8", "Rbbp4", "C8b", "Onecut1", "Il1r1", "Mcm10",
    "Cyp2u1", "Ahcyl2", "Ces3b", "Hal", "Keg1", "Mindy3", "F2r", "Cib3",
    "Asap2", "Pklr", "Nek6", "Gfod2", "Amdhd1", "Ahnak", "Rxfp1", "Pdcd4",
    "Cpb2", "Igf2bp2", "Tcaim", "Lifr", "Hhex", "Platr4", "Map4k5", "Ikzf4",
    "Zc3h13", "Eps8l2", "Mup21", "Stk19", "Sfxn1", "Grb10", "Cish", "Dock1",
    "Gls2", "Nat8f6", "Aox3", "Apobec1", "Adam22", "Acadsb", "Ugt2b1", "Abcg2",
    "Cyp4a32", "Asap3", "Bcl6", "Saa4", "E2f8", "1810064F22Rik", "C9", "Arsa",
    "Bdnf", "Slc3a1", "Hcn3", "Slc22a28", "Tbc1d4", "Sult2a8", "Mup20", "Gpc1",
    "Sgce", "Ston1", "Il1rap", "2610035D17Rik", "Taf2", "Wdr37", "Moxd1",
    "Foxa2", "Mup7", "Fabp5", "Crtc3", "Ldah", "Col6a6", "Grm8", "Acsl1", "Clu",
    "Hopx", "Gldc", "Comt", "Spaca6", "Zfp809", "Rrbp1", "Cpne8", "Lasp1",
    "Fam124a", "Ugt2b38", "Nhsl1", "Phf8", "Galc", "Gm2788", "Mup3", "Myef2",
    "2810459M11Rik", "Fga", "1110059G10Rik", "Klhl2", "Serpina3m", "Tedc2",
    "Smagp", "Tlr5", "B3galt1", "C2", "Gtpbp4", "Camkk2", "Gm30551", "Adck5",
    "Etnppl", "Nrp1", "Cobl", "Aqp11", "Serpine2", "Gas2", "Cxadr", "Nox4",
    "Tbc1d30", "Hsd17b13", "Rbpms2", "Foxp2", "2310001H17Rik", "Saa2", "Acox3",
    "Myo6", "4930512B01Rik", "Dcakd", "Ugt2b35", "Arhgap32", "Cyp4f17",
    "Gm16159", "Tiam2", "Slc35e3", "Dtnb", "Sdr9c7", "Tspan33", "Serpina1c",
    "Tmem51", "Abhd5", "Mreg", "Qsox1", "Ephx2", "Tmem30a", "G0s2", "Adamts7",
    "Fmo5", "St13", "Zfhx2", "Hp", "Dop1b", "Ugt3a1", "2610507B11Rik", "Tuft1",
    "Prpf6", "Hdac1", "Fgb", "Zkscan7", "Foxa3", "Fgg", "Gm29684", "Saa1",
    "Ces3a", "Nrep", "Dalrd3", "Nsd3", "Kcnh7", "Col5a3", "Rlf", "Jph1", "Vmp1",
    "Rarres1", "Kif13b", "Zfp207", "Nat8f5", "4930452B06Rik", "Nufip2", "Edem3",
    "Bmyc", "Tshz3", "Rai14", "Adam2", "Nat8f1", "Gclm", "Cadm4", "Arhgap42",
    "Glyat", "Pax8", "Dennd1b", "Epb41l4b", "C230037L18Rik", "Synrg", "Xbp1",
    "Tmem120b", "Lncppara", "Senp6", "Atf5", "Dclk3", "Wdtc1", "Pdilt",
    "R3hdm1", "Mal2", "Foxq1", "Ddah1", "Tdo2", "3110082I17Rik", "Lrfn3",
    "Stx1b", "Pop1", "Cbfb", "Tceanc2", "Celsr1", "Malt1", "Crybg3", "C1galt1",
    "Zfp952", "Prok1", "Cela1", "Slc15a5", "Simc1", "Podn", "Sri", "Maml3",
    "Pms1", "Pkd2", "Uckl1", "Atxn1", "Gm32461", "Gm4876", "Tacc1", "Gm10658",
    "Rtn4rl1", "Dazap1", "Dusp11", "Gm20319", "Pde4a", "Dst", "Retsat",
    "1700018L02Rik", "Rprd1a", "Fam126b", "Tmem28", "Col4a2", "Gpr146", "Morc3",
    "Slc37a1", "Gm12602", "Itih2", "Gspt1", "Cyp4f14", "Srsf10", "Nf2", "Tmie",
    "Mbl1", "4921524J17Rik", "C1ra", "Apcs", "Gm26876", "Mug1", "Hip1r",
    "Arhgap15os", "Hc", "Car14", "Fzd7", "Slc39a11", "Zkscan3", "2610020C07Rik",
    "Cmtm6", "Elovl3", "Larp4b", "Ces2b", "Trit1", "Stip1", "Agap1", "Plxna2",
    "Xpnpep3", "Gnai1", "C77080", "Foxa1", "Kdm5b", "Srd5a1", "Zfp697",
    "Echdc2", "Hsd3b3", "Gpr155", "Bard1", "Rida", "Col12a1", "Gstp2", "Abca1",
    "Chrm3", "Serpina10", "Alas2", "Kcp", "Rpa3", "Snx5", "Inhba", "Id2",
    "Dnmt3b", "Slc31a1", "Nek1", "Zmpste24", "Gm5535", "Rtp3", "Cep350",
    "Farp1", "Atg16l2", "Akr1c13", "Nlrp12", "Kctd15", "St3gal3", "Gata4",
    "Cldn1", "C730036E19Rik", "Slc25a40", "Kmt5a", "Odf3b", "Nipsnap1",
    "Zbtb8os", "Plekha1", "Ptch1", "Cux1", "Tmem59", "Qprt", "2810410L24Rik",
    "Cenpu", "Map7d1", "Sox9", "Rab3gap1", "Sh3rf1", "Psmb7", "Pip4k2b", "Ddx6",
    "Nectin3", "Gm15638", "Cyp4a10", "Ranbp9", "Ranbp3l", "Atxn7l1", "Cyp4v3",
    "Fahd1", "Tab2", "Aldh1a1", "Pank3", "Cebpa", "Plbd2", "Dph6",
    "C730027H18Rik", "Ces1f", "Unc13b", "Adh7", "Zc3h14", "Rhot1", "Carmil1",
    "Prdm2", "Olfml1", "Rrn3", "Tcea3", "Tgoln1", "Usp10", "Stat2", "Nutf2",
    "Ptpn3", "Polr2b", "Zfp36l2", "Sft2d3", "Zfp563", "Skap2", "Ctsc", "Bik",
    "Ldha", "Zdhhc14", "Abcb10", "Snhg11", "Mmp15", "Ccdc149", "Hsd17b2",
    "Pdik1l", "Adh4", "S1pr1", "Sort1", "St5", "Ak2", "Nsun2", "Kmo", "Cebpd",
    "Fh1", "Pld1", "Uroc1", "Camk2n1", "Neat1", "Zcchc24", "Chpt1", "Abhd17b",
    "Ngrn", "Mgst1", "H6pd", "Hsd3b7", "Cyp2d10", "Cfap20", "B3gnt9", "Akr1c12",
    "Gtf2ird1", "Ccdc25", "Midn", "Cps1", "Zfp872", "Vps13c", "Dbf4", "C4b",
    "Slc16a1", "Tdrd7", "Rbm34", "Gm19619", "Pim1", "Gtf2i", "Socs7",
    "D1Ertd622e", "Zfp750", "Pparg", "Shtn1", "Cyp4f13", "Ftcd", "Zfyve1",
    "Atp8b1", "Arid5b", "Ak4", "D230025D16Rik", "Ankrd13a", "Mn1", "Nectin2",
    "Ccdc117", "Otulin", "Dand5", "Peli1", "Bhlhe41", "Cys1", "Pde9a", "Fyb2",
    "Smad9", "Nr1h4", "Baz2a", "Skp2", "Tbc1d2b", "Mfhas1", "Dnase2a",
    "Slc45a3", "Eif4ebp3", "2300009A05Rik", "Pycrl", "Rasal2", "Fam135a",
    "Plcxd2", "Tff3", "Mxd4", "Il1rl2", "Aldh7a1", "Pfkfb2", "Scarb2", "Fkbp4",
    "Hsp90ab1", "Tpst2", "Slc25a17", "Saa3"
  )
)


# ==============================================================================
# HELPERS
# ==============================================================================
banner <- function(text) {
  line <- paste(rep("=", 70), collapse = "")
  message("\n", line)
  message(text)
  message(line)
}

paired_wilcox <- function(pos, neg) {
  tryCatch(
    wilcox.test(pos, neg, paired = TRUE, exact = TRUE)$p.value,
    warning = function(w) {
      wilcox.test(pos, neg, paired = TRUE, exact = FALSE)$p.value
    }
  )
}

# Within-animal label permutation
perm_test <- function(meta, feat, sx, n_perm = 10000, min_cells = 10) {
  d <- meta %>%
    mutate(sex = tolower(as.character(sex)),
           status = as.character(Ascl1_status)) %>%
    filter(sex == sx, !is.na(.data[[feat]]), !is.na(status), !is.na(sample))

  ok <- d %>% count(sample, status) %>%
    pivot_wider(names_from = status, values_from = n) %>%
    filter(!is.na(`Ascl1+`), !is.na(`Ascl1-`),
           `Ascl1+` >= min_cells, `Ascl1-` >= min_cells) %>%
    pull(sample)
  d <- d %>% filter(sample %in% ok)
  if (length(ok) < 2) return(NULL)

  obs_diff <- d %>% group_by(sample, status) %>%
    summarise(m = mean(.data[[feat]]), .groups = "drop") %>%
    pivot_wider(names_from = status, values_from = m) %>%
    mutate(diff = `Ascl1+` - `Ascl1-`)
  obs <- mean(obs_diff$diff)

  null <- replicate(n_perm, {
    d %>% group_by(sample) %>%
      mutate(status = sample(status)) %>% ungroup() %>%
      group_by(sample, status) %>%
      summarise(m = mean(.data[[feat]]), .groups = "drop") %>%
      pivot_wider(names_from = status, values_from = m) %>%
      mutate(diff = `Ascl1+` - `Ascl1-`) %>%
      pull(diff) %>% mean()
  })

  list(sex = sx, n_mice = length(ok), obs = obs,
       null_mean = mean(null), null_sd = sd(null),
       z = (obs - mean(null)) / sd(null),
       p_perm = (sum(abs(null) >= abs(obs)) + 1) / (n_perm + 1),
       per_mouse = obs_diff, null = null)
}


# ==============================================================================
# STEP 1: LOAD DATA, DEFINE Ascl1+ AND Ascl1- CELLS
# ==============================================================================
banner("STEP 1: Load h5ad, pre-geriatric hepatocytes, Ascl1 status")

seurat_obj <- schard::h5ad2seurat(H5AD_PATH)

# age may be spelled "pre_geriatric" or "pre-geriatric"
age_std <- gsub("-", "_", as.character(seurat_obj$age))
data <- subset(
  seurat_obj,
  cells = colnames(seurat_obj)[which(seurat_obj$cell_type == "hepatocyte" &
                                     age_std == "pre_geriatric")]
)
rm(seurat_obj); gc()

DefaultAssay(data) <- "RNA"
print(table(data$cell_type, data$age))

data$Ascl1_expr <- FetchData(
  data,
  vars = "Ascl1",
  layer = "data"
)[, 1]

data$Ascl1_status <- ifelse(
  data$Ascl1_expr > 0,
  "Ascl1+",
  "Ascl1-"
)

data$Ascl1_status <- factor(
  data$Ascl1_status,
  levels = c("Ascl1+", "Ascl1-")
)

# Check cell numbers
print(table(data$sex, data$Ascl1_status))


# ==============================================================================
# STEP 2: MODULE SCORE, STATISTICS AND PLOTS PER GENE SET
# ==============================================================================
for (out_tag in names(GENE_SETS)) {

  banner(sprintf("STEP 2: %s", out_tag))
  gene_list_combined <- GENE_SETS[[out_tag]]
  out_file <- function(suffix) {
    file.path(OUTPUT_DIR, sprintf("Ascl1_module_score_%s_%s", out_tag, suffix))
  }

  # -----------------------------------
  # Module score
  # -----------------------------------
  genes_present <- gene_list_combined[gene_list_combined %in% rownames(data)]
  cat(sprintf("Genes found: %d of %d (%s)\n",
              length(genes_present), length(gene_list_combined),
              paste(genes_present, collapse = ", ")))

  data <- AddModuleScore(
    object   = data,
    features = list(genes_present),
    name     = "Combined_ModuleScore"
  )

  feat <- "Combined_ModuleScore1"

  # -----------------------------------
  # Sample-level analysis
  # -----------------------------------
  df_sample <- data@meta.data %>%
    mutate(sample_id = .data[["sample"]]) %>%
    filter(!is.na(sample_id), !is.na(.data[[feat]]),
           !is.na(Ascl1_status), !is.na(sex)) %>%
    mutate(
      status = recode(as.character(Ascl1_status),
                      "Ascl1-" = "neg", "Ascl1+" = "pos"),
      sex = factor(tolower(as.character(sex)), levels = c("female", "male"))
    ) %>%
    filter(status %in% c("neg", "pos")) %>%
    group_by(sample_id, sex, status) %>%
    summarise(mean_score = mean(.data[[feat]], na.rm = TRUE),
              n_cells = n(), .groups = "drop") %>%
    filter(n_cells >= 10) %>%
    pivot_wider(names_from = status, values_from = c(mean_score, n_cells)) %>%
    filter(!is.na(mean_score_neg), !is.na(mean_score_pos)) %>%
    mutate(difference = mean_score_pos - mean_score_neg)

  stats_primary <- df_sample %>%
    group_by(sex) %>%
    summarise(
      n_samples            = n(),
      total_cells_neg      = sum(n_cells_neg),
      total_cells_pos      = sum(n_cells_pos),
      mean_neg             = mean(mean_score_neg),
      mean_pos             = mean(mean_score_pos),
      mean_difference      = mean(difference),
      n_positive_direction = sum(difference > 0),
      cohens_dz = if (n() >= 2 && sd(difference) > 0) {
        mean(difference) / sd(difference)
      } else NA_real_,
      pval = if (n() >= 3) paired_wilcox(mean_score_pos, mean_score_neg) else NA_real_,
      .groups = "drop"
    )

  cat("\nSample-level statistics:\n")
  print(as.data.frame(stats_primary))

  # -----------------------------------
  # Cell-level data (descriptive plots)
  # -----------------------------------
  df_all <- data@meta.data %>%
    filter(!is.na(.data[[feat]])) %>%
    mutate(
      Ascl1_status = factor(as.character(Ascl1_status),
                            levels = c("Ascl1-", "Ascl1+")),
      sex = factor(tolower(as.character(sex)), levels = c("female", "male"))
    ) %>%
    filter(!is.na(sex), !is.na(Ascl1_status))

  cat("\nCells per sex x status:\n")
  print(table(df_all$sex, df_all$Ascl1_status))

  cell_totals <- df_all %>%
    filter(Ascl1_status == "Ascl1+") %>%
    count(sex, name = "n_pos_total") %>%
    complete(sex, fill = list(n_pos_total = 0)) %>%
    mutate(sex = as.character(sex))

  write.csv(stats_primary, out_file("stats.csv"), row.names = FALSE)

  # -----------------------------------
  # Within-animal label permutation
  # -----------------------------------
  set.seed(42)
  res_f <- perm_test(data@meta.data, feat, "female", N_PERM)
  res_m <- perm_test(data@meta.data, feat, "male",   N_PERM)

  for (r in list(res_f, res_m)) {
    if (is.null(r)) { cat("\nInsufficient data for one sex\n"); next }
    cat(sprintf("\n=== %s ===\nmice: %d\nobserved: %+.4f\nnull: %.4f +/- %.4f\nz = %.2f\nP = %.4f\n",
                r$sex, r$n_mice, r$obs, r$null_mean, r$null_sd, r$z, r$p_perm))
    print(as.data.frame(r$per_mouse))
  }

  perm_summary <- bind_rows(lapply(list(res_f, res_m), function(r) {
    if (is.null(r)) return(NULL)
    data.frame(sex = r$sex, n_mice = r$n_mice, observed = r$obs,
               null_mean = r$null_mean, null_sd = r$null_sd,
               z = r$z, p_perm = r$p_perm)
  }))
  write.csv(perm_summary, out_file("permutation.csv"), row.names = FALSE)

  # -----------------------------------
  # Plots (animal-level annotation)
  # -----------------------------------
  ymax <- max(df_all[[feat]], na.rm = TRUE)
  ymin <- min(df_all[[feat]], na.rm = TRUE)
  yrng <- ymax - ymin

  perm_p <- if (!is.null(res_f)) res_f$p_perm else NA_real_
  perm_z <- if (!is.null(res_f)) res_f$z      else NA_real_

  ann <- data.frame(sex = levels(df_all$sex), stringsAsFactors = FALSE) %>%
    left_join(stats_primary %>% mutate(sex = as.character(sex)), by = "sex") %>%
    left_join(cell_totals, by = "sex") %>%
    mutate(
      lab = ifelse(
        is.na(pval),
        sprintf("Ascl1+ cells: %d\n(insufficient)", n_pos_total),
        sprintf("n = %d mice\nperm P = %.4f\ndz = %.2f\n%d/%d same direction",
                n_samples, perm_p, cohens_dz, n_positive_direction, n_samples)
      ),
      col = ifelse(!is.na(pval) & !is.na(perm_p) & perm_p < 0.05, "#d62728", "grey30"),
      sex = factor(sex, levels = levels(df_all$sex)),
      x = 1.5, y = ymax + 0.04 * yrng
    ) %>%
    select(sex, lab, col, x, y)

  fills <- c("Ascl1-" = "#89d9e1", "Ascl1+" = "#fb9a99")

  base <- function(g) {
    g +
      facet_wrap(~ sex, labeller = labeller(sex = tools::toTitleCase)) +
      scale_fill_manual(values = fills) +
      geom_text(data = ann, aes(x = x, y = y, label = lab),
                colour = ann$col, inherit.aes = FALSE, size = 3.2, lineheight = 0.95) +
      coord_cartesian(ylim = c(ymin, ymax + 0.26 * yrng)) +
      labs(y = "Module score", x = NULL) +
      theme_classic() +
      theme(legend.position = "none",
            strip.text = element_text(size = 13, face = "bold"),
            strip.background = element_blank(),
            plot.title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 13),
            axis.text  = element_text(size = 12, colour = "black"))
  }

  g0 <- ggplot(df_all, aes(Ascl1_status, .data[[feat]], fill = Ascl1_status))

  p1 <- base(g0 +
    geom_violin(trim = FALSE, alpha = 0.8, scale = "width", colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", lwd = 0.4)) +
    ggtitle("Violin + box")

  p2 <- base(g0 +
    geom_boxplot(width = 0.55, outlier.size = 0.3, outlier.alpha = 0.25, lwd = 0.5)) +
    ggtitle("Boxplot")

  set.seed(1)
  df_j <- df_all %>%
    group_by(sex, Ascl1_status) %>%
    group_modify(~ dplyr::slice_sample(.x, n = min(800, nrow(.x)))) %>%
    ungroup()

  p3 <- base(g0 +
    geom_boxplot(width = 0.5, outlier.shape = NA, lwd = 0.5, alpha = 0.85) +
    geom_jitter(data = df_j, aes(x = Ascl1_status, y = .data[[feat]]),
                inherit.aes = FALSE, colour = "grey20",
                width = 0.16, size = 0.25, alpha = 0.25)) +
    ggtitle("Box + jitter")

  p4 <- ggplot(df_all, aes(.data[[feat]], fill = Ascl1_status)) +
    geom_density(alpha = 0.6, colour = NA) +
    facet_wrap(~ sex, labeller = labeller(sex = tools::toTitleCase)) +
    scale_fill_manual(values = fills, name = NULL) +
    labs(x = "Module score", y = "Density") +
    theme_classic() +
    theme(strip.text = element_text(size = 13, face = "bold"),
          strip.background = element_blank(),
          plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 13),
          axis.text  = element_text(size = 12, colour = "black"),
          legend.position = "top") +
    ggtitle("Density")

  # ---- paired animal-level ----
  pm <- res_f$per_mouse %>%
    rename(pos = `Ascl1+`, neg = `Ascl1-`) %>%
    mutate(direction = ifelse(diff > 0, "up", "down"))

  pm_long <- pm %>%
    pivot_longer(c(neg, pos), names_to = "status", values_to = "score") %>%
    mutate(status = factor(status, levels = c("neg", "pos")))

  dz_val <- mean(pm$diff) / sd(pm$diff)

  p5 <- ggplot(pm_long, aes(status, score)) +
    geom_line(aes(group = sample, colour = direction), lwd = 0.8, alpha = 0.85) +
    geom_point(aes(fill = status), shape = 21, size = 3.5, stroke = 0.4) +
    scale_fill_manual(values = c(neg = "#89d9e1", pos = "#fb9a99"), guide = "none") +
    scale_colour_manual(values = c(up = "#fb9a99", down = "grey60"), guide = "none") +
    scale_x_discrete(labels = c(neg = "Ascl1-", pos = "Ascl1+")) +
    annotate("text", x = 1.5, y = max(pm_long$score) * 1.05,
             label = sprintf("perm P = %.4f | dz = %.2f\n%d/%d animals same direction",
                             res_f$p_perm, dz_val, sum(pm$diff > 0), nrow(pm)),
             size = 3.3, lineheight = 0.95) +
    coord_cartesian(ylim = c(min(pm_long$score) * 0.98, max(pm_long$score) * 1.14)) +
    labs(y = "Mean module score", x = NULL,
         title = "Animal-level paired (female, n = 4)") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 12))

  # ---- permutation null ----
  p6 <- ggplot(data.frame(null = res_f$null), aes(null)) +
    geom_histogram(bins = 60, fill = "grey85", colour = NA) +
    geom_vline(xintercept = res_f$obs, colour = "#d62728", lwd = 1) +
    annotate("text", x = res_f$obs, y = Inf, vjust = 2,
             hjust = ifelse(res_f$obs < 0, 1.1, -0.1),
             colour = "#d62728", size = 3.4,
             label = sprintf("observed\n%+.4f", res_f$obs)) +
    labs(x = "Permuted within-animal difference", y = "Count",
         title = sprintf("Label permutation (z = %.2f, P = %.4f)",
                         res_f$z, res_f$p_perm)) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 11))

  # ---- save ----
  ggsave(out_file("box_jitter.pdf"), p3, width = 7.5, height = 5.5)

  p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
    plot_annotation(
      title = sprintf("ASCL1 status and module score (%s)", out_tag),
      caption = "Cell-level panels are descriptive only; all statistics are animal-level (within-animal label permutation, n = 4 mice).",
      theme = theme(plot.title = element_text(size = 15, face = "bold"),
                    plot.caption = element_text(size = 9, colour = "grey35"))
    )

  ggsave(out_file("ALLPANELS.pdf"), p_all, width = 11, height = 14)
  message("  [OK] ", out_file("box_jitter.pdf"))
  message("  [OK] ", out_file("ALLPANELS.pdf"))
}

banner("Done")
