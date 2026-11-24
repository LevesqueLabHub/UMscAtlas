# ============================================================================
# Victoria's UM scAtlas colour palettes - VERBATIM copy of `colors.R` shipped
# with the October 2026 reprocessed dataset (Switch Drive "UM scAtlas", 2026-10-01).
# Supersedes config_colors_2026.R (July). Do NOT hand-edit: re-copy from Switch
# Drive when Victoria updates it. Aliases needed by the GloScope report are added
# in the Rmd, not here.
# ============================================================================

# location colors: defines colors for different anatomical locations in the dataset (healthy iris, healthy choroid, choroid, ciliary body, iris, ciliochoroidal, brain, liver, thyroid, subcutaneous)
location_colors <- c("non-malignant iris" = "#F3BC6AFF", "non-malignant choroid" = "#BB5137FF", "non-malignant ciliary body"  = "#C7782BFF", "iris" = "#00778BFF", "choroid" = "#211747FF", "ciliary body" = "#B6E2E0FF", "ciliochoroidal" = "#6CACE4FF","brain" = "#661246FF", "thyroid" = "#AE1357FF", "liver" = "#E2ADDDFF", "subcut" = "#F06292FF")


# primary_location_colors colors defines colors for intra ocular loaction of primary tumors (pUM)
primary_location_colors <- c("iris" = "#00778BFF", "choroid" = "#211747FF", "ciliary body" = "#B6E2E0FF", "ciliochoroidal" = "#6CACE4FF")


# met_location_colors colors defines colors for anatomical loaction of metastatic tumors (mUM)
met_location_colors <- c("brain" = "#661246FF", "thyroid" = "#AE1357FF", "liver" = "#E2ADDDFF", "subcut" = "#F06292FF")


# healthy_location_colors defines colors for healthy iris and healthy choroid samples
healthy_location_colors <- c("non-malignant iris" = "#F3BC6AFF", "non-malignant choroid" = "#BB5137FF", "non-malignant ciliary body"  = "#C7782BFF")


# origin_colors defines disease status (healthy tissue, primary tumor, metastatic tumor)
origin_colors <- c("primary" = "royalblue3", "non-malignant" = "sienna3", "metastasis" = "violetred4")


# organ_colors defines colors of organ from which sample was derived (eye, brain, thyroid, liver, skin)
organ_colors <- c("eye" = "#5498BAFF", "brain" = "#661246FF", "thyroid" = "#AE1357FF", "liver" = "#E2ADDDFF", "skin" = "#F06292FF")


# second_mutations_col defines colors for 3 secondary mutations defined in tumors (BAP1, SF3B1, EIF1AX)
second_mutations_col <- c("BAP1" = "#1F867BFF", "SF3B1" = "#EAC541FF", "EIF1AX" = "#7B1C5DFF", "unknown" = "#E1E1E1FF", "none" = "#343434FF")


# primary_mutations_col defines colors for 4 primary mutations defined in tumors (GNAQ, GNA11, PLCB4, CYSLTR2)
primary_mutations_col <- c("GNAQ" = "#70395CFF", "GNA11" = "#B565A7FF", "PLCB4" = "#D5B77DFF", "CYSLTR2" = "#A89E5EFF", "unknown" = "#E1E1E1FF", "none" = "#343434FF")
# old GNA11 color: "GNA11" = "#846D86FF"

# gender_col colors to discriminate between M & F from metadata
gender_col <- c("M" = "#58827EFF", "F" = "#611554FF", "unknown" = "#E1E1E1FF")

# ageGroup_col colors to discriminate between Young, Middle and Old age groups for GloScope Analysis
ageGroup_col <- c("Young" = "#9386A6FF", "Middle" = "#BFCDD9FF", "Old" = "#7F8C72FF")

# majority_celltype_col defines colors for major cell types (stored in majority_celltype)
majority_celltype_col <- c("Glial cells" = "#5A5156", "B/Plasma cells" = "lightslateblue", "Myeloid cells" = "lightblue", "T/NK cells" = "#325A9B", "Melanoma cells" = "#B00068", "Melanocytes" = "#C075A6", "Endothelial cells" = "#006400", "Epithelial cells" = "#00A878", "Fibroblasts" = "#ADFF2F", "Pericytes" = "#6B8E23", "Smooth muscle cells" = "#1CFFCE", "Photoreceptor cells" = "#EFEFCFFF", "Retinal Pigmented Epithelial cells" = "#D5B77DFF", "Proliferating cells" = "#ADD0B5FF", "Granulocytes" = "#184860FF")


# tcell_col defines colors for t cell subtyping (stored in Subtyping)
tcell_col <- c("Proliferating T cells" = "#001889FF", "CD4_N" = "#BED7FBFF", "NK_REST" = "#004B84FF", "CD8_EM" = "#C6DBEFFF", "CD4_EM" = "#00366CFF", "CD8_N" = "#A4C7F3FF", "CD8_RM" = "#00619FFF", "gdT" = "#88B4E8FF","CD8_EMRA" = "#0077BCFF", "CD8_EX" = "#90A5C9FF", "CD4_EX" = "#408DCEFF", "MAIT" = "#486078FF", "CD4_REG" = "#68A1DCFF", "NK_CYTO" = "#7797CAFF")


# otherImmune_color defines colors for other immune cells mostly identified by Azimuth (stored in Subtyping)
otherImmune_color <- c("other Immune cells" = "#6B0077FF", "MAIT" = "#486078FF", "cDC1" = "#6C207EFF", "cDC2" = "#714C94FF", "pDC" = "#6E3587FF", "Plasma cells" = "#786AC6FF", "B cells" = "#3B2D7BFF", "Platelets" = "#091626FF")


# macro_col defines colors for macrophage subtyping (stored in Subtyping)
macro_col <- c("activated tissue-infiltrating Macrophages" = "#023FA5FF", "Intermediate TAMs" = "#485CA8FF", "ISG Macrophages" = "#6A76B2FF", "M1-like mono-derived TAMs" = "#878FBDFF", "M1-like TAMs" = "#A1A6C8FF", "M2-like TAMs" = "#B7BBD1FF", "Mature TAMs" = "#CBCDD9FF", "Mono-derived Macrophages" = "#DADBDFFF", "proliferating TAMs" = "#E2E2E2FF")


# TME_SubTyp_col combining tcell_col, otherImmune_col, macro_col and other cells from TME that were not subtyped
TME_SubTyp_col <- c(tcell_col, otherImmune_color, macro_col, "Glial cells" = "#5A5156", "Endothelial cells" = "#006400", "Epithelial cells" = "#00A878", "Fibroblasts" = "#ADFF2F", "Pericytes" = "#6B8E23", "Smooth muscle cells" = "#1CFFCE", "Photoreceptor cells" = "#EFEFCFFF", "Retinal Pigmented Epithelial cells" = "#D5B77DFF", "Proliferating T cells" = "#001889FF", "Granulocytes" = "#184860FF", "Myeloid cells" = "lightblue") 


# fullData_SubTyp_col combines TME_SubTyp_col with Melanoma and Melanocyte colors to show subtyping of TME in also full Dataset
fullData_SubTyp_col <- c(TME_SubTyp_col, "Melanoma cells" = "#B00068", "Melanocytes" = "#C075A6")


# iCNV_colors defines colors to color UMAP by iCNV status stored in Malignant column (levels: Helathy, Malignant, Reference, NA -> Other cells)
iCNV_colors <- c("non-Malignant" = "#ED8B58FF", "Malignant" = "#D84860FF", "Reference" = "#070707FF", "Other cells" = "#666879FF")

healthy_col <- c("hCB54" = "#FEA433", "hCB55" = "#662506", "hCB56" = "#FFFACF", "hCh18" = "#812C04", "hCh20" = "#FEF5B9", "hCh21" = "#9C3503", "hCh22" = "#FEEBA2", "hCh23" = "#B74202", "hCh24" = "#FEDE88", "hCh25" = "#D05004", "hCh28" = "#FECE65", "hCh29" = "#E1640E", "hIr19" = "#FEBB47", "hIr51" = "#EF7818", "hIr52" = "#FFFFE5","hIr53" = "#F98E23")

primary_col <- c("pUM01" = "#97D586", "pUM02" = "#1E1D3B", "pUM03" = "#006D71", "pUM04" = "#2A5783", "pUM05" = "#5FA598", "pUM06" = "#3E6AB0", "pUM09" = "#61A8C7", "pUM10" = "#A4CEE8", "pUM11" = "#BFE8EE", "pUM12" = "#24693D", "pUM14" = "#0A9396", "pUM15" = "#49964F", "pUM16" = "#055E8F", "pUM17" = "#003F5C", "pUM19" = "#8BC1AF", "pUM20" = "#C9E9D5", "pUM21" = "#40729C", "pUM22" = "#69B15E", "pUM23" = "#77BBCF", "pUM24" = "#2D8884", "pUM25" = "#B7DDC9", "pUM26" = "#328445", "pUM27" = "#A5D2BE", "pUM29" = "#D8F3F4", "pUM32" = "#09622A")


mets_col <- c("mUM01" = "#2D204C", "mUM02" = "#B86897", "mUM04" = "#3F305B", "mUM05" = "#DFD3E4", "mUM06" = "#D2A3C1", "mUM07" = "#D8C3D8", "mUM08" = "#66507A", "mUM09" = "#D4B2CC", "mUM11" = "#7D5B87", "mUM13" = "#D385AC", "mUM14" = "#926390", "mUM15" = "#D394B8", "mUM16" = "#523F6C", "mUM18" = "#A56693", "mUM19" = "#E5E5F0", "mUM20" = "#CC74A0")

sample_col <- c(healthy_col, primary_col, mets_col)

patient_cols <- c("P1" = "#97D586", "P2" = "#1E1D3B", "P3" = "#006D71", "P4" = "#2A5783", "P5" = "#5FA598", "P6" = "#3E6AB0", "P7" = "#61A8C7", "P8" = "#A4CEE8", "P9" = "#BFE8EE", "P10" = "#24693D", "P11" = "#0A9396", "P12" = "#49964F", "P13" = "#055E8F", "P14" = "#003F5C", "P15" = "#812C04", "P16" = "#8BC1AF", "P17" = "#C9E9D5", "P18" = "#40729C", "P19" = "#69B15E", "P20" = "#77BBCF", "P21" = "#2D8884", "P22" = "#B7DDC9", "P23" = "#2D204C", "P24" = "#B86897", "P25" = "#3F305B", "P26" = "#DFD3E4", "P27" = "#D2A3C1", "P28" = "#D8C3D8", "P29" = "#66507A", "P30" = "#D4B2CC", "P31" = "#E5E5F0", "P32" = "#926390", "P33" = "#D394B8", "P34" = "#523F6C", "P35" = "#A56693", "P36" = "#CC74A0", "P37" = "#EF7818", "P38" = "#FFFFE5", "P39" = "#F98E23", "P40" = "#FEA433", "P41" = "#FFFACF", "P42" = "#328445", "P43" = "#A5D2BE", "P44" = "#FECE65", "P45" = "#D8F3F4", "P46" = "#09622A")

batch_col <- c("Batch01"= "#4F3855FF", "Batch02" = "#846D86FF", "Batch03" = "#EFEFCFFF", "Batch04" = "#D5B77DFF", "Batch05" = "#A89E5EFF")
