


product_path <- "product/02.1_summarize_KeBOB_data"
if (!dir.exists(product_path)) {
  cat("Creating ", product_path, "\n", sep = "")
  dir.create(product_path)
}

date_code <- "20250817"


load("intermediate/01.1_gather_KeBOB/KeBOB_references.Rdata")
load("intermediate/01.3_gather_chembl/chembl_KeBOB_documents.Rdata")


plot_data <- KeBOB_references |>
    dplyr::mutate(
        in_chembl = DOI %in% chembl_KeBOB_documents$doi,
        class = dplyr::case_when(
            in_chembl ~ "KeBOB+ChEMBL",
            TRUE ~ "KeBOB Only")) |>
    dplyr::select(
        class,
        DOI,
        year = Year,
        journal_abbrev) |>
    dplyr::bind_rows(
        chembl_KeBOB_documents |>
        dplyr::transmute(
            class = "ChEMBL Only",
            DOI = doi,
            year,
            journal_abbrev = journal))

plot_data <- plot_data |>
    dplyr::count(class, year) |>
    dplyr::group_by(class) |>
    dplyr::arrange(year) |>
    dplyr::mutate(
        cumulative_n = cumsum(n)) |>
    dplyr::ungroup()
        



plot <- ggplot2::ggplot(data = plot_data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = year,
            y = sqrt(cumulative_n),
            color = class,
            group = class),
        linewidth = 1.5) +
    ggplot2::scale_y_continuous(
        "Cumulative N References",
        breaks = sqrt(c(0, 20, 100, 300, 600, 1000, 1600)),
        labels =      c(0, 20, 100, 300, 600, 1000, 1600)) +
    ggplot2::scale_x_continuous(
        "Year",
        breaks = c(1960, 1980, 2000, 2015, 2025),
        labels = c(1960, 1980, 2000, 2015, 2025)) +
    ggplot2::scale_color_discrete("") +
    ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8))

ggplot2::ggsave(
    filename = paste0(product_path, "/cumumative_n_references_", date_code, ".pdf"),
    width = 4.5,
    height = 3.5,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = paste0(product_path, "/cumumative_n_references_", date_code, ".png"),
    width = 4.5,
    height = 3.5)


plot_data |>
    readr::write_tsv(
        file = paste0(product_path, "/cumumative_n_references_", date_code, ".tsv"))


chembl_KeBOB_documents |>
    dplyr::count(doc_type, journal) |>
    dplyr::arrange(dplyr::desc(doc_type), dplyr::desc(n)) |>
    data.frame()

#       doc_type                     journal   n
# 1  PUBLICATION                  J Med Chem 789
# 2  PUBLICATION        Bioorg Med Chem Lett 401
# 3  PUBLICATION             Bioorg Med Chem 128
# 4  PUBLICATION              Eur J Med Chem  76
# 5  PUBLICATION           ACS Med Chem Lett  51
# 6  PUBLICATION                  J Nat Prod  36
# 7  PUBLICATION                 Medchemcomm  12
# 8  PUBLICATION                RSC Med Chem   8
# 9  PUBLICATION Antimicrob Agents Chemother   2
# 10 PUBLICATION                 J Biol Chem   2
# 11 PUBLICATION                Med Chem Res   2
# 12 PUBLICATION    Proc Natl Acad Sci U S A   2
# 13 PUBLICATION               Nat Chem Biol   1
# 14 PUBLICATION                  Nat Commun   1
# 15 PUBLICATION                     Science   1
# 16      PATENT                        <NA>  36
# 17     DATASET                        <NA>  56

KeBOB_references |>
    dplyr::count(journal_abbrev) |>
    dplyr::arrange(dplyr::desc(n)) |>
    data.frame()

                                                           journal_abbrev  n
1                                                                    <NA> 47
2                                                           J. Med. Chem. 24
3                                                               Molecules 16
4                                                     ACS Chem. Neurosci. 12
5                                                  British J Pharmacology 10
6                                           Proc. Natl. Acad. Sci. U.S.A. 10
7                                                       Neuropharmacology  9
8                                                              Nat Commun  8
9                                                                  Nature  8
10                                                               Peptides  8
11              The Journal of Pharmacology and Experimental Therapeutics  7
12                                       European Journal of Pharmacology  6
13                                                                   IJMS  6
14                                                   J. Chem. Inf. Model.  6
15                                                 Molecular Pharmacology  6
16                                                                   Cell  5
17                                European Journal of Medicinal Chemistry  5
18                                                      Front. Pharmacol.  5
19                                                          Nat Chem Biol  5
20                                                     Psychopharmacology  5
21                                               Biochemical Pharmacology  4
22                                                            ChemMedChem  4
23                                        Journal of Biological Chemistry  4
24                                                                Sci Rep  4
25                                     Trends in Pharmacological Sciences  4
26                                                           Biochemistry  3
27                                   Bioorganic &amp; Medicinal Chemistry  3
28                                                    Biophysical Journal  3
29                                                                 Neuron  3
30                                                   Neuropsychopharmacol  3
31                                                           Sci. Signal.  3
32                                                                Science  3
33                                                         ACS Cent. Sci.  2
34                                                   ACS Med. Chem. Lett.  2
35  American Journal of Physiology-Lung Cellular and Molecular Physiology  2
36                                                      Angew Chem Int Ed  2
37                                                           Arch Toxicol  2
38                                                           Biomolecules  2
39                                              Drug Testing and Analysis  2
40                                            Drug and Alcohol Dependence  2
41                                                       Front. Neurosci.  2
42                                                       J. Phys. Chem. B  2
43                                             Medicinal Research Reviews  2
44                                                    Nat Struct Mol Biol  2
45                                                               PLoS ONE  2
46                                                         RSC Med. Chem.  2
47                                                              Structure  2
48                                                              ACS Omega  1
49                                                 Addiction Neuroscience  1
50                                                       Advanced Science  1
51                                                    Annu. Rev. Biochem.  1
52                                         Annu. Rev. Pharmacol. Toxicol.  1
53                                                   Archiv der Pharmazie  1
54                                                                   BCCS  1
55                                                               BMC Biol  1
56                                             Behavioural Brain Research  1
57                    Biochemical and Biophysical Research Communications  1
58                     Biochimica et Biophysica Acta (BBA) - Biomembranes  1
59                                                           Biomedicines  1
60                           Bioorganic &amp; Medicinal Chemistry Letters  1
61                                                   Bioorganic Chemistry  1
62                                                            Biopolymers  1
63                                                         Brain Research  1
64                                                Brain Research Bulletin  1
65                                         British Journal of Anaesthesia  1
66                                                                    CMC  1
67                                                    Cellular Signalling  1
68                                                             Chem. Rev.  1
69                                                            ChemBioChem  1
70                                                          ChemPhotoChem  1
71                                                           ChemPlusChem  1
72                                                 Chemistry A European J  1
73                                        Chemistry and Physics of Lipids  1
74                                                              Chirality  1
75                                                    Clinical Toxicology  1
76                                             Clinical Translational Sci  1
77                                                            Commun Biol  1
78                                        Current Opinion in Cell Biology  1
79                                                                   DDDT  1
80                                                               Database  1
81                                                   Drug Discovery Today  1
82                                             Eur J Nucl Med Mol Imaging  1
83                                               European J Oral Sciences  1
84                                                          Exp Brain Res  1
85                                                            Experientia  1
86                                                           FEBS Letters  1
87                                                           Front. Chem.  1
88                                                 Front. Neural Circuits  1
89                                                      Future Med. Chem.  1
90                                                 Helvetica Chimica Acta  1
91                                                             Il Farmaco  1
92                                                    Int J Pept Res Ther  1
93                                                            J Mol Model  1
94                                                    J Neuroinflammation  1
95                                                             J Nucl Med  1
96                                                           J. Neurosci.  1
97                                                   J. Phys. Chem. Lett.  1
98                                                       J. Toxicol. Sci.  1
99                                                                    JPR  1
100                                      Journal of Chemical Neuroanatomy  1
101                                           Journal of Chromatography A  1
102                  Journal of Enzyme Inhibition and Medicinal Chemistry  1
103                                          Journal of Ethnopharmacology  1
104                                          Journal of Molecular Biology  1
105                                             Journal of Neurochemistry  1
106                                     Journal of Translational Research  1
107                                                                  Life  1
108                                                         Life Sciences  1
109                                                                  MBoC  1
110                                                               Methods  1
111                                                        Molecular Cell  1
112                                                          Nat Neurosci  1
113                                           Neurogastroenterology Motil  1
114                                                  Neuroscience Letters  1
115                                                    Org. Biomol. Chem.  1
116                                                      PLoS Comput Biol  1
117                                                             Pain,PAIN  1
118                                                       Pharmaceuticals  1
119                                                         Pharmaceutics  1
120                                              Pharmacological Research  1
121                                               Pharmacological Reviews  1
122                                Pharmacology Biochemistry and Behavior  1
123                                               Phys. Chem. Chem. Phys.  1
124                                                              Proteins  1
125                                              Psychiatry International  1
126                                             Sig Transduct Target Ther  1
127                                                     The FASEB Journal  1
128                                                      The FEBS Journal  1
129                                        The Journal of Clinical Pharma  1
130                                                    Toxicology Letters  1
131                                               Trends in Neurosciences  1
132                                                              iScience  1
