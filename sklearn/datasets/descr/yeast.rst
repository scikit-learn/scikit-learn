Yeast Dataset
=============

Notes
-----
Data Set Characteristics:
    :Number of Instances: 1484

    :Number of Attributes: 1 string attribute, 8 numeric attributes and the class

    :Attribute Information:
        - Sequence Name: Accession number for the SWISS-PROT database
        - mcg: McGeoch's method for signal sequence recognition.
        - gvh: von Heijne's method for signal sequence recognition.
        - alm: Score of the ALOM membrane spanning region prediction program.
        - mit: Score of discriminant analysis of the amino acid content of the N-terminal
        region (20 residues long) of mitochondrial and non-mitochondrial proteins.
        - erl: Presence of "HDEL" substring (thought to act as a signal for retention in
        the endoplasmic reticulum lumen). Binary attribute.
        - pox: Peroxisomal targeting signal in the C-terminus.
        - vac: Score of discriminant analysis of the amino acid content of vacuolar and
        extracellular proteins.
        - nuc: Score of discriminant analysis of nuclear localization signals of nuclear
        and non-nuclear proteins.

    :Summary Statistics: #TODO

    :Missing Attribute Values: None

    :Class Distribution: 463 - CYT, 5 - ERL, 35 - EXC, 44 - ME1, 51 - ME2,
                         163 - ME3, 244 - MIT, 429 - NUC, 20 - POX, 30 - VAC

    :Creator: Kenta Nakai

    :Donor: Paul Horton

    :Date: September, 1996

Predicted Attribute: Localization site of protein. ( non-numeric ).

References
----------
Paul Horton & Kenta Nakai, "A Probablistic Classification System for Predicting the Cellular Localization Sites of Proteins", Intelligent Systems in Molecular Biology, 109-115. St. Louis, USA 1996.

The references below describe a predecessor to this dataset and its development. They also give results (not cross-validated) for classification by a rule-based expert system with that version of the dataset:

Kenta Nakai & Minoru Kanehisa, "Expert Sytem for Predicting Protein Localization Sites in Gram-Negative Bacteria", PROTEINS: Structure, Function, and Genetics 11:95-110, 1991.

Kenta Nakai & Minoru Kanehisa, "A Knowledge Base for Predicting Protein Localization Sites in Eukaryotic Cells", Genomics 14:897-911, 1992.
