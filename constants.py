SCAFFOLD_LENGTH = 300
C1P_BASEDIST = 10.7
# P_BASEDIST = 20.0 ?
MINOR_DIST = 11.7
MAJOR_DIST = 17.2
BP_TWIST = 34.6
BP_DIST = 3.4
P_DIST = 6.6
C1P_DIST = 4.9

WC_DICT = {"DC": "DG", "DG": "DC", "DT": "DA", "DA": "DT",
           "C": "G", "G": "C", "T": "A", "A": "T",
           "CYT": "GUA", "GUA": "CYT", "THY": "ADE", "ADE": "THY"}
WC_HBONDS = {"DA": ("N6", "N1"), "A": ("N6", "N1"), "ADE": ("N6", "N1"),
             "DT": ("O4", "N3"), "T": ("O4", "N3"), "THY": ("O4", "N3"),
             "DG": ("O6", "N1", "N2"), "G": ("O6", "N1", "N2"), "GUA": ("O6", "N1", "N2"),
             "DC": ("N4", "N3", "O2"), "C": ("N4", "N3", "O2"), "CYT": ("N4", "N3", "O2")}
(N6O4, N1N3_AT) = (2.95, 2.88)
(O6N4, N1N3_GC, N2O2) = (2.85, 2.85, 2.85)
WC_HBONDS_DIST = {"DA": (N6O4, N1N3_AT), "A": (N6O4, N1N3_AT), "ADE": (N6O4, N1N3_AT),
                  "DT": (N6O4, N1N3_AT), "T": (N6O4, N1N3_AT), "THY": (N6O4, N1N3_AT),
                  "DG": (O6N4, N1N3_GC, N2O2), "G": (O6N4, N1N3_GC, N2O2), "GUA": (O6N4, N1N3_GC, N2O2),
                  "DC": (O6N4, N1N3_GC, N2O2), "C": (O6N4, N1N3_GC, N2O2), "CYT": (O6N4, N1N3_GC, N2O2)}
