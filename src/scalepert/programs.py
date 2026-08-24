PROGRAM_NAMES = ["InnateSensing", "OxidativeMetabolic", "DDRRepair", "DNASensingBridge"]

PROGRAMS = {
    "InnateSensing": [
        "TLR2", "TLR4", "TLR7", "TLR8", "ITGAM", "MYD88", "IRAK1", "IRAK4",
        "TRAF6", "NLRP3", "CCL5", "CCL2", "CXCL10", "CXCL8", "IL1B", "IL6",
        "TNF", "IL18", "RELA", "NFKB1", "IKBKB",
    ],
    "OxidativeMetabolic": [
        "BAX", "BAK1", "BID", "RICTOR", "RHOA", "ROCK1", "ROCK2", "NINJ1",
        "TNFRSF1A", "TNFRSF1B", "HIF1A", "PPARGC1A", "HK2", "LDHA", "PKM",
        "SLC2A1", "SLC2A3", "PDK1", "G6PD", "SOD2", "CAT", "GPX4", "NFE2L2",
        "HMOX1", "NOX4", "CYBB", "DNM1L", "MFN1", "MFN2", "OPA1", "VCAM1",
        "ICAM1", "JUN", "FOS", "MAPK8", "MAPK14",
    ],
    "DDRRepair": [
        "ATM", "ATR", "CHEK1", "CHEK2", "TP53", "H2AFX", "MDC1", "MRE11",
        "RAD50", "NBN", "BRCA1", "BRCA2", "RAD51", "FANCD2", "FANCI", "PALB2",
        "BARD1", "RPA1", "RPA2", "TOPBP1", "PARP1", "PARP2", "XRCC5", "XRCC6",
        "PRKDC", "BLM", "WRN", "EXO1", "CLSPN", "ERCC1", "ERCC2", "ERCC3",
        "ERCC4", "ERCC5", "XPA", "XPC", "DDB1", "DDB2", "LIG1", "LIG3", "LIG4",
        "POLB", "POLD1", "POLE", "FEN1", "UNG", "APEX1", "OGG1", "MUTYH",
        "SMC1A", "SMC3", "RAD17", "RAD9A", "HUS1",
    ],
    "DNASensingBridge": [
        "IFI16", "AIM2", "CGAS", "STING1", "TBK1", "DDX41", "ZBP1", "PYCARD",
    ],
}

GENE_ALIASES = {
    "STING1": "TMEM173",
    "TMEM173": "STING1",
    "H2AFX": "H2AX",
    "CGAS": "MB21D1",
    "RAGE": "AGER",
}

LR_PAIRS = [
    ("CCL5", "CCR5"),
    ("CCL2", "CCR2"),
    ("CXCL10", "CXCR3"),
    ("IL1B", "IL1R1"),
    ("IL18", "IL18R1"),
    ("TNF", "TNFRSF1A"),
    ("IFNG", "IFNGR1"),
    ("IL6", "IL6R"),
    ("HMGB1", "TLR2"),
    ("HMGB1", "TLR4"),
    ("HMGB1", "AGER"),
    ("VCAM1", "ITGA4"),
    ("ICAM1", "ITGAL"),
    ("VEGFA", "FLT1"),
    ("VEGFA", "KDR"),
    ("IL1A", "IL1R1"),
    ("CXCL8", "CXCR1"),
    ("CXCL8", "CXCR2"),
]

DEFAULT_TARGETS = [
    "TLR2", "ITGAM", "CCL5", "HIF1A", "VCAM1", "JUN", "AIM2", "IFI16", "PYCARD",
]

HUB_PROGRAM = {
    "TLR2": "P1",
    "ITGAM": "P1",
    "CCL5": "P1",
    "HIF1A": "P2",
    "VCAM1": "P2",
    "JUN": "P2",
    "AIM2": "P4",
    "IFI16": "P4",
    "PYCARD": "P4",
}


def resolve_gene(gene, var_names):
    names = set(var_names)
    if gene in names:
        return gene
    alias = GENE_ALIASES.get(gene)
    if alias and alias in names:
        return alias
    return None


def resolve_programs(var_names):
    resolved = {}
    for name, genes in PROGRAMS.items():
        mapped = []
        for g in genes:
            r = resolve_gene(g, var_names)
            if r:
                mapped.append(r)
        resolved[name] = sorted(set(mapped))
    return resolved
