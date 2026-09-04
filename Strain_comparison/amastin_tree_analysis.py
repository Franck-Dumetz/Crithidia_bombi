"""
Amastin gene family: sequence extraction, alignment, and bootstrap-supported
NJ tree for CbWHA1 (Crithidia bombi), vs. Leishmania major subfamily-defining
references (alpha/beta/gamma/delta, per Jackson 2010, Mol Biol Evol 27:33-45).

Pipeline is identical in structure to gp63_tree_analysis.py — see that file
for extended comments. Requires: biopython, pyfamsa, matplotlib.
"""

import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from Bio.Seq import Seq
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import Consensus
from pyfamsa import Aligner, Sequence

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
GFF = "CbWHA1_full_annotation.final_numbered.TATE_reclassified.GAPDHcorrected.gff3"
PEP = "transdecoder_collapsed.renumbered.GAPDHcorrected.pep"
GENOME = "CbWHA1_assembly.final_numbered.telopatched2.fasta"

OUT_PNG = "CbWHA1_amastin_tree_bootstrap.png"
OUT_NWK = "CbWHA1_amastin_tree_bootstrap.nwk"

# ---------------------------------------------------------------------------
# 1. Subfamily-labeled L. major references (UniProt)
# ---------------------------------------------------------------------------
REFERENCES = {
    "alpha_LMJF28_1390_Q4Q8C7": "MTMQGQGNNVPQRHQQSQLDSDDESYTGSSGSYSERVSEAPPRQQQSAAVAERRLLPKAE"
    "QGSRPRSLQQTPSPVKRGDVQGTQQEEGGETAQKKGRAARGKAKARTMMNAAWEAPVARV"
    "RAVTGAVSGRDSQRVTIYVFFCIAWVHLLFVILSSALSQIDVVGGGCYTFWGYKANCDTV"
    "SYTRRTVLLQSCRRVRSALQAGSAFSIISVLASAATLVSSWVLCCRLREADRHVRSPSRY"
    "ANMGEVALAQEPDCNDGIENGRQADAPVYDAGNLKKVIMIVVAFSLACELTCWALIAAII"
    "TERYCDEIYLWSTTATYGAGFGLGLTAWFAELIAYIGFVTVV",
    "alpha_LMJF28_1400_Q4Q8C6": "MSQVVRRKKMFDDVSSFSEEDDVDDEAMMQHNVKSPPVMRSAHHVDAAASYLTALSATAA"
    "PQADSDDDSNAFVPKPASSFQQVTKPEGAGRENAGHVTMVARVSVSLTPKDSPRAASPAE"
    "PSAQAQVLVQYPTATMQALQQPESAGAGDLPWSRAVLGSASAYSAPAAVSQRQSWPEEVE"
    "TAVNPEHPGADKASLQKGRRRAYTPCVHSSRANERGGSAAATGGEGAIPAQHPEQMETEA"
    "GAFATTENRGEVKSTRLSPDESDDEMYEQYTPASMNRATAAHIEYSGEEYATHYGAHTGL"
    "RQGSHSAVSKDKDFPFASYIFPRLNPALADYQGARTAVAAVKQGEVMPPNAVETCCAYLM"
    "VADIRVAVYLIVLFVSVALALVSIPTSQLDVVGAACFTYWGYKDNCDNSTYTIPRPLYPT"
    "PYIRRHLGVGAAFSIVTLIVLLVNFAAAVIAVCCLTQAPHTISLNSRIVLGTLGCVGALT"
    "QLISWAVVARIYNSGHYAAGELAYGAGFGLNLSSWVMHLLGVALVFAAPSHFVNRRQQRG",
    "beta_LMJF30_0860_Q4Q7L7": "MGFEALRGRMDVALSMLCSCIVFMFLVTSAPISQFRGRGMNVGGASKLSCVTVWGLKNDC"
    "TANNYDYRPTSIGCARSKQLFQVSEAFYIVAVIVSFLSCLMSGLYFMGIKAKVLLVLLAV"
    "LEVGFALIPWVCMTAVWHGNYCGGSTVKINTSSGKADGVPLGSVLRESFKASAGYGLTVA"
    "AWCTQVIGLVLLIIM",
    "beta_LMJF30_0870_Q4Q7L6": "MANKKSFYNQEYSKHVGAVILFIVSFVAMTFTVCGTPLGMLMIRSWGEDLSGSSAELELN"
    "PCFTLWGLHSDCSKPDYSLRITDSPIVNCSDMHVRFEAAEAFSIVAIFSLVGLFGASWYM"
    "ICGSKIKKAVMLLAVFAIGSTTVPWAIVTAFYYTPFCGLDFLTNTHTRFGAGYALLVTSF"
    "VLQIVGLILFVIFEPNTSKKLEENAKGAASEVWSSTASALR",
    "gamma_LMJF24_1280_Q4QAK4": "MPACLAYYTGAMCFTIIHFLAWAFALVATPTAQFQTPGHGCYTMWGYRKFCGDVPYDLTG"
    "DAAFGCARRTSTMRCGAAFGVMASACGFAGLISAIVLNTQIQIPVIVPFVLAAVCIPCTM"
    "ISWACVASVYNLTMCGDRFGSKYPYTAGFALMVASWGLEIIAVVILACTNWTRPPKEEEN"
    "AADAKH",
    "gamma_LMJF24_1250_Q4QAK7": "MPACLAYYTGAMCFTIIQFFAWAFALVATPTAQFQTPGHGCYTMWGYRKFCGDVPYDLTG"
    "DAAFGCARRTSTMRCGAAFGVMASACGFAGLISAIVLNTQIQIPVIVPFVLAAVCIPCTM"
    "ISWACVASVYNLTMCGDRFGSKYPYTAGFALMVASWGLEIIAVVILACTNWTRPPKEEEN"
    "AADAKH",
    "protodelta_LMJF34_0970_Q4Q3A8": "MSHSFCRVGIAIYCLLQLIAFIFILVGTLIDQFRVQNVDALSNDPCLTIWGFKDKCISLK"
    "WSVRTKDLWKGCPQRLQRFNAAEALSIAAVLISALACLIGFVMLCCCRCLRWLCLILNIL"
    "ATFCGCAVTALMTDAFYNNHEEGLQQYNNSCYALRQNGSVIHPSAIADGNPVATHYKYGA"
    "GFAIYIVGWGLCFINIFFLMLPC",
    "delta_LMJF34_1600_Q4Q341": "MKRSIPLVVYVVVQFVAFLLVLVGTPLDMFRAHNRPGVAQCLTLFGFKLDCESLEYLETV"
    "DTQWVDCPARITRFRLAQAFAIISIFVYGAAFVLGLVLLYGCTIHRWVCLALNIVGAVTL"
    "FIVWAAMVVTYNKDDGQKCLKVRDTGYRLGAGFALLVVAWILDILNIIFLLLPCR",
    "delta_LMJF08_0760_Q4QI93": "MACKLGVAIYVVLQLIAFVAVMVGTGVDMFYNKPEHSSGARVCITLWGLKTDCRKPKITD"
    "SSSVRWALCPIRLKNFRLCQVFAIISILVYGAAFLFGFLLLYCCSGFRWLCLALNIVGAA"
    "TACVVWAVMVVTYRLPEPKCLELSDGYDFGVGFGLFVLAWILDIVDIIFLMLPWQIGEFV"
    "EGEEPSEKEMEKSKNAAQE",
    "delta_LMJF36_4140_Q4Q101": "MLIIRFVLLVLIFVFFVIALVGTVSLPLYSNRITSYNQDGKVEVSLWKIVVGKVTVTGEN"
    "STNVPASTPVRINYAGCEEFRATFRAMQAFAIGGTVFGFYALLVSCLQCFCRLKVKLPLF"
    "LFLFLAMLCELCVVFIGGAAYSKEFCKNLEKNGNLTTIIFKGAGYKLDTAFILQVVALVG"
    "YAICTIITPFTQQLWCGKC",
}
FRIENDLY_NAME = {
    "alpha_LMJF28_1390_Q4Q8C7": "α-amastin L. major (LMJF28.1390)",
    "alpha_LMJF28_1400_Q4Q8C6": "α-amastin L. major (LMJF28.1400)",
    "beta_LMJF30_0860_Q4Q7L7": "β-amastin L. major (LMJF30.0860)",
    "beta_LMJF30_0870_Q4Q7L6": "β-amastin L. major (LMJF30.0870)",
    "gamma_LMJF24_1280_Q4QAK4": "γ-amastin L. major (LMJF24.1280)",
    "gamma_LMJF24_1250_Q4QAK7": "γ-amastin L. major (LMJF24.1250)",
    "protodelta_LMJF34_0970_Q4Q3A8": "proto-δ-amastin L. major (LMJF34.0970)",
    "delta_LMJF34_1600_Q4Q341": "δ-amastin L. major (LMJF34.1600)",
    "delta_LMJF08_0760_Q4QI93": "δ-amastin L. major (LMJF08.0760)",
    "delta_LMJF36_4140_Q4Q101": "δ-amastin L. major (LMJF36.4140)",
}

# ---------------------------------------------------------------------------
# 2. Pull CbWHA1 amastin sequences  (same helpers as the GP63 script)
# ---------------------------------------------------------------------------
def load_fasta_seq(path, target):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name == target:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                if name == target:
                    buf.append(line.strip())
        if name == target:
            seqs[name] = "".join(buf)
    return seqs.get(target, "")


def load_pep(path):
    seqs, header, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if header:
            seqs[header] = "".join(buf)
    return seqs


def get_amastin_genes(gff_path):
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if "\tCDS\t" not in line or "amastin" not in line.lower():
                continue
            fields = line.rstrip("\n").split("\t")
            parent = re.search(r"Parent=([^;]+)", fields[8]).group(1)
            base = re.match(r"(CbWHA1\.\d+\.\d+)", parent).group(1)
            genes.setdefault(base, []).append(
                (parent, fields[0], int(fields[3]), int(fields[4]), fields[6])
            )
    return genes


def build_cb_sequences(gff_path, pep_path, genome_path):
    pep_seqs = load_pep(pep_path)
    genes = get_amastin_genes(gff_path)
    genome_cache = {}

    def get_chr(chrom):
        if chrom not in genome_cache:
            genome_cache[chrom] = load_fasta_seq(genome_path, chrom)
        return genome_cache[chrom]

    result = {}
    for base, entries in genes.items():
        protein_id, chrom, start, end, strand = entries[0]
        if protein_id in pep_seqs:
            result[base] = pep_seqs[protein_id].rstrip("*")
        else:
            sub = get_chr(chrom)[start - 1 : end]
            if strand == "-":
                sub = str(Seq(sub).reverse_complement())
            result[base] = str(Seq(sub).translate()).rstrip("*")
    return result


# ---------------------------------------------------------------------------
# 3. Align, 4. NJ tree + bootstrap  (identical helpers to gp63 script)
# ---------------------------------------------------------------------------
def align_sequences(seq_dict):
    records = [Sequence(name.encode(), seq.encode()) for name, seq in seq_dict.items()]
    aligner = Aligner(guide_tree="upgma")
    msa = aligner.align(records)
    return {s.id.decode(): s.sequence.decode() for s in msa}


def write_fasta(seq_dict, path):
    with open(path, "w") as out:
        for name, seq in seq_dict.items():
            out.write(f">{name}\n{seq}\n")


def build_bootstrap_tree(aligned_fasta_path, n_replicates=500):
    aln = AlignIO.read(aligned_fasta_path, "fasta")
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, "nj")
    target_tree = constructor.build_tree(aln)
    replicate_trees = list(Consensus.bootstrap_trees(aln, n_replicates, constructor))
    support_tree = Consensus.get_support(target_tree, replicate_trees)
    support_tree.ladderize()
    return support_tree


# ---------------------------------------------------------------------------
# 5. Plotting — color by nearest-neighbor subfamily / genomic cluster
# ---------------------------------------------------------------------------
def color_for(name):
    if name is None:
        return "black"
    if name.startswith("CbWHA1.29."):
        return "#1b9e77"  # delta/proto-delta — chr29 tandem array
    if name.startswith("CbWHA1.21."):
        return "#d95f02"  # gamma — chr21
    if name.startswith("CbWHA1.30."):
        return "#7570b3"  # beta — chr30
    if name.startswith("CbWHA1.32.") or name.startswith("CbWHA1.12.") or name.startswith("CbWHA1.14.001300"):
        return "#e7298a"  # alpha-like grouping (chr32 pair + weakly related singletons)
    if name.startswith("CbWHA1."):
        return "#666666"  # unresolved / divergent singleton
    if name.startswith("alpha_"):
        return "#e7298a"
    if name.startswith("beta_"):
        return "#7570b3"
    if name.startswith("gamma_"):
        return "#d95f02"
    if name.startswith(("protodelta_", "delta_")):
        return "#1b9e77"
    return "black"


def plot_tree(support_tree, out_png, title):
    label_color = {}

    def label_func(clade):
        if not clade.name or clade.name.startswith("Inner"):
            return ""
        lbl = FRIENDLY_NAME.get(clade.name, clade.name)
        label_color[lbl] = color_for(clade.name)
        return lbl

    for clade in support_tree.find_clades():
        if clade.name and not clade.name.startswith("Inner"):
            clade.color = color_for(clade.name)

    def branch_label(clade):
        if not clade.is_terminal() and clade.confidence is not None:
            return f"{clade.confidence:.0f}"
        return None

    fig = plt.figure(figsize=(15, 17))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(
        support_tree,
        axes=ax,
        do_show=False,
        label_func=label_func,
        label_colors=lambda lbl: label_color.get(lbl, "black"),
        branch_labels=branch_label,
        show_confidence=False,
    )

    for t in ax.texts:
        txt = t.get_text()
        if not txt:
            continue
        if any(ch.isalpha() for ch in txt):
            t.set_fontweight("bold")
            t.set_fontsize(9)
        else:
            t.set_fontsize(7.5)
            t.set_color("#c00000")
            t.set_fontweight("bold")
            offset = mtransforms.offset_copy(t.get_transform(), fig=fig, x=0, y=4, units="points")
            t.set_transform(offset)

    ax.set_title(title, fontsize=10)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    cb_seqs = build_cb_sequences(GFF, PEP, GENOME)
    print(f"CbWHA1 amastin sequences: {len(cb_seqs)}")

    all_seqs = {**cb_seqs, **REFERENCES}
    write_fasta(all_seqs, "combined.fa")

    aligned = align_sequences(all_seqs)
    write_fasta(aligned, "aligned.fa")

    tree = build_bootstrap_tree("aligned.fa", n_replicates=500)
    Phylo.write(tree, OUT_NWK, "newick")
    # See note in gp63_tree_analysis.py: plot from the in-memory tree, not a
    # reloaded .nwk, because Bio.Phylo's writer merges confidence into node
    # names on write.

    title = (
        "C. bombi amastin family vs. L. major subfamily-labeled references\n"
        "NJ tree, p-distance, 500 bootstrap replicates (red = support at internal nodes)\n"
        "green=δ/proto-δ-like  purple=β-like  orange=γ-like  pink=α-like"
    )
    plot_tree(tree, OUT_PNG, title)
    print("done")
