"""
GP63/leishmanolysin gene family: sequence extraction, alignment, and
bootstrap-supported NJ tree for CbWHA1 (Crithidia bombi).

Pipeline:
  1. Pull all genes annotated as GP63/leishmanolysin from the GFF3
     (gene_product tag), get their protein sequence from the .pep file
     (falling back to direct genome extraction + translation for genes
     TransDecoder didn't emit a clean ORF for).
  2. Combine with two reference sequences (C. fasciculata, L. major).
  3. Align with FAMSA.
  4. Build a neighbor-joining tree (p-distance) and assess node support
     with 500 bootstrap replicates (Bio.Phylo.Consensus).
  5. Plot with tips colored by genomic cluster and bootstrap values
     shown at internal nodes.

Requires: biopython, pyfamsa, matplotlib
    pip install biopython pyfamsa matplotlib --break-system-packages
"""

import re
import json
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
# 0. Paths — edit these for your own setup
# ---------------------------------------------------------------------------
GFF = "CbWHA1_full_annotation.final_numbered.TATE_reclassified.GAPDHcorrected.gff3"
PEP = "transdecoder_collapsed.renumbered.GAPDHcorrected.pep"
GENOME = "CbWHA1_assembly.final_numbered.telopatched2.fasta"

# Genes to exclude (e.g. partial/divergent alignments, likely gene-model
# artifacts — see accompanying notes)
EXCLUDE = {"CbWHA1.31.003560"}

OUT_PNG = "CbWHA1_GP63_tree_bootstrap.png"
OUT_NWK = "CbWHA1_GP63_tree_bootstrap.nwk"

# ---------------------------------------------------------------------------
# 1. Reference sequences (UniProt)
# ---------------------------------------------------------------------------
REFERENCES = {
    "GP63_CRIFA_Q06031": (  # Crithidia fasciculata, CFAC1_040011900
        "MHAPPTATRRSGPRRTHGIMARLVRLAAGVLVVTLVIGALTALSADDAKTHPHKVCIHDE"
        "LQQSLLDSVAQQGLAPQRVSRVGLPYVASATAAPAAQVGGVDFALAGDSAPDVTRSAEWG"
        "ELRITVSAEELTDPAYHCATVGQVISNHIDDYVTCTADDIMTAEKLDILMNYLIPEALQM"
        "HKDRLQVQQVQGTWKVARMTSYCGRFKVPEEHFTTGLSNTDFVLYVASVPTSPGVLAWAN"
        "TCQVFSNDQPAVGVINIPAATITERYDHLMVHAVTHEIAHSLGFSNAFFTNTGIGQFVTG"
        "VRGNPDTVPVINSPTVVAKAREHYGCDDVTYVELEDAGGSGTMGSHWKIRNAQDELMAGI"
        "SGVAYYTSLTLSAFEDLGYYKANYSNAETMKWGKDVGCAFLTGKCVVDNVTQFPSMYCDK"
        "DENVYRCHTARLNLGSCEVTDYTFDLPDYLQYFTVPSVGGSADYYDYCPYIVRSPIGSCT"
        "QAASSASPFVSAFNTFSMASRCIDGTFTPKSTGGATVTAHLGMCTNVACNTADKTYSIQV"
        "YGNGAYIPCTPGATISLDTVSDAFEAGGNITCPPYLEVCQSNVKGAMDYESMTNSGSGSS"
        "RPAPVEPSGSGSGSSAATTAPSPTRDGSAAADRIAPRTAAVALLALAVAAACV"
    ),
    "GP63_LEIMA_P08148": (  # Leishmania major, LmjF.10.0460
        "MSVDSSSTHRRRCVAARLVRLAAAGAAVTVAVGTAAAWAHAGALQHRCVHDAMQARVRQS"
        "VADHHKAPGAVSAVGLPYVTLDAAHTAAAADPRPGSARSVVRDVNWGALRIAVSTEDLTD"
        "PAYHCARVGQHVKDHAGAIVTCTAEDILTNEKRDILVKHLIPQAVQLHTERLKVQQVQGK"
        "WKVTDMVGDICGDFKVPQAHITEGFSNTDFVMYVASVPSEEGVLAWATTCQTFSDGHPAV"
        "GVINIPAANIASRYDQLVTRVVTHEMAHALGFSGPFFEDARIVANVPNVRGKNFDVPVIN"
        "SSTAVAKAREQYGCDTLEYLEVEDQGGAGSAGSHIKMRNAQDELMAPAAAAGYYTALTMA"
        "IFQDLGFYQADFSKAEVMPWGQNAGCAFLTNKCMEQSVTQWPAMFCNESEDAIRCPTSRL"
        "SLGACGVTRHPGLPPYWQYFTDPSLAGVSAFMDYCPVVVPYSDGSCTQRASEAHASLLPF"
        "NVFSDAARCIDGAFRPKATDGIVKSYAGLCANVQCDTATRTYSVQVHGSNDYTNCTPGLR"
        "VELSTVSNAFEGGGYITCPPYVEVCQGNVQAAKDGGNTAAGRRGPRAAATALLVAALLAV"
        "AL"
    ),
}
FRIENDLY_NAME = {
    "GP63_CRIFA_Q06031": "GP63 C. fasciculata (CFAC1_040011900)",
    "GP63_LEIMA_P08148": "GP63 L. major (LmjF.10.0460)",
}

# ---------------------------------------------------------------------------
# 2. Pull CbWHA1 GP63 sequences
# ---------------------------------------------------------------------------
def load_fasta_seq(path, target):
    """Return one sequence (by header, first whitespace-delimited token) from a FASTA file."""
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


def get_gp63_genes(gff_path):
    """base_gene_id -> [(mRNA/protein_id, seqid, start, end, strand), ...] from CDS lines."""
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if "\tCDS\t" not in line:
                continue
            if "gp63" not in line.lower() and "leishmanolysin" not in line.lower():
                continue
            fields = line.rstrip("\n").split("\t")
            parent = re.search(r"Parent=([^;]+)", fields[8]).group(1)
            base = re.match(r"(CbWHA1\.\d+\.\d+)", parent).group(1)
            genes.setdefault(base, []).append(
                (parent, fields[0], int(fields[3]), int(fields[4]), fields[6])
            )
    return genes


def build_cb_sequences(gff_path, pep_path, genome_path, exclude=frozenset()):
    pep_seqs = load_pep(pep_path)
    genes = get_gp63_genes(gff_path)
    genome_cache = {}

    def get_chr(chrom):
        if chrom not in genome_cache:
            genome_cache[chrom] = load_fasta_seq(genome_path, chrom)
        return genome_cache[chrom]

    result = {}
    for base, entries in genes.items():
        if base in exclude:
            continue
        protein_id, chrom, start, end, strand = entries[0]
        if protein_id in pep_seqs:
            result[base] = pep_seqs[protein_id].rstrip("*")
        else:
            # TransDecoder didn't emit a peptide record for this one
            # (usually a 5'-partial ORF call) — translate the annotated
            # CDS span directly from the genome instead.
            sub = get_chr(chrom)[start - 1 : end]
            if strand == "-":
                sub = str(Seq(sub).reverse_complement())
            result[base] = str(Seq(sub).translate()).rstrip("*")
    return result


# ---------------------------------------------------------------------------
# 3. Align with FAMSA
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


# ---------------------------------------------------------------------------
# 4. NJ tree + bootstrap support
# ---------------------------------------------------------------------------
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
# 5. Plotting
# ---------------------------------------------------------------------------
def color_for(name):
    if name is None:
        return "black"
    if name.startswith("CbWHA1.07."):
        return "#1b9e77"  # chr07 tandem array — closest to C. fasciculata
    if name.startswith(("CbWHA1.24.", "CbWHA1.27.")):
        return "#d95f02"  # dispersed, more divergent clade — closest to L. major
    if name.startswith("GP63_CRIFA"):
        return "#1b9e77"
    if name.startswith("GP63_LEIMA"):
        return "#d95f02"
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

    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(
        support_tree,
        axes=ax,
        do_show=False,
        label_func=label_func,
        label_colors=lambda lbl: label_color.get(lbl, "black"),  # NB: keyed by
        # the exact rendered label string, which Bio.Phylo prefixes with a
        # leading space — build label_color from label_func's return value,
        # not from clade.name, or every lookup silently falls back to black.
        branch_labels=branch_label,
        show_confidence=False,
    )

    for t in ax.texts:
        txt = t.get_text()
        if not txt:
            continue
        if any(ch.isalpha() for ch in txt):
            t.set_fontweight("bold")
            t.set_fontsize(10)
        else:
            # bootstrap value: bold red, nudged up so it clears the branch line
            t.set_fontsize(8)
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
    cb_seqs = build_cb_sequences(GFF, PEP, GENOME, exclude=EXCLUDE)
    print(f"CbWHA1 GP63 sequences: {len(cb_seqs)}")

    all_seqs = {**cb_seqs, **REFERENCES}
    write_fasta(all_seqs, "combined.fa")

    aligned = align_sequences(all_seqs)
    write_fasta(aligned, "aligned.fa")

    tree = build_bootstrap_tree("aligned.fa", n_replicates=500)
    Phylo.write(tree, OUT_NWK, "newick")
    # NB: Bio.Phylo's Newick writer merges confidence values into internal
    # node *names* on write (e.g. "Inner22" + 100.0 -> "Inner22100.00"), so a
    # tree reloaded from this .nwk will have clade.confidence == None again.
    # Plot from the in-memory `tree` object (as done below), not a reloaded one.

    title = (
        "C. bombi GP63/leishmanolysin family vs. C. fasciculata & L. major references\n"
        "NJ tree, p-distance, 500 bootstrap replicates (red = support at internal nodes)\n"
        "CbWHA1.31.003560 excluded (partial/divergent alignment)"
    )
    plot_tree(tree, OUT_PNG, title)
    print("done")
