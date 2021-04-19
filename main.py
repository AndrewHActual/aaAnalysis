# coding=utf-8
# An attempt at extracting aa information from all the different sources

from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import re

# Snagging the translated gene sequences and locations from the genbank file
from pandas import DataFrame


def parse_genbank(genbank_file):
    # Give us that delicious genbank file as a SeqRecord
    with open(genbank_file) as handle:
        gbRecord = SeqIO.read(handle, "genbank")

    gene_aas = []

    # Pulling out only the CDS features and the associated qualifiers/locations we want
    if gbRecord.features:
        for feature in gbRecord.features:
            if feature.type != "CDS":
                continue
            gene_aas.append([feature.qualifiers["gene"][0],
                             feature.location.start,feature.location.end,
                             feature.qualifiers["translation"][0]])

    # Name the columns as we like
    parseddf = pd.DataFrame(gene_aas,
                            columns=['gene_name','nt_start','nt_end','translation'])

    return parseddf

#   Light parsing of a provided stripped GFF file into a dataframe
#   from a file path
def parse_strippedGFF(strippedGFF_file):
    with open(strippedGFF_file) as handle:
        df = pd.read_csv(handle,
                         sep="\t")
    return df

# IRMA/Dais ribosome output parser
#   Take in dais-ribosome output
#   Pull in ref aa seq for each gene from the genbank record
#   Parse .ins and .del files
#   Generate aa subs from aligned aa qry seq
#   Return table of aa changes
def irma_aavariants(dirpath):

    # Parse out the dais ribosome output into dfs
    seqheader = ['ID', 'C_type', 'Ref_ID', 'Protein', 'VH', 'AA_seq', 'AA_aln', 'CDS_id',
                 'Insertion', 'Shift_Insert', 'CDS_seq', 'CDS_aln', 'Query_nt_coordinates',
                 'CDS_nt_coordinates']
    insheader = ['ID','C_type','Ref_ID','Protein','Upstream_aa','Inserted_nucleotides',
                 'Inserted_residues','Upstream_nt','Codon_shift']
    delheader = ['ID','C_type','Ref_ID','Protein','VH','Del_AA_start','Del_AA_end',
                 'Del_AA_len','In_frame','CDS_ID','Del_CDS_start','Del_CDS_end','Del_CDS_len']

    qryseqdf = dir_to_dataframe_plus(dirpath,".seq","\t",".",seqheader)
    qryinsdf = dir_to_dataframe_plus(dirpath,'.ins','\t',".",insheader)
    qrydeldf = dir_to_dataframe_plus(dirpath,'.del','\t','.',delheader)

    # Parse the reference genbank file into a df
    refdf = parse_genbank("ref/NC_045512.gb")

    # Limit ourselves to these genes
    onlythesegenes = set(gene.casefold() for gene in ("S"))

    variants_df = pd.DataFrame(columns=["seq_name","gene_name","RV_ref_aa","RV_aa_pos","RV_qry_aa"])

    #   Extract substitutions from qryseqdf
    #   Go through each seq + gene in qryseqdf
    #       Protein matching gene from list
    #       C_type not "UNRECOGNIZABLE", as happens with bad assemblies
    #       AA_seq not "\\N", again couldn't come up with enough to translate
    #   Get the reference sequence and aligned aa seq for that sequence's gene/protein
    #       And do a literal string comparison for substitutions
    #   Output a dataframe with:
    #       Seq_name, gene name, ref aa, aa pos, and qry aa
    #   TODO: This is long enough to be its own method
    for index, qryCDSrec in qryseqdf.iterrows():
        if not qryCDSrec['Protein'].casefold() in onlythesegenes:
            continue
        if qryCDSrec['C_type'] == 'UNRECOGNIZABLE':
            continue
        if qryCDSrec['AA_seq'] == '\\N':
            continue

        # Get the query aligned aa seq for this gene
        qryseqname = qryCDSrec['seq_name']
        qrygenename = qryCDSrec['Protein']
        qryseq = qryCDSrec['AA_aln']

        # To get the ref aa seq, do a case insensitive search for our query gene name,
        #   and return a Boolean Series
        refaabools = refdf['gene_name'].str.contains(qrygenename, case=False)
        # Pandas just wants to keep returning series and dataframes all the way down unless you use this
        refseq = refdf[refaabools]['translation'].iloc[0]

        # String comparison between two
        mutation_list = mutation_list_strict(refseq,qryseq)

        # Turn list of mutations (ref aa, aa pos, qry aa) into a df
        # TODO: This returns a lot of substitutions that are X's (ambiguity aa). Should at least offer
        #   a flag to ignore Xs
        mutation_df = pd.DataFrame(mutation_list,
                                   columns=["RV_ref_aa","RV_aa_pos","RV_qry_aa"])
        # Then brute force add the sequence name and gene/protein
        mutation_df.insert(0,"gene_name",qrygenename)
        mutation_df.insert(0,"seq_name",qryseqname)

        # Append that to our variants df
        variants_df = variants_df.append(mutation_df)

    #  Extract insertions from qryinsdf
    #   Output a df with:
    #   Seq_name, gene name, ref aa (or -), upstream aa, pos, qry aa
    #   Ugh, we're fully committed to separating out each aa in a multi-aa insertion too
    #   Or not... We're just reporting them directly.
    #   TODO: Not tested even a little
    for index, qryinsrec in qryinsdf.iterrows():
        if not qryinsrec['Protein'].casefold() in onlythesegenes:
            continue
        ins_df = qryinsrec[['seq_name','Protein','Upstream_aa','Inserted_residues']]
        ins_df.insert(2,'RV_ref_aa','-')
        ins_df.columns = ["seq_name","gene_name","RV_ref_aa","RV_aa_pos","RV_qry_aa"]

        # Append that to our variants df
        variants_df = variants_df.append(ins_df)

    # Extract deletions from qrydeldf
    #   Output a df with:
    #   Seq_name, gene name, ref aa, upstream aa, pos, qry aa (-)
    #   Break up longer deletions into single aa deletions
    #   TODO: This is long enough to be its own method
    for index, qrydelrec in qrydeldf.iterrows():
        if not qrydelrec['Protein'].casefold() in onlythesegenes:
            continue

        # First, get the ref aa seq, do a case insensitive search for our query gene name,
        #   and return a Boolean Series
        refaabools = refdf['gene_name'].str.contains(qrydelrec['Protein'], case=False)
        # Pandas just wants to keep returning series and dataframes all the way down unless you use this
        refseq = refdf[refaabools]['translation'].iloc[0]

        # Iterate over the query's deletion range to output
        #   each ref aa, position, and qry gap
        del_list = []
        delrange = range(qrydelrec['Del_AA_start'], qrydelrec['Del_AA_end'] + 1)
        for pos in delrange:
            del_list.append([qrydelrec['seq_name'],qrydelrec['Protein'],refseq[pos-1],pos,"-"])

        del_df = pd.DataFrame(del_list)
        del_df.columns = ["seq_name", "gene_name", "RV_ref_aa", "RV_aa_pos", "RV_qry_aa"]

        # Append that to our variants df
        variants_df = variants_df.append(del_df)

    # mutationdf = mutation_list_strict("MAIVBOIF","MAQV.PIFRT")

    return variants_df

#   Replicate udf-bioutils:Mutation_List_Strict used by RVB to analyze Dais-Ribosome output
#       Takes two strings as input
#   Steps
#       Gets shorter length of the two strings
#       Sets seq1, seq2, and a string buffer (for changes?)
#       For indexes 0 to shorter length
#           if seq1[i] != seq2[i]
#               if seq1[i] != seq2[i] as upper cases (whyyyyy)
#                   if seq1[i] != '.' && seq2[i] != '.' (dais-ribo's missing alignment data, non-standard)
#                       if buffer.length > 0
#                           add seq1[i], i+1 (pos), seq2[i]
#                       else
#                           prep buffer...
#       returns buffer
#   What could we do?
#       seq1/seq2 to upper
#       minlen = len(seq1) if len(seq1) < len(seq2) else len(seq2)
#       changes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
#   Input: reference aa seq (as string), query aa seq (as string)
#   Output: df of just aa changes from ref
#               ref aa, pos, qry aa
def mutation_list_strict(refaaseq, qryaaseq):
    # TODO: check for aa content in both
    #   Kick error if not...

    # Convert both seq to upper case, just in case
    rseq = refaaseq.upper()
    qseq = qryaaseq.upper()

    # Figure out which sequence is shorter based for range function
    minlen = len(rseq) if len(rseq) < len(qseq) else len(qseq)

    # Compare differences in strings and return those aas and positions
    #   Also adds 1 to the position to shift from 0 indexed to 1 indexed for aa
    changes = [(rseq[i],i+1,qseq[i]) for i in range(minlen) if rseq[i] != qseq[i]]

    return changes

# iVar variants
#   Pull in table to link gene names and nt coords to CDS-ids
#   Pull in ivar variant files from directory
#   Convert from nt coord to aa coord
#   Add gene names
#   Return table of aa changes
def ivar_aavariants(dirpath):
    rawdf = dir_to_dataframe_plus(dirpath, ".variants.tsv", "\t", ".")
    CDSdf = parse_strippedGFF("ref/MN908947.3.tsv")

    #   Inner join to get annotations in
    annotdf = pd.merge(rawdf[["seq_name","POS","GFF_FEATURE","REF_AA","ALT_AA"]], CDSdf,
                       left_on="GFF_FEATURE", right_on="CDS_id")
    # And do a dirty calculation of aa position
    #   Just the position of the aa change minus the start position (from ref) of its gene
    #   divided by 3 and plus 1 for the Hebrew g-d
    annotdf["aa_pos"] = (annotdf["POS"] - annotdf["start_nt_pos"])/3
    annotdf["aa_pos"] = annotdf["aa_pos"].apply(np.floor)
    annotdf["aa_pos"] = annotdf["aa_pos"] + 1

    # Take only the columns we need
    annotdf = annotdf[["seq_name","gene_name","REF_AA","aa_pos","ALT_AA"]]
    annotdf.columns = ["seq_name", "gene_name", "IV_ref_aa", "IV_aa_pos", "IV_qry_AA"]

    return annotdf

# NextClade variants
#   Input: directory path
#   Output: dataframe with samples and aa variants by gene names
#   Pull in all nextclade files from directory
#   Nextclade has a file per sample
#       Each sample is in a single row with all of the aa sub/del/ins in one cell
#       Parse out sample name, and aa subs/ins/del parsed into gene name, ref aa, aa pos, and qry aa
def nextclade_aavariants(dirpath):
    rawdf = dir_to_dataframe_plus(dirpath,".csv",";","_")
    parsedvariantslist = []

    for index, row in rawdf.iterrows():
        # Splitting out the substitutions and deletions, should they exist, into a list for each sequence
        subsdelslist = []
        if type(row["aaSubstitutions"]) == str:
            subsdelslist = subsdelslist + row["aaSubstitutions"].split(",")
        if type(row["aaDeletions"])  == str:
            subsdelslist = subsdelslist + row["aaDeletions"].split(",")
            # Fab note: Using append instead of + in the above
            #   yielded weird behavior when iterating below

        # Extract the sequence name from the sequence name column that NextClade assigns
        seqName = row["seq_name"]

        # Going through list of subsdels and parse them into tuples,
        #   then into a list with the sequence name
        for subdel in subsdelslist:
            # Grouping out the line for
            #   gene name (alphanumeric), ref aa (incl symbols), position integer, qry aa
            matches = re.match(r'(\w*):(\D*)(\d*)(\D*)', subdel)
            parsed = [seqName] + list(matches.groups())
            parsedvariantslist.append(parsed)

    # TODO: Probably need to make this column def a passed argument
    parseddf = pd.DataFrame(parsedvariantslist,
                            columns=["seq_name", "gene_name", "NC_ref_aa", "NC_aa_pos", "NC_qry_aa"])
    return(parseddf)

#   Extract lightly processed .csv containing s gene aa variants from the SC2 Overarching Sample file
#       As of 4/18/21
#   Input CSV Rows will have the sample ID, lineage call (out of date), and amino acid substitutions
#       each followed by another row that is just the aa substitutions which can be discarded
#   Returns df with one aa variant per line
#       NOTE: Follows a different column order from the other aa variant dfs
#       seq name, (s)gene name (static), aa pos, ref aa, qry aa
def rvb_aavariants(dirpath):
    raw_df = pd.read_csv(dirpath,
                         sep = "\t")
    # Strip off those null named rows which were just part of the terrible excel processing
    raw_df = raw_df[raw_df['Strain_ID'].notnull()]

    # Break down into table with strain ID and parsed out s-gene substitutions
    #   strain ID, aa pos, ref_aa_rvb, qry_aa_rvb
    #   We're ignoring the hell out of that lineage
    #   Every single gene in these list are the S genes
    parsed_list = []
    genename = 'S'

    for index, row in raw_df.iterrows():
        subst_raw = row['s_subst']
        
        # Drop a mostly empty row for samples without substitutions
        if subst_raw == "None":
            parsed_list.append([row['Strain_ID'], genename, 0, None, None])
            continue
            
        # For each subst in s_subst...
        #   Split out the subst and parse each out into a fresh row
        subst_list = subst_raw.split(",")
        for subst in subst_list:
            # Added something to get rid of white space characters because Excel sucks
            matches = re.match(r'(\D-*)(\d*)(\D-*)', subst.strip())
            matches_array = list(matches.groups())

            # We put the aa pos right after the name to make the final table look nicer
            parsed = [row['Strain_ID'], genename, matches_array[1], matches_array[0], matches_array[2]]
            parsed_list.append(parsed)

    parsed_df: DataFrame = pd.DataFrame(parsed_list,
                             columns=["seq_name", "gene_name", "aa_pos_RVB", "ref_aa_RVB", "qry_aa_RVB"])
    return parsed_df

# Directory to not so raw dataframe, but better
#   Inputs: directory path string
#           file suffixes (string or tuple of strings),
#           separator symbol (string)
#           character to trim beyond for sequence name
#               like with ".", it trims off suffixes after the first "."
#           column header names list [optional]
#   output: pandas dataframe with added first column as sequence name
#   Generate raw dataframe from directory of headered files of the specified suffixes
#   Header detection relies on pandas' read_csv default behavior
def dir_to_dataframe_plus(dirpath, suffixes, sep, trimcharacter, columnnames=None):
    # Avoiding that mutable default argument error...
    if columnnames is None:
        columnnames = []

    rowlist = []

    with os.scandir(dirpath) as it:
        for direntry in it:
            if not direntry.name.endswith(suffixes):
                continue
            if not direntry.is_file():
                continue

            if columnnames:
               rows = pd.read_csv(direntry,
                                  sep=sep,
                                  names=columnnames)
            else:
                rows = pd.read_csv(direntry,
                                   sep=sep)

            # TODO: Kind of a terrible hack for now to get the sequence name
            #   Just trims off everything after the trimcharacter
            seqname = direntry.name.split(trimcharacter)[0]
            rows.insert(0,"seq_name",seqname)

            rowlist.append(rows)

    #   Pandas magic that turns a list into a df
    df = pd.concat(rowlist)

    return df

if __name__ == '__main__':
#    parse_genbank("ref/NC_045512.gb")
#    parse_strippedGFF("ref/MN908947.3.tsv")
#    dir_to_dataframe_plus("M347-21-011/ivar_variant/",".variants.tsv","\t",".")

    nextclade_df = nextclade_aavariants("M3235-21-013/nextclade/")
    ivar_df = ivar_aavariants("M3235-21-013/ivar_variants/")
    irma_df = irma_aavariants("M3235-21-013/dais-ribo/")
    rvb_df = rvb_aavariants("ref/SC2_Overarching_lineage_sgenesubst.txt")

    # Fix aa column name for easy merging
    nextclade_df.columns = ["seq_name", "gene_name", "ref_aa_NC", "aa_pos", "qry_aa_NC"]
    ivar_df.columns = ["seq_name", "gene_name", "ref_aa_IV", "aa_pos", "qry_aa_IV"]
    irma_df.columns = ["seq_name", "gene_name", "ref_aa_DR", "aa_pos", "qry_aa_DR"]
    rvb_df.columns = ["seq_name", "gene_name", "aa_pos", "ref_aa_RVB", "qry_aa_RVB"]

    # TODO: Guess we have to convert the dtypes of all of our columns for this to be okay???
    #   https://stackoverflow.com/questions/39582984/pandas-merging-on-string-columns-not-working-bug/49912193
    #   Super gross typing problem. Perhaps unavoidable
    nextclade_df = nextclade_df.astype({"seq_name": str, "gene_name": str, "ref_aa_NC": str, "aa_pos": int, "qry_aa_NC":str})
    ivar_df = ivar_df.astype({"seq_name": str, "gene_name": str, "ref_aa_IV": str, "aa_pos": int, "qry_aa_IV":str})
    irma_df = irma_df.astype({"seq_name": str, "gene_name": str, "ref_aa_DR": str, "aa_pos": int, "qry_aa_DR":str})
    rvb_df = rvb_df.astype({"seq_name": str, "gene_name": str, "aa_pos": int, "ref_aa_RVB": str, "qry_aa_RVB":str})

    # Bunch of filters to just get us the S gene
    #   IRMA is already S gene filtered, but we have to filter out X ambiguity aa's
    nextclade_Sonly_df = nextclade_df[nextclade_df.gene_name == 'S']
    ivar_Sonly_df = ivar_df[ivar_df.gene_name == 'S']
    irma_Sonly_df = irma_df[irma_df.qry_aa_DR != 'X']

    # Actual 3 way merge
    merged_df = nextclade_Sonly_df.merge(right=ivar_Sonly_df,
                         how='outer',
                         on=['seq_name','gene_name','aa_pos'])

    merged_df = merged_df.merge(right=irma_Sonly_df,
                                how='outer',
                                on=['seq_name', 'gene_name', 'aa_pos'],
                                sort=True)

    # One more merge with the table of actual aa values from RVB
    #   But first we have to convert the seq_name into something that matches the RVB table
    merged_df['seq_name_short'] = merged_df['seq_name'].str.rpartition('-')[0]

    merged_df = merged_df.astype({"seq_name_short": str})

    merged_df = merged_df.merge(right=rvb_df,
                                how='left',
                                left_on=['seq_name_short', 'gene_name', 'aa_pos'],
                                right_on=['seq_name', 'gene_name', 'aa_pos'],
                                sort=True)
    print("bleh")
