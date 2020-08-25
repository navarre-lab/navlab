from Bio import SeqIO
from Bio.SeqUtils import GC


def convert_gbk_to_fna(filename):
    """
    Takes Genbank file (.gbff) and converts it to FASTA nucleic acid file (.fna)
    """

    input_handle = open(filename, "r")
    output_handle = open(filename.split(".")[0] + ".fna", "w")

    index = 0
    rna_index = 0
    seq_accession_index = ""

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        print("Dealing with GenBank record %s" % seq_record.id)
        seq_accession_prefix = seq_record.id.split("_")[0]

        if seq_accession_prefix == "NC":
            index = 0
        elif seq_accession_prefix == "NZ":
            new_seq_accession_index = seq_record.id[:7]
            if new_seq_accession_index != seq_accession_index:
                index = 0
        # setting a base locus number in case the gbk file does not automatically have one for its CDS features.
        locus_number = 5
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                index += 1
                seq_feature_sequence = line_format(str(seq_feature.extract(seq_record).seq))
                GC_content = str(GC(seq_feature.extract(seq_record).seq))
                location = str(seq_feature.location)
                if "locus_tag" in seq_feature.qualifiers:
                    locus_tag = seq_feature.qualifiers["locus_tag"][0]
                else:
                    title = (seq_record.annotations["accessions"][0])
                    locus_tag = (title + "_" + str(locus_number).zfill(5))
                    seq_feature.qualifiers["locus_tag"] = [locus_tag]
                    locus_number += 5
                if "gene" in seq_feature.qualifiers:
                    gene = seq_feature.qualifiers["gene"][0]
                else:
                    gene = "none"
                if "protein_id" in seq_feature.qualifiers:
                    protein_id = seq_feature.qualifiers["protein_id"][0]
                    protein_accession = "_" + protein_id
                if "product" in seq_feature.qualifiers:
                    product = seq_feature.qualifiers["product"][0]
                if "pseudo" in seq_feature.qualifiers:
                    protein_id = "PSEUDOGENE"
                    protein_accession = ""
                output_handle.write(
                    ">lcl|%s_cds%s%s [gene=%s] [locus_tag=%s] [protein=%s] [protein_id=%s] [location=%s] [GCpct=%s] [gbkey=CDS]\n%s\n" % (
                        seq_record.id,
                        protein_accession,
                        str("_" + str(index)),
                        str(gene),
                        locus_tag,
                        product,
                        protein_id,
                        location,
                        '{:.4}'.format(GC_content),
                        seq_feature_sequence))

            if seq_feature.type in ["tRNA", "rRNA", "ncRNA", "tmRNA"]:
                rna_index += 1
                seq_feature_sequence = line_format(str(seq_feature.extract(seq_record).seq))
                if "locus_tag" in seq_feature.qualifiers:
                    locus_tag = seq_feature.qualifiers["locus_tag"][0]
                if "product" in seq_feature.qualifiers:
                    product = seq_feature.qualifiers["product"][0]
                if "gene" in seq_feature.qualifiers:
                    gene = seq_feature.qualifiers["gene"][0]
                else:
                    gene = "none"
                location = str(seq_feature.location)
                GC_content = str(GC(seq_feature.extract(seq_record).seq))
                protein_id = product
                protein_accession = ""
                output_handle.write(
                    ">lcl|%s_%s%s [gene=%s] [locus_tag=%s] [product=%s] [location=%s] [GCpct=%s] [gbkey=%s]\n%s\n" % (
                        seq_record.id,
                        seq_feature.type.lower(),
                        str("_" + str(rna_index)),
                        str(gene),
                        locus_tag,
                        product,
                        location,
                        '{:.4}'.format(GC_content),
                        seq_feature.type,
                        seq_feature_sequence))

    output_handle.close()
    return


def convert_fna_to_faa(filename):
    """
    Takes a nucleotide FASTA file 'filename' and converts it to a regular Protein FASTA file while maintaining the
    info in the CDS header.
    """
    output_handle = open(filename.split(".")[0] + ".faa", "w")
    for seq_record in SeqIO.parse(filename, "fasta"):
        if get_gb_key(seq_record.description) == "CDS":
            output_handle.write(">" + seq_record.description[4:] + "\n")
            # the [:-1] removes the * for the stop codon from the sequence record
            output_handle.write(line_format(str(seq_record.seq.translate()[:-1])) + "\n")
    output_handle.close()
    print("Done")


def get_gb_key(description_line):
    """
    Retrieves Genbank key from string "description_line".
    """
    # remove the brackets at the end of the line
    descriptors = str(description_line[1:-1])
    # further split the tags into individual items by splitting at the "] [" between each tag
    descriptor_list = descriptors.split("] [")
    for item in descriptor_list:
        if "gbkey=" in str(item):
            descriptor = (str(item)[6:])
            return descriptor
    descriptor = "none"
    return descriptor


def get_locus_id(seq_record):
    """
    This functions retrieves the locus id from the sequence record seq_record.
    """
    for seq_feature in seq_record.features:
        if seq_feature.type == "CDS" and "locus_tag" in seq_feature.qualifiers:
            locus_id = seq_feature.qualifiers["locus_tag"][0]
            return locus_id


def line_format(sequence):
    """
    This function returns a FASTA sequence with 80 characters per line
    """
    newline_list = []
    while sequence:
        newline = sequence[0:80]
        newline_list.append(newline)
        sequence = sequence[80:]
    output = "\n".join(newline_list)
    return output


def dictmaker(filename):

    input_handle = open(filename, "r")
    locus_tag_set = set()
    locus_descriptor_list = []
    dict_entry = ""

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        description = seq_record.description
        accession = (" | ".join(seq_record.dbxrefs))
        locus_id = get_locus_id(seq_record)
        locus_tag = get_locus_tag(locus_id)
        if (locus_tag is not None) and (locus_tag not in locus_tag_set):
            locus_tag_set.add(locus_tag)
            dict_entry = str(locus_tag + " | " + accession + " | " + description)
            locus_descriptor_list.append(dict_entry)
        else:
            dict_entry = accession + " | no locus_tag found"
    return dict_entry
