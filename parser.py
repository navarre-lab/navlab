"""
MODULE FOR PARSING XML TO TAB
REQUIRES A DICT FILE AND AN XML BLAST OUTPUT
the functions below extract info from the "description" lines of the XML file
...called from the BioPython parser by tags like "hsp.hit_description"
"""

"""
----- HELPER FUNCTIONS -----
"""


def get_descriptor_list(line) -> list:
    """
    This function removes the brackets at the the end of str "line" and further splits the tags into individual
    items by splitting at the "] [" between each tag
    """
    descriptors = str(line[1:-1])
    descriptor_list = descriptors.split("] [")
    return descriptor_list


def get_gene(line):
    """

    """
    descriptor_list = get_descriptor_list(line)
    for i in descriptor_list:
        if "gene=" in str(i):
            descriptor = (str(i)[5:])
            return descriptor
    return "none"


def get_location(line):
    """

    """
    descriptor_list = get_descriptor_list(line)
    for i in descriptor_list:
        if "location=" in str(i):
            descriptor = (str(i)[9:])
            return descriptor
    return "none"


def get_GC(line):
    """


    """
    descriptor_list = get_descriptor_list(line)
    for i in descriptor_list:
        if "GCpct=" in str(i):
            descriptor = (str(i)[6:])
            return descriptor
    return "none"


def get_locus_id(line):
    """
    This will give you the full locus id of the object from the description line [locus_id=XXXX####]
    """
    descriptor_list = get_descriptor_list(line)
    for i in descriptor_list:
        if "locus_tag=" in str(i):
            locus_id = (str(i)[10:])
            return locus_id
    locus_id = "no_locus_id"
    return locus_id


def get_locus_tag(line):
    """
    This will give you the locus tag at the front of the locusid.  There are two general formats of locus tags.
    """
    import re

    # sometimes we pass items like "Y75_RS14370|WP_000220066.1|0.84" and we only want the front part.
    line = str(line.split("|")[0])
    if "_" in line:
        # this works for locus tags formatted like STM14_#######
        locus_tag = str(line.split("_")[0])
        return locus_tag
    else:
        # this works for locus tags formatted like STM#### (no hyphen...letters before a string of #s)
        locus = re.match("\D+", line)
        locus_tag = locus.group(0)
        return locus_tag


def get_protein(x):
    descriptor_list = get_descriptor_list(x)
    for i in descriptor_list:
        if "protein=" in str(i):
            descriptor = (str(i)[8:])
            return descriptor
        if "product=" in str(i):
            descriptor = (str(i)[8:])
            return descriptor
        if "gbkey=tRNA" in str(i):
            return "tRNA"
    return "none"


def get_protein_id(x):
    descriptor_list = get_descriptor_list(x)
    for i in descriptor_list:
        if "protein_id=" in str(i):
            descriptor = (str(i)[11:])
            return descriptor
        if "pseudo=true" in str(i):
            descriptor = "PSUEDOGENE"
            return descriptor
        if "product=" in str(i):
            descriptor = (str(i)[8:])
            return descriptor
        if "gbkey=tRNA" in str(i):
            return "tRNA"
    return "none"


def parse_cluster(infile):
    from Bio import SeqIO
    import nav

    with open(infile) as f:
        lines = [line.rstrip() for line in f]
        i = 1  # start an index.  Skip the very first line in the file.
        member_list = []
        cluster_list = []
        while i < len(lines):
            first_line = lines[i]
            second_line = lines[i + 1]
            if second_line[0] == ">":  # if there are two > in a row it means a new cluster has started
                i += 1  # advance the "frame" by one so you get back on track.
                cluster_list.append(member_list)  # add previously assembled cluster to a growing list of clusters.
                member_list = []  # reset the member_list to get ready to add sequences to.
            else:  # otherwise the first line should be the definition line and the second line should be the sequence.
                member_list.append(
                    first_line + "\n" + nav.line_format(second_line))  # append the sequence to the growing cluster.
                i += 2  # skip ahead two lines
        cluster_list.append(member_list)  # the last cluster in the cue now needs to be added to the cluster list

        # now sort the cluster_list by size.  The most abundant proteins appear first in the list.
        cluster_list = sorted(cluster_list, key=len, reverse=True)  # sort the cluster list by the size of each cluster.

        output_handle = open("cluster_summary.faa", "w")
        output_handle2 = open("representative_proteins.faa", "w")
        output_handle3 = open("faa/CLUSTER.faa", "w")

        # now we iterate through the list of clusters and rename them
        p = 1
        output_handle.write("#Number of clusters: " + str(len(cluster_list)) + "\n")
        for item in cluster_list:
            gene_set = set()
            protein_set = set()
            output_handle.write("#CLUSTER_" + str(p) + " [members= " + str(len(item)) + "]\n")
            output_handle.write((str("\n".join(item)) + "\n"))  # write the members of each cluster

            # now we change the headers to incorporate *all* the info for a given cluster
            # because some proteins in the cluster have a different annotation than others
            for cluster_member in item:
                member_fragments = cluster_member.split(" [")
                for n, x in enumerate(member_fragments):
                    if "gene=" in x:
                        gene_set.add(str(x[5:-1].split("_")[0]))
                    if "protein=" in x:
                        protein_set.add(str(x[8:-1]))
            if len(gene_set) > 1:
                gene_set.discard("none")
            gene_names = ", ".join(gene_set)
            protein_names = ", ".join(protein_set)
            representative = item[0]
            representative_header = representative.split("\n")[0]
            representative_sequence = representative.split("\n")[1:]
            rep_header_fragments = representative_header.split(" [")
            for n, x in enumerate(rep_header_fragments):
                if "gene=" in x:
                    rep_header_fragments[n] = "gene=" + gene_names + "]"
                if "protein=" in x:
                    rep_header_fragments[n] = "protein=" + protein_names + "]"
            new_representative_header = " [".join(rep_header_fragments)
            output_handle2.write((str(new_representative_header) + "\n"))
            output_handle2.write((str("\n".join(representative_sequence) + "\n")))
            p += 1
        output_handle.close()
        output_handle2.close()

        # AT THIS POINT WE READ THE REPRESENTATIVE PROTEINS AND GIVE THEM CLUSTER IDS
        q = 1
        for member in SeqIO.parse("representative_proteins.faa", "fasta"):
            fragment = member.id.split("_cds_")[1]
            fragment = fragment.split("_")[0:2]
            wp_name = '_'.join(fragment)
            member.id = str("CLUSTER_" + str(q)) + "_cds_" + wp_name
            description_list = member.description.split(" ")[1:]
            for n, x in enumerate(description_list):
                if "locus_tag=" in x:
                    description_list[n] = "[locus_tag=CLUSTER_" + str(q) + "]"
            description_line = " ".join(description_list)
            description_line2 = (">" + member.id + " " + description_line)
            output_handle3.write(description_line2 + "\n")
            output_handle3.write(nav.line_format(str(member.seq) + "\n"))
            q += 1

    output_handle3.close()


def parse_xml(infile, outfile, strain_dictionary, pct_id_threshold):
    import numpy
    from Bio import SearchIO

    output_tab_file = open(outfile, "w")

    # Below we set the percent identity required for inclusion in our list.  The default value is set to 70%.
    try:
        pct_id_threshold = float(pct_id_threshold) / float(100)
    except TypeError:
        pct_id_threshold = float(0.7)  # set the default to 70%

    # Here we read the "dictionary" file and make our list.
    # The format should be: Locus tag | BioSampleID | species and strain

    strain_list = []
    strain_locus_tag_list = []

    with open(strain_dictionary, "r") as f:
        for line in f:
            strain_list.append(line.rstrip())
            dict_locus_tag = line.split(" | ")[0]
            strain_locus_tag_list.append(dict_locus_tag)

    strains = "\t".join(strain_list)
    # HERE WE IMPORT DATA FROM THE FILE WITH THE NAME GIVEN FROM THE COMMAND LINE
    q_results = SearchIO.parse(infile, 'blast-xml')

    output_tab_file.write("LOCUS ID\tACCESSION\tGENE NAME\tDNA ACCESSION\tLOCATION\tPROTEIN ID\tGENE SIZE\tPCT_GC\t"
                          "NUM OF HITS\tSTRAINS WITH HITS\tAVG % ID\tDESCRIPTION\t" + strains + "\n")

    for q_result in q_results:
        query_accessions = str(q_result.id)
        query_DNA_accession = str(query_accessions.split("_cds_")[0])  # retrieves only the nucleotide portion
        query_locus_id = get_locus_id(q_result.description)
        query_protein_id = str(get_protein_id(q_result.description))
        query_location = str(get_location(q_result.description))
        query_gene_name = str(get_gene(q_result.description))
        query_GC = str(get_GC(q_result.description))
        percent_id_list = []
        locus_list = set()
        descriptions = set()
        hit = str()
        hsp = str()

        for hit in q_result:
            for hsp in hit:
                pct_id = float(hsp.ident_num) / float(hit.seq_len)
                percent_id = str(pct_id)[0:4]
                if float(pct_id) > pct_id_threshold:
                    hit_accession = str(hsp.hit_id)
                    percent_id_list.append(pct_id)
                    locus_id = get_locus_id(hsp.hit_description)
                    protein = get_protein(hsp.hit_description)
                    protein_id = get_protein_id(hsp.hit_description)
                    hit_info = str(hit_accession) + "|" + str(locus_id) + "|" + str(protein_id) + "|" + percent_id
                    locus_list.add(hit_info)
                    descriptions.add(protein)
        if len(percent_id_list) >= 1:
            average_percent_id = str(numpy.mean(percent_id_list))[0:5]
        else:
            average_percent_id = "N/A"

        if hit == "":
            descriptions.add("no hits that exceeded e-value")
        else:
            descriptions.add(get_protein(hsp.fragment.query.description))

        description_list = list(descriptions)
        description_list.sort()
        list_of_descriptions = ", ".join(description_list)

        hit_locus_list = list(locus_list)
        hit_locus_list.sort()

        strain_hits = []
        for i in strain_locus_tag_list:
            match_list = []
            for j in hit_locus_list:
                j_locus = str(j.split("|")[1])
                j_locus_tag = get_locus_tag(j_locus)
                if j_locus_tag == i:
                    match_list.append(j)
            if len(match_list) >= 1:
                matches = "[" + str(len(match_list)) + "]|" + str(", ".join(match_list) + "\t")
            else:
                matches = str("*\t")
            strain_hits.append(matches)
        number_of_strains_with_hits = len(strain_hits) - strain_hits.count("*\t")
        list_of_hit_loci = "".join(strain_hits)
        output_tab_file.write(str(query_locus_id) + "\t" + str(query_accessions) + "\t" + query_gene_name + "\t"
                              + query_DNA_accession + "\t" + query_location + "\t" + query_protein_id + "\t"
                              + str(q_result.seq_len) + "\t" + query_GC + "\t" + str(len(hit_locus_list)) + "\t"
                              + str(number_of_strains_with_hits) + "\t" + average_percent_id + "\t"
                              + list_of_descriptions + "\t" + list_of_hit_loci + "\n")
    return


def pseudo_remover(infile):
    from Bio import SeqIO
    import nav

    # make a new file to write out to.
    output_handle = open("faa/pseudofree.faa", "w")

    # define a new list variable, object_list, where we will put the sequence objects.
    object_list = []
    for seq_record in SeqIO.parse(infile, "fasta"):
        if len(seq_record.seq) > 30:  # we will only keep proteins larger than 60 AA
            object_list.append(seq_record)

    for item in object_list:
        if "PSEUDOGENE" not in item.description:
            output_handle.write(">" + item.description + "\n")
            output_handle.write(nav.line_format(str(item.seq) + "\n"))

    output_handle.close()
    print("Done")
