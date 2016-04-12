# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:21:00 2015

@author: tyler
FRI-Failure Analysis Pipeline
Copyright Tyler Camp 2015
module: EvolutionaryAnalysis
"""
from Bio import SeqIO
from Bio import AlignIO
import os
import re
import sys
import subprocess
from multiprocessing import Process
from bioutils import *
from decimal import *
from collections import namedtuple
from Mutation import Mutation
from Analysis import Analysis
from AnalysisDefaults import *

class EvolutionaryAnalysis(Analysis):
    """EvolutionaryAnalysis: A subclass of Analysis used to perform evolution
    and stability based analyses.

    Unique instance variables:

    ancestor_seq_filepath: location on disk of ancestor sequence file
    sample_lst: list of sample sequence file paths"""

    def __init__(self, input_dir, output_dir):
        if input_dir[-1] == "/":
            self.input_dir = input_dir
        else:
            self.input_dir = input_dir + "/"
        if output_dir[-1] == "/":
            self.output_dir = output_dir
        else:
            self.output_dir = output_dir + "/"

        # Default values for instance variables
        self.ancestor_seq_path = None
        self.sample_lst = list()
        self.tmp_dir = self.init_temp_dir(type(self).__name__)
        print "tmp_dir:", self.tmp_dir
        self.template_dir = self.tmp_dir+"templates/"
        self.alignment_dir = self.tmp_dir+"alignments/"


        # Cutoff score for Sanger reads
        self.phred_cutoff = 10

        # Location of reference files
        self.mobile_ele_dir = MOB_ELE_DIR

        return

    def start_analysis(self):
        """Start the analysis for an Evolutionary study"""

        """Identify sequences to be analyzed"""
        file_list = list()
        for d, sd, fl in os.walk(self.input_dir):
            for f in fl:
                if ".gb" in f:
                    self.ancestor_seq_path = self.input_dir + f
                elif ".fasta" in f or ".ab1" in f:
                    self.sample_lst.append(self.input_dir + f)
                else:
                    print "Ignoring file: ", f

        """Create miscellaneous temporary directories"""
        os.mkdir(self.template_dir)
        os.mkdir(self.alignment_dir)

        # Begin work flow

        # Create template files for alignment software
        self.seqproc()

        # Align the templates
        self.create_alignment_driver()

        # Analyze the alignments for mutations
        self.mutation_analysis()

        print "Success."
        return

    def cleanup(self):
        """Relinquish any resources acquired"""
        return None



    """input_dir: director containing input"""

    _input_dir = None

    def get_input_dir(self):
        return self._input_dir

    def set_inputdirs(self, val):
        self._input_dir = val

    input_dir = property(get_input_dir, set_inputdirs)

    """output_dir: directory containing output"""

    _output_dir = None

    def get_output_dir(self):
        return self._output_dir

    def set_output_dir(self, val):
        self._output_dir = val

    output_dir = property(get_output_dir, set_output_dir)

    """log file"""
    @property
    def log_file(self):
        return None

    def write_log(self, data):
        return None

    """Analysis metadata"""

    @property
    def start_time(self):
        return None

    @property
    def end_time(self):
        return None

    @property
    def core_allocation(self):
        return 8

    """Miscellaneous"""

    MobileElementInfo = namedtuple('MobileElementInfo',
                               ['sangerID', #id of Sanger read
                                'elementID', #id of mobile element
                                'startIndex', #index the mobile element inserted
                                'validity', #proportion of correctly matched n.t.
                                'mutations']) #list of differences bt. read & mob
    """Declare flags."""
    NO_CALL = 1
    NO_MUT = 2

    """Ancestor sequence (Biopython Sequence object)"""
    _ancestor_seq_obj = None
    @property
    def ancestor_seq_obj(self):
        return self._ancestor_seq_obj

    # Utility functions

    def seqproc(self):
        """Create templates that will be later aligned with MAFFT"""

        # Create list containing mobile elements
        mob_list = list()
        for dirName, subdirList, fileList in os.walk(self.mobile_ele_dir):
            for f in fileList:
                mob_seq = SeqIO.read(dirName+f, "genbank")

                mob_rc = mob_seq.reverse_complement()
                mob_rc.id = mob_seq.id + "-reverse_complement"
                mob_list.append(mob_seq)
                mob_list.append(mob_rc)

        read_list = list()
        plasmid = SeqIO.read(self.ancestor_seq_path, "genbank") # object containing plasmid info
        plasmid_id = plasmid.id

        # Create a list of trimmed samples
        for read_file in self.sample_lst:
            out_file = SeqIO.read(read_file, "abi")
            trimmed_file = self.trim_n(out_file, self.phred_cutoff)  # trim n's
            read_list.append(trimmed_file)

        # At this point we have a list of all Sanger reads corresponding to our current plasmid file.
        # Ends have been trimmed of N's.

        for i, sequence in enumerate(read_list):
            sample_id = sequence.id
            for j, mob_seq in enumerate(mob_list):
                fname = "rm-"+sample_id+"-"+mob_seq.id
                temp_lst = [sequence, mob_seq]
                SeqIO.write(temp_lst, self.template_dir+fname, "fasta")
                del temp_lst

        for i, sequence in enumerate(read_list):
            sample_id = sequence.id
            fname = "pr-"+sample_id+"-"+plasmid_id
            temp_lst = [plasmid, sequence]
            SeqIO.write(temp_lst, self.template_dir+fname, "fasta")
            del temp_lst
        return

    def trim_n(self, seq, phred_cutoff):
        """Trim the 5' and 3' ends of Sanger reads.
        :param seq: sequence object to trim
        :param phred_cutoff: cutoff score for Sanger reads"""

        # Trim 5' end
        i = 0
        skip = False
        done = False

        while (i < len(seq.seq)) and not(done):
            if (((seq.letter_annotations.values()[0][i]+
                  seq.letter_annotations.values()[0][i+1]+
                  seq.letter_annotations.values()[0][i+2])/3) > 30):
                seq = seq[i:]
                done = True
                break
            i += 1

        # Trim 3' end
        j = i
        while (j < len(seq.seq)) and not(done):
            if (((seq.letter_annotations.values()[0][j]+
                  seq.letter_annotations.values()[0][j+1]+
                  seq.letter_annotations.values()[0][j+2])/3) < phred_cutoff):
                seq = seq[i:j]
                break
            j += 1
        return seq

    def create_alignment(self, lst):
        """Create several alignments file_list, and output them to outfile_path.

        This function handles the generation of alignments with MAFFT, and options
        should be changed here. This function runs concurrently on different cores
        to achieve higher throughput."""
        for f in lst:
            with open(self.alignment_dir+f, "w") as out_file:
                check = subprocess.call(["/usr/lib/mafft/bin/mafft", "--quiet",
                                        self.template_dir+f], stdout=out_file)
                sys.stdout.flush()
                if (check):
                    err_str = "MAFFT error on input", f
                    self.write_log(err_str)
        return

    def create_alignment_driver(self):
        """Partition templates and send to separate cores for alignment"""

        lst = list()
        for d, sd, fl in os.walk(self.template_dir):
            pass

        lst = partition_lst(fl, self.core_allocation)

        proc_list = list()

        # Create a new process for each partition of the file list
        for template_list in lst:
            proc = Process(target=self.create_alignment,
                           args=(template_list,))
            proc_list.append(proc)
            proc.start()

        """Join the processes."""
        for proc in proc_list:
            proc.join()

        return

    def mutation_analysis(self):
        """Gathers mobile element homology by stepping through the alignment
        path and analyzing each alignment for mutations. Inside the loop,
        map a key (which is a sample ID) onto a named tuple that contains
        information about the mobile element with which the sample shares
        the highest degree of homology."""
        print "Gathering mobile element homology:"
        mob_evidence_dict = dict()
        for dirName, subdirList, fileList in os.walk(self.alignment_dir):
            """Build necessary file lists for analysis."""
            plasmid_read_filelist = list()
            read_mobile_filelist = list()
            for f in fileList:
                if "pr" in f:
                    plasmid_read_filelist.append(f)
                elif "rm" in f:
                    read_mobile_filelist.append(f)

        self.analyze_mobile_homology(read_mobile_filelist, mob_evidence_dict)

        print "Analyzing alignments of plasmids/reads:"
        mutation_lst = self.parse_samples(plasmid_read_filelist, mob_evidence_dict)

        return mutation_lst

    def analyze_mobile_homology(self, fileList, mob_evidence_dict):
        """Gather data relating to homology between samples and mobile elements.

        This function uses a namedtuple to store information about homology. After calling
        helper function mob_info(), update the mob_evidence_dict (a dictionary) with keys
        (sample names) and values (namedtuple containing homology data)."""
        file_num = len(fileList)
        total_complete = 0
        for f in fileList:
            """Display current progress."""
            percent_complete = Decimal(total_complete)/Decimal(file_num)
            percent_complete = int(percent_complete * 100)
            # print str(percent_complete) + "% complete...\r",
            """Here, we determine mobile element homology to inform our
            analysis later in the work flow."""
            rm_alignments = AlignIO.read(self.alignment_dir+f, "fasta")

            """Generate a named tuple that contains information about
            the strongest (if any) mobile element insertion."""
            mobile_evidence_tuple = EvolutionaryAnalysis.mob_info(rm_alignments)

            if not(mobile_evidence_tuple is None):
                key = mobile_evidence_tuple.sangerID
                value = mobile_evidence_tuple
                mob_evidence_dict[key]= value

            total_complete += 1
        print "100% complete."

        return

    def parse_samples(self, fileList, mob_evidence_dict):
        """Traverse the list of alignments containing plasmids and reads and
        analyze them for mutations.

        This function utilizes a helper function analyze_sample to perform the
        analysis. Then, parse_samples outputs the results to a file named for
        the current sample."""

        """Variables used to display progress to the user"""
        file_num = len(fileList)
        total_complete = 0
        """Build a list of files that contain alignments between plasmids
        and reads."""

        for f in fileList:
            """Display current progress."""
            percent_complete = Decimal(total_complete)/Decimal(file_num)
            percent_complete = int(percent_complete * 100) # Estimate percentage complete
            # pytprint str(percent_complete) + "% complete...\r", # Move cursor to beginning
            """Define names for file I/O"""
            mutation_list = list()

            outfile_name = self.analyze_sample(self.alignment_dir+f, mob_evidence_dict, mutation_list)

            """At this point, we have the list of mutations. Now we write them to disk."""
            if not outfile_name:
                print "ERROR: no outfile name"
                continue

            with open(self.output_dir+outfile_name+ ".gd", "w") as outfile:
                outfile.write("#=GENOME_DIFF 1.0\n")
                for mutation in mutation_list:
                    outfile.write(mutation + "\n")
            total_complete += 1
        return

    def analyze_sample(self, infile_path, mob_evidence_dict, mutation_list):
        """This function performs the actual analysis of samples against their
        respective template plasmids (i.e., this is where we generate most
        of the most important information to our future analyses).

        For each alignment file, we first check if we had previously identified
        homology between the sample and some mobile element. If so, we
        truncate analysis at the start of the insertion.

        This function assigns two of its arguments, outfile_name and mutation_list,
        to new values. The caller expects them to have a certain state. Modify these
        values with care."""

        """Begin analysis of sample and template alignment."""
        pr_alignment = AlignIO.read(infile_path, "fasta") # Read sample file
        outfile_name = pr_alignment[1].id

        """Get indices at which to begin and terminate analysis.
        To determine the start, we use the maximum of the two
        sequences; for the end, we use the minimum of the two.
        This prevents mistakenly identifying mutations in the
        portions of the alignment where the sequences do not
        overlap."""
        plasmid_indices = self.get_indices(pr_alignment[0])
        read_indices = self.get_indices(pr_alignment[1])
        start = max(plasmid_indices[0], read_indices[0])
        stop = min(plasmid_indices[1], read_indices[1])

        """Determine whether or not the sample exhibits substantial
        homology to mobile elements. If so, we can use the start
        of the insertion to truncate analysis, thereby preventing
        erroneous mutations from appearing in our final analysis."""
        key = pr_alignment[1].id
        if (key in mob_evidence_dict):
            strand = 1
            """Compute the index at which to terminate analysis. This
            should be based on two things:
            1) The index (on the sample) of the insertion
            2) The offset at which the sample first appears in this
            (current) alignment file
            Taking these two factors into account yields the index
            at which the mobile element insertion occurs in the current
            alignment.

            However, when documenting any mutations, we are
            interested in where those mutations occur in relation to
            the original template sequence only. To this end, we must
            subtract from the index obtained above the count of insertions
            that occur in our sample prior to said mobile element insertion."""
            stop = start + mob_evidence_dict[key].startIndex
            mob_ele_id = mob_evidence_dict[key].elementID
            self.analyze_mobile_ins(mob_ele_id, pr_alignment[0], pr_alignment[1], start, stop, mutation_list)
        else:
            validity = self.analyze_seq(pr_alignment[0], pr_alignment[1], start, stop, mutation_list)

        return outfile_name

    def analyze_mobile_ins(self, mob_ele_id, template, sample, start, stop, mutation_list):
        """Analyze a sequence, compensating for mobile element insertions.
        This function finalizes parameters for analysis and proceeds with the analyze function.
        This function performs no I/O; it is left to the caller."""
        strand = 1
        if ("reverse_complement" in mob_ele_id):
            temp = re.split("-", mob_ele_id)
            mob_ele_id = temp[0]
            strand = -1
        """Compensate for insertions before mobile element insertion to properly map
        to reference seq"""
        validity = EvolutionaryAnalysis.analyze_seq(template, sample,
                                                    start, stop, mutation_list)
        ins_count = 0
        for mutation in mutation_list:
            tmp = re.split("\t", mutation)
            if ("INS" in tmp[0]):
                ins_count += len(tmp[5]) # add length of any insertions

        mutation_string = ("MOB\t.\t.\t"+template.id+"\t"+str(stop-ins_count) + "\t" +
                           mob_ele_id + "\t" + str(strand))
        mutation_list.append(mutation_string)

        return

    def analyze_seq(self, template, target, start_index, stop_index, mutation_list):
        """Analyze an alignment between template and target
        between the given indices."""
        template_seq = template.seq
        target_seq = target.seq
        template_id = template.id
        if ("reverse_complement" in template_id):
            template_id = (re.split("-", template_id))[0]


        ins_count = 0 # total number of nts inserted
        del_count = 0 # total number of nts deleted
        skip = 0

        true_count = 0
        for x in range(start_index, stop_index):
            # calculate "true index" to map mutations onto reference sequence
            # only things that affect index into reference are start index and total number of inserted bases
            ref_index = x - ins_count
            if (skip):
                skip -= 1
                continue
            template_nt = template.seq[x]
            target_nt = target.seq[x]

            """Evaluate the alignment at index x, identifying any mutations
            and returning them."""
            mutObject = self.eval_nt(template, target, x, ref_index)

            """Did a mutation occur?"""
            if (mutObject == self.NO_MUT):
                true_count += 1
                continue

            """Was there a base call for the target sequence at this location?"""
            if (mutObject == self.NO_CALL):
                continue

            mut_type = mutObject.get_mutation_type()
            """Do we need to skip nucleotides following this mutation?
            We can calculate this for in/dels by subtracting 1 from the
            length of the in/del."""
            if (mut_type == "DEL" or mut_type == "INS"):
                skip = mutObject.get_length() - 1

            """Update del_count and ins_count"""
            if (mut_type == "DEL"):
                del_count += 1
            elif (mut_type == "INS"):
                ins_count += 1

            """Add the Mutation object's string field to the mutation list."""
            mutation_list.append(mutObject.get_string())

            """Delete the object."""
            del(mutObject)

        validity = Decimal(true_count)/Decimal(stop_index-start_index) # percent nts that match in analysis frame
        return validity

    def get_seq_list(self, seq_path):
        """Return a list of Seq objects from GenBank files located in seq_path."""
        seq_list = []
        for dirName, subdirList, fileList in os.walk(seq_path):
            for f in fileList:
                seq_list.append(SeqIO.read(dirName+f, "genbank"))
        return seq_list
    def get_indices(self, seq):
        """Determine start and stop indices for analysis.

        Obtain start and stop indices for analyzing a single alignment
        within a multiple sequence alignment file. This sequence should be
        smaller than your template."""
        start = 0
        stop = len(seq.seq)-1
        while ((seq.seq[stop] == "-") and (stop > 0)):
            stop -= 1
        if (stop == 0):
            print "ERROR: stop_index = 0 for sequence", seq.id
            return
        while (seq.seq[start] == "-") and (start < len(seq.seq)):
            start += 1
        if (start >= len(seq.seq)):
            print "ERROR: start >= length for sequence", seq.id
            return
        if (start > stop):
            print "ERROR: start > stop for sequence", seq.id
            return
        if (stop > len(seq.seq)):
            print "ERROR: stop > len(seq.seq) for sequence", seq.id
        #print "returning", start, stop, "for seq", seq.id
        return [start, stop]
    def get_sub_rep(self, seq_frag):
        """Return the length of the smallest repeated sequence in the
           sequence fragment, and how many times the sequence is repeated."""
        repeat_num = 1
        for x in range(1, len(seq_frag)):
            i = x
            j= 0
            while (i<len(seq_frag) and str(seq_frag[j:i]) == str(seq_frag[i:i+x])):
                repeat_num += 1
                i += x
                j += x
            if (repeat_num * x == len(seq_frag)):
                return [seq_frag[:x], repeat_num]
            else:
                repeat_num = 1

        return [seq_frag, 1]

     # determine number of repeats, if any
    def build_repeat_string(self, seq, x, mut_length, ins, tandem_seq):
        """Determine number of tandem repeats in seq, starting at x."""

        """Starting from seq[x], determine the number of repeated units
        of length mut_length directly downstream or upstream of
        seq[x:x+mut_length].
        arg: ins is boolean (1 -> insertion, 0 -> deletion"""
        mut_str = ""
        seq_str = str.replace("-", "", str(seq))
        sub_rep = self.get_sub_rep(seq_str[x:x+mut_length])
        repeat_seq = sub_rep[0]
        delta_units = sub_rep[1] # number of repeated units deleted/inserted
        # determine if read sequence contains repeats
        # now we must find the bounds of the read repeat sequence
        repeat_len = len(repeat_seq)
        start_index = x
        stop_index = x+repeat_len
        skip = mut_length-1

        """This is an error condition and will cause the process to exit."""
        if (repeat_len <= 0):
            sys.exit("ERROR: repeat_len <= 0")

        """First, we move toward the 5' end by repeat_len each iteration"""
        while (str(seq_str[start_index:start_index+repeat_len])
        == str(seq_str[start_index-repeat_len:start_index]) and start_index>0):
            start_index -= repeat_len

        """Now, we move toward the 3' end in the same manner"""
        while (str(seq_str[stop_index-repeat_len:stop_index])
        == str(seq_str[stop_index:stop_index+repeat_len]) and
               (stop_index+repeat_len)<len(seq_str)):
            stop_index += repeat_len

        repeat_ref_num = (stop_index-start_index)/repeat_len
        if (repeat_ref_num == 1 and delta_units == 1):
            """Not repeat mediated."""
            return
        elif (repeat_len * repeat_ref_num >= 5):
            """Meets GenomeDiff criteria for RMD/RMI"""
            repeat_new_copies = 0
            if (ins):
                repeat_new_copies = repeat_ref_num + delta_units
            else:
                repeat_new_copies = repeat_ref_num - delta_units

            tandem_seq = repeat_seq
            mut_str = ("\trepeat_seq="+ repeat_seq + "\trepeat_len="+str(repeat_len)+
                      "\t"+"repeat_ref_num="+str(repeat_ref_num)+
                      "\trepeat_new_copies="+str(repeat_new_copies))
            return mut_str
        else:
            return


    def identify_disjoint_rmd(self, template_seq, x, del_length):
        """Determine whether or not a disjoint repeat mediated deletion
        occurred at this location in the template.

        There are two ways for the sequences to align if this occurs.
        The alignment program can either place the deleted region first,
        or it can align the complement of the repeat sequence first.
        The ASCII art below gives examples of both cases. Repeat sequences
        are bracketed.

        Case 1: Deletion first
        Template: [ctagag]-----[ctagag]tag
        Sample:    ------ ----- ctagag tag

        Case 2: One copy of repeat sequence first
        Template: [ctagag]-----[ctagag]tag
        Sample:   [ctagag]----- ------ tag

        This function will account for both possibilities: if the deletion
        occurs first, we build the repeat sequence starting at the 5' end
        and compare to the region just downstream of the deletion region,
        building the length of the repeat sequence one at a time up to a
        reasonable limit. If a copy of the repeat sequence is aligned first,
        we build the sequence starting at the 3' end of the deletion region,
        and we compare to the region just upstream of the deletion region."""

        repeat_length = 3
        repeat_seq = None

        """If the deletion length is less than our minimum repeat length,
        this algorithm will fail with false positives, so we terminate if
        the deletion length is too small."""
        if (del_length < 3):
            return
        """Case 1: deletion aligned first. Start building from the 5' end
        (inside the deletion region) and compare to the region just downstream
        of the deletion region."""

        while (repeat_length < 10):
            upstr_index = x
            downstr_index = x + repeat_length
            candidate_seq = template_seq[upstr_index:downstr_index]
            compare_upstr_index = x + del_length
            compare_downstr_index = compare_upstr_index + repeat_length
            compare_seq = template_seq[compare_upstr_index:compare_downstr_index]
            if str(candidate_seq) == str(compare_seq):
                repeat_seq = candidate_seq
                repeat_length += 1
            else:
                break

        if not(repeat_seq is None):
            return repeat_seq


        """Case 2: a copy of the repeat sequence aligned first. Start building
        from the 3' end (inside the deletion region) and compare to the region
        just upstream of the deletion region."""


        """The repeat sequence should be, at minimum, 3 nucleotides in length."""
        repeat_length = 3
        while (repeat_length < 10):
            downstr_index = x+del_length
            upstr_index = downstr_index-repeat_length
            candidate_seq = template_seq[upstr_index:downstr_index]
            compare_upstr_index = x-repeat_length
            compare_downstr_index = x
            compare_seq = template_seq[compare_upstr_index:compare_downstr_index]
            if str(candidate_seq) == str(compare_seq):
                repeat_seq = candidate_seq
            else:
                break

        if not(repeat_seq is None):
            return repeat_seq

        return


    def eval_del(self, mutObject, template, target, x, ref_index):
        """Evaluate a previously identified deletion.

        At entry, the caller has identified that there is a deletion at
        index x in the alignment. Assign a Mutation object with
        information based on what kind of deletion it is."""

        mutObject.set_mutation_type("DEL")

        template_seq = template.seq
        template_id = template.id
        target_seq = target.seq

        """Determine the length of the deletion and store it in the
        Mutation object."""
        z = x
        while (z < len(target_seq) and target_seq[z] == "-"):
            z += 1
        del_length = z - x
        mutObject.set_length(del_length)

        """Build the initial output string."""
        mutation_str = ("DEL\t.\t.\t"+ template_id + "\t" +
                        str(ref_index) + "\t" + str(del_length))

        """Determine if this is a RMD, which can be one of two cases:
        1) tandem repeats
        2) between disjoint repeats
        These cases are mutually exclusive."""

        """Determine if this is an RMD caused by tandem repeats."""
        tandem_seq = None

        """tandem_rmd_str holds a GenomeDiff mutation string if there is
        a tandem RMD here."""
        tandem_rmd_str = self.build_repeat_string(template_seq, x, del_length, 0,
                                            tandem_seq)
        """If tandem repeats, we append the new information."""
        if (tandem_rmd_str):
            mutation_str += tandem_rmd_str
        """Regardless of status, set the mutation string of the Mutation
        object."""
        mutObject.set_string(mutation_str)

        """Determine if this is a disjoint repeat mediated deletion.
        Mutually exclusive with tandem repeats."""
        if not(tandem_rmd_str):
            disj_rmd_seq = self.identify_disjoint_rmd(template_seq, x, del_length)
            if (disj_rmd_seq):
                errata = "\tbetween="+str(disj_rmd_seq)
                mutObject.set_string(mutation_str+errata)
        else:
            mutObject.set_string(mutation_str)

        return 0

    def eval_ins(self, mutObject, template, target, x, ref_index):
        """Evaluate a previously identified insertion.

        At entry, the caller has identified that there is an insertion at
        index x in the alignment. Assign a  Mutation object with
        information based on what kind of insertion it is."""


        mutObject.set_mutation_type("INS")

        template_seq = template.seq
        template_id = template.id
        target_seq = target.seq
        z = x
        ignore_count = 0
        while (z < len(template_seq) and template_seq[z] == "-"):
            z += 1
        ins_length = z - x

        mutation_str = ("INS\t.\t.\t" + template_id + "\t" +
                           str(ref_index) + "\t" +
                           str(target_seq[x:x+ins_length]))
        tandem_seq = None
        errata = None

        """tandem_rmi_str holds a GenomeDiff mutation string if there is
        a tandem RMI here."""
        tandem_rmi_str = self.build_repeat_string(target_seq, x, ins_length, 1,
                                     tandem_seq)
        """If tandem repeats, we append the new information."""
        if (tandem_rmi_str):
            mutation_str += tandem_rmi_str
        """Regardless of status, set the mutation string of the Mutation
        object."""
        mutObject.set_string(mutation_str)

        return 0


    def eval_snp(self, mutObject, template, target, x, ref_index):
        """Evaluate a previously identified SNP.

        At entry, the caller has identified that there is a SNP at
        index x in the alignment. Assign a  Mutation object with
        information based on the SNP's location."""

        mutObject.set_mutation_type("SNP")
        template_id = template.id
        target_nt = target.seq[x]
        mutObject.set_string("SNP\t.\t.\t"+template_id+"\t"+str(ref_index)+"\t"+target_nt)

        return 0


    def eval_nt(self, template, target, x, ref_index):
        """Evaluate whether or not a mutation occurred here,
        and call helpers to identify them.

        This function returns a Mutation object that contains
        information relevant to the identified mutation. All
        information needed by the caller should be returned
        in this object."""

        template_nt = template.seq[x]
        target_nt = target.seq[x]
        """Some kind of discrepancy here."""
        if (target_nt != template_nt):

            """Define Mutation object that will be passed to evaluation
            functions."""
            mutObject = Mutation()
            """Sequencing could not reach a consensus for any nucleotide here,
            so we don't count these as mutations."""
            if (target_nt == "n"):
                del(mutObject)
                return self.NO_CALL
            elif (target_nt == "-" and template_nt != "-"):
                """Deletion"""
                check = self.eval_del(mutObject, template, target, x, ref_index)
                return mutObject
            elif (template_nt == "-" and target_nt != "-"):
                """Insertion"""
                check = self.eval_ins(mutObject, template, target, x, ref_index)
                return mutObject
            else:
                check = self.eval_snp(mutObject, template, target, x, ref_index)
                return mutObject

        else:
            return self.NO_MUT




    def mob_info(self, alignments):
        """Returns a list of info about an alignment analysis based on
        sequence homology to mobile elements."""

        read = alignments[0] # This is the Sanger read

        mob_evidence_list = list()
        for r in range(1, len(alignments)):
            mob = alignments[r] # this is the current mobile element we want to analyze
            # Analyze the alignment only in the region of the mobile element
            mob_indices = self.get_indices(mob)
            read_indices = self.get_indices(read)
            # We want indices that capture only portions of the alignment where both sequences
            # are aligned. This may overlook edge insertions/deletions.
            start = max(mob_indices[0], read_indices[0])
            stop = min(mob_indices[1], read_indices[1])
            if ((stop - start) < 50): # not long enough for any meaningful analysis
                return
            # call analyze_seq on our read/mobile element to get info about sequence homology
            current_element_mutation_list = list()
            validity = self.analyze_seq(mob, read, start, stop, current_element_mutation_list)
            current_mobile_element_info = self.MobileElementInfo(read.id, mob.id, start, validity,
                                                            current_element_mutation_list)
            mob_evidence_list.append(current_mobile_element_info)
            del current_element_mutation_list
        curr_info_element = ""
        curr_len = 9999999 # do not use this value for "big data" alignments/analyses

        """Shouldn't attempt to access items in empty list"""
        if not(mob_evidence_list):
            return
        """The goal is to return the strongest possible evidence of a mobile
        element insertion. To this end, we must do two things:
        1) Identify the collection of information which has the fewest differences
           between its sample and mobile element
        2) Verify that this pair actually exhibits a significant degree of homology"""
        for info_element in mob_evidence_list: # find the mobile element with the least number of mutations
            if (len(info_element.mutations) < curr_len):
                curr_len = len(info_element.mutations)
                curr_info_element = info_element
        if (curr_info_element and curr_info_element.validity >= .45): # is the mobile element actually homologous?
            del(mob_evidence_list)
            return curr_info_element

        return

    """Build a dictionary of information relating to mobile element
    homology, such that we can determine if a given sample exhibits
    substantial homology to any mobile elements (note: only
    informatino about the mobile element which is most similar to
    the sample is stored in the helper function mob_info).

    The dictionary that this function returns can be indexed by
    sample ID, i.e. each sample ID is a key in the dictinoary,
    which stores a MobileElementInfo namedtuple."""
    def build_mob_evidence_dict(self):
        mob_evidence_dict = dict()
        for dirName, subdirList, fileList in os.walk(self.alignment_dir):
            for f in fileList:
                if "rm" in f:
                    # read/mobile element alignment
                    # need to determine if any homology to do nt analysis
                    rm_alignments = AlignIO.read(self.alignment_dir+f, "fasta")
                    print ("Gathering mobile element homology data for read " +
                    rm_alignments[0].id + " and " + rm_alignments[1].id + "..."),
                    mobile_element_evidence = self.mob_info(rm_alignments)
                    if not(mobile_element_evidence is None):
                        key = mobile_element_evidence['sangerID']
                        value = mobile_element_evidence
                        mob_evidence_dict[key]= value
                    print "done."
        return mob_evidence_dict
