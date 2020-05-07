#
# shift_reads.py
#
# This script shifts aligned ATAC-seq reads for transcription factor footprinting.
# Reads aligning to the positive strand are shifted by +4 nt; reads aligning to the
# negative strand are shifted by -5 nt.
#
# This is required for the ATAC-seq protocol in Buenrostro et al. (2013) Nat Methods
# 10(12):1213-1218.

import os.path
import argparse
import re
import pysam

parser = argparse.ArgumentParser(description = "This script shifts aligned ATAC-seq reads for transcription factor footprinting. Reads aligning to the positive strand are shifted by +4 nt; reads aligning to the negative strand are shifted by -5 nt. This is required for the ATAC-seq protocol in Buenrostro et al. (2013) Nat Methods 10(12):1213-1218.")
parser.add_argument("infile", metavar = "<input file>", help = "This is assumed to be in BAM format, and must end with a .bam extension.")
parser.add_argument("-o", dest = "outfile", metavar = "<output file>", required = False, help = "If this is not specified, output will be written to a file where the .bam extension of <input file> is replaced with .shifted.bam")
parser.add_argument("-mf", dest = "maxfragment", metavar = "[max frag size]", type = int, required = False, help = "Only include fragments <= a maximum size, e.g. sub-nucleosomal fragments <= 100 nt per Mo et al. (2015) Neuron 86(6):1369-1384")

args = parser.parse_args()

infile = ""
outfile = ""

if (not re.search(r"\.bam$", args.infile, re.IGNORECASE)):
    # If input file does not end with .bam (case insensitive), then exit with error
    parser.print_help()
    parser.exit(status = 1, message = "\nERROR: Input file does not end with .bam\n")
elif (not args.outfile):
    # Construct output file name from input file name
    regex = re.compile(r"\.bam$", re.IGNORECASE)
    infile = args.infile
    outfile = regex.sub(r".shifted.bam", args.infile)
else:
    # User has specified both input and output file names, go with those
    infile = args.infile
    outfile = args.outfile

# Check to see that input file exists before proceeding with processing
if(os.path.isfile(infile)):
    in_file = pysam.AlignmentFile(infile, "rb")
    out_file = pysam.AlignmentFile(outfile, "wb", template=in_file)

    print "Processing %s." % infile
    
    if(args.maxfragment):
        # Exclude fragments above a specified size maxfragment
        # Loop over all reads in input file
        for read in in_file.fetch(until_eof = True):
            if(abs(read.template_length) <= args.maxfragment):
                # The fragment is <= maximum size specified by user
                if read.is_reverse:
                    # Shift read by -5 bp if on negative strand
                    read.reference_start -= 5
                else:
                    # Shift read by +4 bp if on positive strand
                    read.reference_start += 4
                out_file.write(read)
  
    else:
        # Include all fragments
        # Loop over all reads in input file
        for read in in_file.fetch(until_eof = True):
            if read.is_reverse:
                # Shift read by -5 bp if on negative strand
                read.reference_start -= 5
            else:
                # Shift read by +4 bp if on positive strand
                read.reference_start += 4
            out_file.write(read)

else:
    # Input file does not exist
    parser.print_help()
    parser.exit(status = 1, message = "\nERROR: Input file %s does not exist\n" % infile)

print "Writing output to %s." % outfile
in_file.close
out_file.close
