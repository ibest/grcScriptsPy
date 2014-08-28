#!/usr/bin/env python

# This file is part of grcScriptsPy, http://github.com/ibest/grcScriptsPy/
# Copyright 2014, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License; see LICENSE.txt.
# Contact: msettles@uidaho.edu


import argparse 
import sys, os, traceback
import time

from grcScriptsPy import (TwoReadIlluminaRun,OneReadIlluminaRun)

from khmer import ()


DEFAULT_DESIRED_COVERAGE = 30

DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MIN_TABLESIZE = 1e6
DEFAULT_N_THREADS = 1

MAX_FALSE_POSITIVE_RATE = 0.8             # see Zhang et al.,
# http://arxiv.org/abs/1309.2975

profile = False

def batchwise(coll, size):
    iter_coll = iter(coll)
    return izip(*[iter_coll] * size)

# Returns true if the pair of records are properly pairs

def validpair(read0, read1):
    return read0.name[-1] == "1" and \
        read1.name[-1] == "2" and \
        read0.name[0:-1] == read1.name[0:-1]

def handle_error(error, , input_name, output_name):
    print >> sys.stderr, '** ERROR:', error
    print >> sys.stderr, '** Failed on %s: ' % input_name
    try:
        os.remove(output_name)
    except:  # pylint: disable=bare-except
        print >> sys.stderr, '** ERROR: problem removing corrupt filtered file'



# pylint: disable=too-many-locals,too-many-branches
def normalize_by_median(input_filename, outfp, htable, args, report_fp=None):

    desired_coverage = args.cutoff
    ksize = htable.ksize()

    # In paired mode we read two records at a time
    batch_size = 1
    if args.paired:
        batch_size = 2

    index = -1
    total = 0
    discarded = 0
    for index, batch in enumerate(batchwise(screed.open(
            input_filename), batch_size)):
        if index > 0 and index % 100000 == 0:
            print '... kept {kept} of {total} or {perc:2}%'.format(
                kept=total - discarded, total=total,
                perc=int(100. - discarded / float(total) * 100.))
            print '... in file', input_filename

            if report_fp:
                print >> report_fp, total, total - discarded, \
                    1. - (discarded / float(total))
                report_fp.flush()

        total += batch_size

        # If in paired mode, check that the reads are properly interleaved
        if args.paired:
            if not validpair(batch[0], batch[1]):
                raise IOError('Error: Improperly interleaved pairs \
                    {b0} {b1}'.format(b0=batch[0].name, b1=batch[1].name))

        # Emit the batch of reads if any read passes the filter
        # and all reads are longer than K
        passed_filter = False
        passed_length = True
        for record in batch:
            if len(record.sequence) < ksize:
                passed_length = False
                continue

            seq = record.sequence.replace('N', 'A')
            med, _, _ = htable.get_median_count(seq)

            if med < desired_coverage:
                htable.consume(seq)
                passed_filter = True

        # Emit records if any passed
        if passed_length and passed_filter:
            for record in batch:
                if hasattr(record, 'accuracy'):
                    outfp.write(
                        '@{name}\n{seq}\n'
                        '+\n{acc}\n'.format(name=record.name,
                                            seq=record.sequence,
                                            acc=record.accuracy))
                else:
                    outfp.write(
                        '>{name}\n{seq}\n'.format(name=record.name,
                                                  seq=record.sequence))
        else:
            discarded += batch_size

    return total, discarded


class khmer_normalize_app:
    """
    application to normalize reads using khmer
    """
    def __init__(fastq_pair1, fastq_pair2, fastq_single, cutoff, ksize, n_tables, min_tablesize, verbose, debug):
        self.verbose=verbose
        self.debug = debug

        self.cutoff = cutoff
        self.ksize= ksize
        self.n_table = n_tables
        self.min_tablesize = min_tablesize

        try:
            ## establish and open the Illumin run
            if fastq_pair1 != None and fastq_pair2 != None:
                self.runPairs = TwoReadIlluminaRun(fastq_pair1,fastq_pair2)
            else:
                self.runPairs = None
            if fastq_single != None:
                self.runSingle = OneReadIlluminaRun(fastq_single)
            else:
                self.runSingle = None
            if self.runPairs == None and self.runSingle == None:
                print("ERROR: input reads not specified, or incorrect pairs")
                raise Exception


        except (KeyboardInterrupt, SystemExit):
            self.clean()
            print("%s unexpectedly terminated" % (__name__))
            return 1
        except:
            self.clean()
            print("A fatal error was encountered.")
            if debug:
                print "".join(traceback.format_exception(*sys.exc_info()))
            return 1

    def start(self, fastq_file1, fastq_file2, fastq_fileU, output_prefix, rdpPath='./classifier.jar', gene='16srrna', batchsize=10000, procs = 1, verbose=True, debug=False):
        """
            Start the normalization process
        """

        check_space(args.input_filenames)

        print 'making k-mer counting table'
        htable = khmer.new_counting_hash(args.ksize, args.min_tablesize, args.n_tables)

        total = 0
        discarded = 0

        for index, input_filename in enumerate(args.input_filenames):
            if args.single_output_filename != '':
                output_name = args.single_output_filename
                outfp = open(args.single_output_filename, 'a')
            else:
                output_name = os.path.basename(input_filename) + '.keep'
                outfp = open(output_name, 'w')

            total_acc = 0
            discarded_acc = 0

            try:
                total_acc, discarded_acc = normalize_by_median(input_filename,
                                                               outfp, htable, args,
                                                               report_fp)
            except IOError as err:
                handle_error(err, input_filename, output_name)
                print >> sys.stderr, '** Exiting!'
                sys.exit(1)
            else:
                if total_acc == 0 and discarded_acc == 0:
                    print 'SKIPPED empty file', input_filename
                else:
                    total += total_acc
                    discarded += discarded_acc
                    print 'DONE with %s; kept %s of %s or %s%' % (input_filename, total - discarded, total, round(int(100. - discarded / float(total) * 100.),2))
                    print 'output in', output_name

        fp_rate = khmer.calc_expected_collisions(htable)
        print 'fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate)

        if fp_rate > MAX_FALSE_POSITIVE_RATE:
            print >> sys.stderr, "**"
            print >> sys.stderr, ("** ERROR: the k-mer counting table is too small"
                                  " for this data set.  Increase tablesize/# "
                                  "tables.")
            print >> sys.stderr, "**"
            print >> sys.stderr, "** Do not use these results!!"
            sys.exit(1)


    def clean(self):
        if self.verbose:
            print("Cleaning up.")
        try:
            self.runSingle.close()
            self.runPairs.close()
        except:
            pass



##
#####################################################################################
##  Master parser arguments
def parseArgs():
    """
    generate main parser
    """
    parser = argparse.ArgumentParser( \
        description = 'normalize_khmer, using khmer normalize coverage in reads', \
        epilog ='For questions or comments, please contact Matt Settles <msettles@uidaho.edu>', add_help=True)
    parser.add_argument('--version', action='version', version="%(prog)s Version " + grcScriptsPy.__version__)
    parser.add_argument('-C', '--cutoff', help='normalized median cutoff value [default: %(default)s]',
                        type=int, dest='cutoff', default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('--ksize', '-k', help='k-mer size to use [default: %(default)s]',
                        type=int, dest='ksize', default=DEFAULT_K)
    parser.add_argument('--n_tables', '-N', help='number of k-mer counting tables to use [default: %(default)s]',
                        type=int, dest='n_tables', default=DEFAULT_N_TABLES)
    parser.add_argument('--min_tablesize', '-x', help='lower bound on tablesize to use [default: %(default)s]',
                        type=float, dest='min_tablesize', default=DEFAULT_MIN_TABLESIZE)
    parser.add_argument('-1', metavar="read1", dest='fastq_pair1', help='read1 of a fastq paired read set',
                        action='store',type=str, default=None, required=False, nargs='+')
    parser.add_argument('-2', metavar="read2", dest='fastq_pair2', help='read2 of a fastq paired read set',
                        action='store',type=str, default=None, required=False, nargs='+')
    parser.add_argument('-U', metavar="single", dest='fastq_single', help='single-end amplicon, typically from joined paired reads',
                        action='store',type=str, default=None, required=False, nargs='+')
    parser.add_argument('-v', '--silent', help='verbose output [default: %(default)s]',
                        action='store_true', dest='verbose', default=False)
    parser.add_argument('--debug', help='show traceback on error [default: %(default)s]',
                        action='store_true', dest="debug", default = False)

    args = parser.parse_args() 

    return args


def main():
    """
    main function
    """
    args = parseArgs()

    # ----------------------- other options ------------
    debug = args.debug
    verbose = not args.verbose

#    app = classifyApp(args.fastq_pair1, args.fastq_pair2, args.fastq_single, args.cutoff, args.ksize, args.n_tables, args.min_tablesize, args.verbose, args.debug)

#    if profile:
#        import cProfile
#       cProfile.runctx('app.start()', globals(), locals())
#        return 255
#    else:
#        return app.start()


    print("C Version")
    grcScriptsPy.hello_world()
    print("Py Version")
    grcScriptsPy.print_hello_world_py()
    
    print(grcScriptsPy.hamming_distance("ABCDEFG","ABDDEFG"))
    print(grcScriptsPy.hamming_distance("ABCDEFG","ABDDEGG"))


if __name__ == '__main__':
    main()
 
