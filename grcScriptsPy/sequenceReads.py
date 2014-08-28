#!/usr/bin/env python

# Copyright 2013, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
sequenceReads.py stores and processes individual DNA sequence reads.
"""

from grcScriptsPy import misc


# ---------------- Class for 2 read sequence data processed with dbcAmplicons preprocess ----------------
class TwoSequenceReadSet:
    """ 
    Class to hold one Illumina two read set, assumed to have already been preprocessed with Barcodes and Primers already
    identified. Class processes a read by defining sample and project ids. Finally class returns a paired read set for output.
    """
    def __init__(self,name_1,read_1,qual_1,name_2,read_2,qual_2):
        """
        Initialize a TwoSequenceReadSet with names, two read sequences and cooresponding quality sequence.
        Barcode and primer sequence is inferred by their placement in the read names.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.barcode = split_name[1].split(":")[3]
            self.sample = self.barcode
            if (len(split_name) == 4):
                self.primer_string1 = split_name[3]
                self.primer_string2 = name_2.split(" ")[3]
                self.primer = split_name[1].split(":")[4]
                self.barcode_string = split_name[2]
            elif (len(split_name) == 3):
                self.primer_string1 = None
                self.primer_string2 = None
                self.primer = None
                self.barcode_string = split_name[2]
            else:
                self.primer_string1 = None
                self.primer_string2 = None
                self.primer = None
                self.barcode_string = None                
            self.read_1 = read_1
            self.qual_1 = qual_1
            self.read_2 = read_2
            self.qual_2 = qual_2
            self.goodRead = False
        except IndexError:
            print 'ERROR:[TwoSequenceReadSet] Read names are not formatted in the expected manner'
            raise
        except:
            print 'ERROR:[TwoSequenceReadSet] Unknown error occured initiating read'
            raise
    def assignRead(self, sTable):
        """
        Given a samplesTable object, assign a sample ID and project ID using the reads barcode and primer designation
        """
        self.project = sTable.getProjectID(self.barcode,self.primer)
        self.sample = sTable.getSampleID(self.barcode,self.primer)
        self.goodRead = False
        if self.project != None:
            self.goodRead = True
        return 0
    def getFastq(self):
        """ 
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        if self.primer != None:
            read1_name = "%s 1:N:0:%s:%s %s %s" % (self.name, self.sample, self.primer, self.barcode_string, self.primer_string1)
            read2_name = "%s 2:N:0:%s:%s %s %s" % (self.name, self.sample, self.primer, self.barcode_string, self.primer_string2)
        else:
            read1_name = "%s 1:N:0:%s %s" % (self.name, self.sample, self.barcode_string)
            read2_name = "%s 2:N:0:%s %s" % (self.name, self.sample, self.barcode_string)
        r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
        r2 = '\n'.join([read2_name, self.read_2,'+',self.qual_2])
        return [r1,r2]
    def getFourReads(self):
        """ 
        Create four line string ('\n' separator included) for the read set, returning a length 4 vector (one for each read)
        """
        try:
            if self.barcode_string is not None:
                bc1 = self.barcode_string.split('|')[0]
                bc2 = self.barcode_string.split('|')[2]
            elif len(self.barcode) is 16: ## assume 8 bp barcodes for now
                bc1 = self.barcode[0:8]
                bc2 = self.barcode[8:16]
            else:
                raise Exception("string in the barcode is not 16 characters")
            r1 = '\n'.join([self.name + ' 1:N:0:', self.read_1, '+',self.qual_1])
            r2 = '\n'.join([self.name + ' 2:N:0:', bc1, '+', 'C' * len(bc1) ]) ## Give barcodes and arbitary quality of Cs
            r3 = '\n'.join([self.name + ' 3:N:0:', bc2, '+', 'C' * len(bc2) ])
            r4 = '\n'.join([self.name + ' 4:N:0:', self.read_2, '+', self.qual_2])
            return [r1,r2,r3,r4]
        except IndexError:
            print 'ERROR:[TwoSequenceReadSet] unable to exract barocode sequence from the read names'
            raise
        except:
            print 'ERROR:[TwoSequenceReadSet] Unknown error occured generating four read set'
            raise
    def getFasta(self):
        """ 
        Create two line string ('\n' separator included) for the read pair, returning a length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        if self.primer != None:
            read1_name = "%s 1:N:0:%s:%s" % (name, self.sample, self.primer)
            read2_name = "%s 2:N:0:%s:%s" % (name, self.sample, self.primer)
        else:
            read1_name = "%s 1:N:0:%s" % (name, self.sample)
            read2_name = "%s 1:N:0:%s" % (name, self.sample)
        r1 = '\n'.join([read1_name, self.read_1])
        r2 = '\n'.join([read2_name, self.read_2])
        return [r1,r2]
    def getJoinedFasta(self):
        """ 
        Create two line string ('\n' separator included) for the read pair, concatenating the two reads into a single returning length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        if self.primer != None:
            read1_name = "%s|%s:%s" % (name, self.sample, self.primer)
        else:
            read1_name = "%s:%s" % (name, self.sample)
        r1 = '\n'.join([read1_name, self.read_1 + misc.reverseComplement(self.read_2)])
        return [r1]
    

# ---------------- Class for 2 read sequence data processed with dbcAmplicons preprocess ----------------
class OneSequenceReadSet:
    """ 
    Class to hold a one Illumina read set, assumes the paired reads produced by dbcAmplicons preprocess have been merged
    """
    def __init__(self,name_1,read_1,qual_1):
        """
        Initialize a OneSequenceReadSet with name, one read sequences and cooresponding quality sequence.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        self.goodRead = False
        self.primer = None
        ## parse the read name for the sample and primer ids
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.sample = split_name[1].split(":")[3]
            if (len(split_name) == 4):
                self.primer = split_name[1].split(":")[4]
        except IndexError:
            print 'ERROR:[OneSequenceReadSet] Read names are not formatted in the expected manner'
            raise
        except:
            print 'ERROR:[OneSequenceReadSet] Unknown error occured initiating read'
            raise            
        self.read_1 = read_1
        self.qual_1 = qual_1
    def getFastq(self):
        """ 
        Create four line string ('\n' separator included) for the read, returning a length 1 vector (one read)
        """
        if self.primer != None:
            read1_name = "%s 1:N:0:%s:%s" % (self.name, self.sample, self.primer)
        else:
            read1_name = "%s 1:N:0:%s" % (self.name, self.sample)
        r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
        return [r1]
    def getFasta(self):
        """ 
        Create two line string ('\n' separator included) for the read, returning a length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        if self.primer != None:
            read1_name = "%s|%s:%s" % (name, self.sample, self.primer)
        else:
            read1_name = "%s|%s" % (name, self.sample)
        r1 = '\n'.join([read1_name, self.read_1])
        return [r1]



