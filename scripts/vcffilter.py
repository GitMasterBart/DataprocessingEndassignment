##! usr/bin/env python3
"""
For questions reach to:
b.engels@st.hanze.nl
"""
import argparse
import sys
import re

___author__ = 'Bart Engels '
___data__ = '2021-09-20'
___version__ = 'v1-d2020'


class VcfFilter:
    """
    This program is used to filter a vcf file. It has 3 methods for filtering.
    on the bases of quality (QUAL), depth (INFO(DP=)) and decomposing.
    More information about each methode is supplied in the docstring of the methode is question.
    """

    def __init__(self, methode, vcf_file, quality, depth, output):
        self.methode = methode
        self.vcf_file = open(vcf_file)
        self.quality = quality
        self.depth = depth
        self.output = output

    def quality_filter(self):
        """
        The quality filter removes instances with a Quality les than the given param
        so if a instance has a quality of 25 and the param is 30 it write this instance to a file.
        :param:
            file : Line of files
            quality : max Quality score
        :return: lines filtered on Quality
        """
        message = ""
        with open(self.output, 'w') as new_file:
            for line in self.vcf_file:
                if line.startswith("#"):
                    new_file.write(line)
                try:
                    if not line.startswith('#'):
                        if int(line.split('\t')[5]) <= int(self.quality):
                            new_file.write(line)
                except ValueError:
                    return 0
                except TypeError:
                    message += "fill in the -q"
            if message != "":
                print(message)
            new_file.close()

    def depth_filter(self):
        """
         The depth filter removes instances with a DP les than the given param
        so if a instance has a DP of 13 and the param is 10 it write this instance to a file.
        :param:
            file : Line of files
            Depth : Depth of filter
        :return: lines filtered on depth
        """
        with open(self.output, 'w') as new_file:
            for line in self.vcf_file:
                if line.startswith("#"):
                    new_file.write(line)

                if not line.startswith('#'):
                    match = re.search("(DP=)(\d+)", str(line.split('\t')[7]))
                    if int(match.group(2)) > int(self.depth):
                        new_file.write(line)
            new_file.close()

    def decomposer_filter(self):
        """
        The decomposer filer makes a new instance for each item that is
        flagged under this param Alt. If the alt param contains 3 items, say:
        [a,agg, aggg] it wil make a new instace for a, agg and aggg.
        :param: vcf file lines
        :return: Decomposed Vcf file lines
        """
        with open(self.output, 'w') as new_file:
            for line in self.vcf_file:
                if line.startswith("#"):
                    new_file.write(line)

                else:
                    list_container = line.split('\t')
                    if len(list_container[4].split(',')) > 1:
                        for i in list_container[4].split(','):
                            list_container[4] = i
                            new_file.write('\t'.join([str(new_string)
                                                      for new_string in list_container]))
                    else:
                        new_file.write(line)
            new_file.close()

        return 0


def main(args):
    """
    main
    """
    inp_args = argparse.ArgumentParser()

    inp_args.add_argument('methode', help='Fill in the methode that you want to use.')
    inp_args.add_argument('-i', '--infile', help='Supply the name of the input-file, '
                                                 'this can be combined with the filepath')
    inp_args.add_argument('-q', '--quality', help='Supply the maximum border of withs the quality is allowed to be.')
    inp_args.add_argument('-d', '--depth', help='Give the minimum depth that is allowed for '
                                                'the instance in the vcf-file')
    inp_args.add_argument('-o', '--outfile', help='Supply the name of the output-file, '
                                                  'this can be combined with the filepath')

    args = inp_args.parse_args()
    methode = args.methode
    vcf_file = args.infile
    quality = args.quality
    depth = args.depth
    output = args.outfile

    run = VcfFilter(methode, vcf_file, quality, depth, output)

    if methode == "Decompose":
        run.decomposer_filter()
    if methode == "Quality":
        run.quality_filter()
    if methode == "Depth":
        run.depth_filter()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
