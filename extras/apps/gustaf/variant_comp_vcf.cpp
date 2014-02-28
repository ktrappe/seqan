// ==========================================================================
//                                 variant_comp_vcf
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "msplazer.h"
#include "msplazer_out.h"
#include <seqan/arg_parse.h>

using namespace seqan;

// This struct stores the options from the command line.
struct VariantCompOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Path to VCF input files
    String<CharString> inPaths;

    // Paths to output files
    CharString stats_out;
    CharString tp_out;
    CharString fp_out;
    CharString fn_out;

    VariantCompOptions() :
        verbosity(1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(VariantCompOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("variant_comp_vcf");
    // Set short description, version, and date.
    setShortDescription(parser, "Comparing two variant files in VCF format.");
    setVersion(parser, "0.1");
    setDate(parser, "March 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIREF VARIANTS VCF FILE\\fP\" \"\\fIPREDICTED VARIANTS VCF FILE\\fP\"");
    addDescription(parser, "Comparing a VCF file with predicted variants to a VCF file with known variants.");
    addDescription(parser, "This simple program takes as input two VCF files with structural variants "
                            "where the first file is supposed two contain known variants and the second file "
                            "to contain predicted variants that should be compared to the known variants.");

    /*
    addDescription(parser, "To prepare the joined mate file for the Gustaf paired-end example, call \n "

    "./gustaf_mate_joining adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa "
        "-rc -o adeno_modified_reads_joinedMates.fa");
    */

    // We require two arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "VCF FILE(S)", true));
    setValidValues(parser, 0, "vcf");
    /*
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA/FASTQ FILE 2"));
    setValidValues(parser, 1, "fasta fa fastq fq");
    */

    addOption(parser, ArgParseOption("o", "outPath", "Set name of output statistic file(s).", ArgParseOption::OUTPUTFILE, "txt"));
    setValidValues(parser, "o", "txt");
    setDefaultValue(parser, "o", "variant_comp.txt");

    addOption(parser, ArgParseOption("tp", "tp_outPath", "Set name of output true positive file.", ArgParseOption::OUTPUTFILE, "vcf"));
    setValidValues(parser, "tp", "vcf");
    setDefaultValue(parser, "tp", "tp_out.vcf");
    addOption(parser, ArgParseOption("fp", "fp_outPath", "Set name of output false positive file.", ArgParseOption::OUTPUTFILE, "vcf"));
    setValidValues(parser, "fp", "vcf");
    setDefaultValue(parser, "fp", "fp_out.vcf");
    addOption(parser, ArgParseOption("fn", "fn_outPath", "Set name of output false negative file.", ArgParseOption::OUTPUTFILE, "vcf"));
    setValidValues(parser, "fn", "vcf");
    setDefaultValue(parser, "fn", "fn_out.vcf");

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    resize(options.inPaths, getArgumentValueCount(parser, 0), Exact());
    for (unsigned i = 0; i < length(options.inPaths); ++i)
        getArgumentValue(options.inPaths[i], parser, 0, i);
    if (length(options.inPaths) < 2)
        return ArgumentParser::PARSE_ERROR;
    getOptionValue(options.stats_out, parser, "o");
    getOptionValue(options.tp_out, parser, "tp");
    getOptionValue(options.fp_out, parser, "fp");
    getOptionValue(options.fn_out, parser, "fn");

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function _importVariants()
// --------------------------------------------------------------------------

// Import variants from a vcf file and write them into breakpoint format
template <typename TSequence, typename TId>
inline bool
_importVariants(CharString const & fileName,
                String<Breakpoint<TSequence, TId> > & variants
                 )
{

    typedef typename Position<TSequence>::Type TPos;
    VcfStream f(toCString(fileName));
    VcfStream o("-", VcfStream::WRITE);
    if (!isGood(f))
    {
        std::cerr << "Failed to open file." << std::endl;
        return false;
    }

    o.header = f.header;

    VcfRecord record;
    while (!atEnd(f))
    {
        readRecord(record, f);
        writeRecord(o, record);
    }

    return true;
}

// --------------------------------------------------------------------------
// Function _writeVariants()
// --------------------------------------------------------------------------

// Writes out variants in vcf format
template <typename TSequence>
int _writeVariants(CharString & outPath
                    )
{
    /*
    SequenceStream seqStream(toCString(outPath), SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "Error: Could not open output file!" << std::endl;
        return 1;
    }
    if (writeAll(seqStream, sIds, seqs, quals) != 0)
    {
        std::cerr << "Error: Could not write to file!" << std::endl;
        return 1;
    }
    */
    return 0;
}

// Writes out sequences and ids in FASTA format
template <typename TSequence>
int _writeStats(CharString & outPath,
                StringSet<TSequence> const & seqs
                )
{
    return 0;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    VariantCompOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
    {
        std::cout << "Error parsing command line, please check correct number and values of input parameters!" << std::endl;
        return res == ArgumentParser::PARSE_ERROR;
    }

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY     \t" << options.verbosity << '\n'
                  << "INPUT FILE 1     \t" << options.inPaths[0] << '\n'
                  << "INPUT FILE 2     \t" << options.inPaths[1] << '\n';
                  std::cout << "OUTPUT FILE     \t" << options.stats_out << '\n';
                  // Make output to files optional
                      std::cout << "TP OUTPUT FILE     \t" << options.tp_out << '\n';
                      std::cout << "TP OUTPUT FILE     \t" << options.fp_out << '\n';
                      std::cout << "TP OUTPUT FILE     \t" << options.fn_out << '\n';
    }

    typedef String<Dna5> TSequence;
    typedef CharString TId;
    String<Breakpoint<TSequence, TId> > ref_variants;

    // Read in reference variants
    _importVariants(options.inPaths[0], ref_variants);

    return 0;
}
