#!/usr/bin/env python
"""Execute the tests for the razers2 program.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os.path
import sys

# Automagically add util/py_lib to PYTHONPATH environment variable.
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                    '..', '..', 'util', 'py_lib'))
sys.path.insert(0, path)

import seqan.app_tests as app_tests

def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for razers2'
    print '==========================='
    print

    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'apps/razers2/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'apps/razers2', 'razers2')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # ============================================================
    # Run Adeno Single-End Tests
    # ============================================================

    # We run the following for all read lengths we have reads for.
    for rl in [36, 100]:
        # Run with default options.
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile('se-adeno-reads%d_1.stdout' % rl),
            args=['--low-memory',
                  ph.inFile('adeno-genome.fa'),
                  ph.inFile('adeno-reads%d_1.fa' % rl),
                  '-o', ph.outFile('se-adeno-reads%d_1.razers' % rl)],
            to_diff=[(ph.inFile('se-adeno-reads%d_1.razers' % rl),
                      ph.outFile('se-adeno-reads%d_1.razers' % rl)),
                     (ph.inFile('se-adeno-reads%d_1.stdout' % rl),
                      ph.outFile('se-adeno-reads%d_1.stdout' % rl))])
        conf_list.append(conf)

        # Allow indels.
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile('se-adeno-reads%d_1-id.stdout' % rl),
            args=['--low-memory', '-id',
                  ph.inFile('adeno-genome.fa'),
                  ph.inFile('adeno-reads%d_1.fa' % rl),
                  '-o', ph.outFile('se-adeno-reads%d_1-id.razers' % rl)],
            to_diff=[(ph.inFile('se-adeno-reads%d_1-id.razers' % rl),
                      ph.outFile('se-adeno-reads%d_1-id.razers' % rl)),
                     (ph.inFile('se-adeno-reads%d_1-id.stdout' % rl),
                      ph.outFile('se-adeno-reads%d_1-id.stdout' % rl))])
        conf_list.append(conf)

        # Compute forward/reverse matches only.
        for o in ['-r', '-f']:
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('se-adeno-reads%d_1-id%s.stdout' % (rl, o)),
                args=['--low-memory', '-id', o,
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      '-o', ph.outFile('se-adeno-reads%d_1-id%s.razers' % (rl, o))],
                to_diff=[(ph.inFile('se-adeno-reads%d_1-id%s.razers' % (rl, o)),
                          ph.outFile('se-adeno-reads%d_1-id%s.razers' % (rl, o))),
                         (ph.inFile('se-adeno-reads%d_1-id%s.stdout' % (rl, o)),
                          ph.outFile('se-adeno-reads%d_1-id%s.stdout' % (rl, o)))])
            conf_list.append(conf)

        # Compute with different identity rates.
        for i in range(90, 101):
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('se-adeno-reads%d_1-id-i%d.stdout' % (rl, i)),
                args=['--low-memory', '-id', '-i', str(i),
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      '-o', ph.outFile('se-adeno-reads%d_1-id-i%d.razers' % (rl, i))],
                to_diff=[(ph.inFile('se-adeno-reads%d_1-id-i%d.razers' % (rl, i)),
                          ph.outFile('se-adeno-reads%d_1-id-i%d.razers' % (rl, i))),
                         (ph.inFile('se-adeno-reads%d_1-id-i%d.stdout' % (rl, i)),
                          ph.outFile('se-adeno-reads%d_1-id-i%d.stdout' % (rl, i)))])
            conf_list.append(conf)

        # Compute with different output formats.
        for suffix in ['razers', 'fa', 'eland', 'gff', 'sam', 'afg']:
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('se-adeno-reads%d_1-id.%s.stdout' % (rl, suffix)),
                args=['--low-memory', '-id',
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      '-o', ph.outFile('se-adeno-reads%d_1-id.%s' % (rl, suffix))],
                to_diff=[(ph.inFile('se-adeno-reads%d_1-id.%s' % (rl, suffix)),
                          ph.outFile('se-adeno-reads%d_1-id.%s' % (rl, suffix))),
                         (ph.inFile('se-adeno-reads%d_1-id.%s.stdout' % (rl, suffix)),
                          ph.outFile('se-adeno-reads%d_1-id.%s.stdout' % (rl, suffix)))])
            conf_list.append(conf)

        # Compute with different sort orders.
        for so in [0, 1]:
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('se-adeno-reads%d_1-id-so%d.stdout' % (rl, so)),
                args=['--low-memory', '-id', '-so', str(so),
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      '-o', ph.outFile('se-adeno-reads%d_1-id-so%d.razers' % (rl, so))],
                to_diff=[(ph.inFile('se-adeno-reads%d_1-id-so%d.razers' % (rl, so)),
                          ph.outFile('se-adeno-reads%d_1-id-so%d.razers' % (rl, so))),
                         (ph.inFile('se-adeno-reads%d_1-id-so%d.stdout' % (rl, so)),
                          ph.outFile('se-adeno-reads%d_1-id-so%d.stdout' % (rl, so)))])
            conf_list.append(conf)

    # ============================================================
    # Run Adeno Paired-End Tests
    # ============================================================

    # We run the following for all read lengths we have reads for.
    for rl in [36, 100]:
        # Run with default options.
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile('pe-adeno-reads%d_2.stdout' % rl),
            args=['--low-memory',
                  ph.inFile('adeno-genome.fa'),
                  ph.inFile('adeno-reads%d_1.fa' % rl),
                  ph.inFile('adeno-reads%d_2.fa' % rl),
                  '-o', ph.outFile('pe-adeno-reads%d_2.razers' % rl)],
            to_diff=[(ph.inFile('pe-adeno-reads%d_2.razers' % rl),
                      ph.outFile('pe-adeno-reads%d_2.razers' % rl)),
                     (ph.inFile('pe-adeno-reads%d_2.stdout' % rl),
                      ph.outFile('pe-adeno-reads%d_2.stdout' % rl))])
        conf_list.append(conf)

        # Allow indels.
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile('pe-adeno-reads%d_2-id.stdout' % rl),
            args=['--low-memory', '-id',
                  ph.inFile('adeno-genome.fa'),
                  ph.inFile('adeno-reads%d_1.fa' % rl),
                  ph.inFile('adeno-reads%d_2.fa' % rl),
                  '-o', ph.outFile('pe-adeno-reads%d_2-id.razers' % rl)],
            to_diff=[(ph.inFile('pe-adeno-reads%d_2-id.razers' % rl),
                      ph.outFile('pe-adeno-reads%d_2-id.razers' % rl)),
                     (ph.inFile('pe-adeno-reads%d_2-id.stdout' % rl),
                      ph.outFile('pe-adeno-reads%d_2-id.stdout' % rl))])
        conf_list.append(conf)

        # Compute forward/reverse matches only.
        for o in ['-r', '-f']:
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('pe-adeno-reads%d_2-id%s.stdout' % (rl, o)),
                args=['--low-memory', '-id', o,
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      ph.inFile('adeno-reads%d_2.fa' % rl),
                      '-o', ph.outFile('pe-adeno-reads%d_2-id%s.razers' % (rl, o))],
                to_diff=[(ph.inFile('pe-adeno-reads%d_2-id%s.razers' % (rl, o)),
                          ph.outFile('pe-adeno-reads%d_2-id%s.razers' % (rl, o))),
                         (ph.inFile('pe-adeno-reads%d_2-id%s.stdout' % (rl, o)),
                          ph.outFile('pe-adeno-reads%d_2-id%s.stdout' % (rl, o)))])
            conf_list.append(conf)

        # Compute with different identity rates.
        for i in range(90, 101):
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('pe-adeno-reads%d_2-id-i%d.stdout' % (rl, i)),
                args=['--low-memory', '-id', '-i', str(i),
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      ph.inFile('adeno-reads%d_2.fa' % rl),
                      '-o', ph.outFile('pe-adeno-reads%d_2-id-i%d.razers' % (rl, i))],
                to_diff=[(ph.inFile('pe-adeno-reads%d_2-id-i%d.razers' % (rl, i)),
                          ph.outFile('pe-adeno-reads%d_2-id-i%d.razers' % (rl, i))),
                         (ph.inFile('pe-adeno-reads%d_2-id-i%d.stdout' % (rl, i)),
                          ph.outFile('pe-adeno-reads%d_2-id-i%d.stdout' % (rl, i)))])
            conf_list.append(conf)

        # Compute with different output formats.
        for suffix in ['razers', 'fa', 'eland', 'gff', 'sam', 'afg']:
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('pe-adeno-reads%d_2-id.%s.stdout' % (rl, suffix)),
                args=['--low-memory', '-id',
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      ph.inFile('adeno-reads%d_2.fa' % rl),
                      '-o', ph.outFile('pe-adeno-reads%d_2-id.%s' % (rl, suffix))],
                to_diff=[(ph.inFile('pe-adeno-reads%d_2-id.%s' % (rl, suffix)),
                          ph.outFile('pe-adeno-reads%d_2-id.%s' % (rl, suffix))),
                         (ph.inFile('pe-adeno-reads%d_2-id.%s.stdout' % (rl, suffix)),
                          ph.outFile('pe-adeno-reads%d_2-id.%s.stdout' % (rl, suffix)))])
            conf_list.append(conf)

        # Compute with different sort orders.
        for so in [0, 1]:
            conf = app_tests.TestConf(
                program=path_to_program,
                redir_stdout=ph.outFile('pe-adeno-reads%d_2-id-so%d.stdout' % (rl, so)),
                args=['--low-memory', '-id', '-so', str(so),
                      ph.inFile('adeno-genome.fa'),
                      ph.inFile('adeno-reads%d_1.fa' % rl),
                      ph.inFile('adeno-reads%d_2.fa' % rl),
                      '-o', ph.outFile('pe-adeno-reads%d_2-id-so%d.razers' % (rl, so))],
                to_diff=[(ph.inFile('pe-adeno-reads%d_2-id-so%d.razers' % (rl, so)),
                          ph.outFile('pe-adeno-reads%d_2-id-so%d.razers' % (rl, so))),
                         (ph.inFile('pe-adeno-reads%d_2-id-so%d.stdout' % (rl, so)),
                          ph.outFile('pe-adeno-reads%d_2-id-so%d.stdout' % (rl, so)))])
            conf_list.append(conf)

    # Execute the tests.
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['razers2'] + conf.args),
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'

    # Cleanup.
    ph.deleteTempDir()

    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    sys.exit(app_tests.main(main))
