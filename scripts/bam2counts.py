#!/usr/bin/env python

import pysam
from tqdm import tqdm
import pyfastx
import os
import pickle
from collections import namedtuple
import sys
import json
from datetime import datetime
import argparse  # New import

acore = namedtuple('acore', ['rst', 'ren', 'score'])

def load_aln(fn, ti_map):
    aln_mat = []
    ctr = -1
    qi_map = dict()
    pri_lst = []
    unmapped = set()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        for brec in tqdm(fh):
            if brec.is_unmapped:
                unmapped.add(brec.query_name)
                continue
            elif brec.is_supplementary:
                continue
            else:
                qname = brec.query_name
                if qname not in qi_map:
                    ctr += 1
                    qi_map[qname] = ctr
                    aln_mat.append(dict())
                    pri_lst.append(-1)
                    qi = ctr
                else:
                    qi = qi_map[qname]
                ti = ti_map[brec.reference_name]
                curr_acore = acore(rst=brec.reference_start, \
                                    ren=brec.reference_end, \
                                    score=int(brec.get_tag("ms")))
                if ti in aln_mat[qi]:
                    if curr_acore.score > aln_mat[qi][ti].score:
                        aln_mat[qi][ti] = curr_acore
                else:
                    aln_mat[qi][ti] = curr_acore

                if not brec.is_secondary:
                    assert pri_lst[qi] == -1
                    pri_lst[qi] = ti
    return aln_mat, pri_lst, qi_map, unmapped               


def build_ti_map(fn):
    fa = pyfastx.Fasta(fn)
    ti_map = dict()
    tlen_lst = []
    ctr = 0
    for obj in fa:
        ti_map[obj.name] = ctr
        tlen_lst.append(len(obj.seq))
        ctr += 1
    return ti_map, tlen_lst

def build_cmpt_mat(aln_mat, qi_map, pri_lst, tlen_lst, \
                filter, fp_thres, tp_thres, surrender, tcov_thres):
    qi_size = len(qi_map)
    cmpt_mat = [dict() for _ in range(qi_size)]
    for qi in tqdm(range(qi_size)):
        pri_ti = pri_lst[qi]
        if pri_ti == -1:
            # set the alingment with the highest score as the primary alignment
            print("FATAL ERROR: no primary aln for a read")
            sys.exit(-1)
        pri_rst, pri_ren, pri_score = aln_mat[qi][pri_ti]
        pri_tlen = tlen_lst[pri_ti]
        if surrender:
            tcov = (pri_ren - pri_rst) / tlen_lst[pri_ti]
            if tcov < tcov_thres:
                continue
        cmpt_mat[qi][pri_ti] = 1.0
        for ti in aln_mat[qi]:
            if ti == pri_ti:
                continue
            acore = aln_mat[qi][ti]

            if filter:
                fp_dist = pri_rst - acore.rst
                if tp_thres == -1:
                    if fp_dist < fp_thres:
                        continue
                tlen = tlen_lst[ti]
                tp_dist = (pri_tlen - pri_ren) - (tlen - acore.ren)
                if fp_dist < fp_thres or tp_dist < tp_thres:
                    continue
            cmpt_score = acore.score / pri_score

            if ti in cmpt_mat[qi]:
                cmpt_mat[qi][ti] = max(cmpt_mat[qi][ti], cmpt_score)
            else:
                cmpt_mat[qi][ti] = cmpt_score
    return cmpt_mat


def main():
    # Argument parsing setup
    parser = argparse.ArgumentParser(description="Process alignment and build compatibility matrix.")
    parser.add_argument('--aln', required=True, help="Path to the alignment BAM file")
    parser.add_argument('--target', required=True, help="Path to the target transcriptome FASTA file")
    parser.add_argument('--out_dir', required=True, help="Directory for output files")
    parser.add_argument('--filter', action='store_true', help="Enable filtering")
    parser.add_argument('--five_prime', type=int, default=500, help="5' end threshold for filtering")
    parser.add_argument('--three_prime', type=int, default=500, help="3' end threshold for filtering")
    parser.add_argument('--surrender', action='store_true', help="Enable surrender for low target coverage")
    parser.add_argument('--target_cover', type=float, default=0.8, help="Target coverage threshold")

    # Parse arguments
    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(args.out_dir, exist_ok=True)

    # Define the file path
    cmd_fn = os.path.join(args.out_dir, "pref_cmd_info.json")

    # Write command-line arguments to the JSON file
    with open(cmd_fn, 'w') as f:
        json.dump(vars(args), f, indent=2)

    print(datetime.now(), f"PROGRESS loading target transcriptome")
    ti_map, tlen_lst = build_ti_map(args.target)

    print(datetime.now(), f"PROGRESS loading query-to-target alignments")
    aln_mat, pri_lst, qi_map, unmapped = load_aln(args.aln, ti_map)

    # Output the unmapped reads
    unm_fn = os.path.join(args.out_dir, "unmapped.txt")
    with open(unm_fn, 'w') as unm_fh:
        for qname in unmapped:
            unm_fh.write(qname + "\n")

    # Save the target information map
    ti_map_fn = os.path.join(args.out_dir, "ti.pkl")
    with open(ti_map_fn, "wb") as f:
        pickle.dump(ti_map, f)

    print(datetime.now(), f"PROGRESS building a compatibility matrix")
    cmpt_mat = build_cmpt_mat(aln_mat, qi_map, pri_lst, tlen_lst, \
                              args.filter, args.five_prime, args.three_prime, \
                              args.surrender, args.target_cover)
    
    # Save the compatibility matrix and score results
    cmpt_fn = os.path.join(args.out_dir, "cmpt_mat.tsv")
    score_fn = os.path.join(args.out_dir, "scores.tsv")
    sur_fn = os.path.join(args.out_dir, "unassigned.txt")

    qnames = sorted(qi_map, key=lambda k: qi_map[k])
    tnames = sorted(ti_map, key=lambda k: ti_map[k])

    with open(cmpt_fn, 'w') as cmpt_fh, open(score_fn, 'w') as score_fh, open(sur_fn, 'w') as sur_fh:
        cmpt_fh.write("#" + str(len(cmpt_mat)) + "\n")
        score_fh.write("#" + str(len(cmpt_mat)) + "\n")
        for qi in range(len(qi_map)):
            qname = qnames[qi]
            if len(cmpt_mat[qi]) == 0:
                sur_fh.write(qname + "\n")
                continue
            cmpt_fh.write(qname)
            score_fh.write(qname)
            for ti in cmpt_mat[qi]:
                cmpt_fh.write("\t" + tnames[ti])
                score_fh.write("\t(" + tnames[ti] + "," + str(cmpt_mat[qi][ti]) + ")")
            cmpt_fh.write("\n")
            score_fh.write("\n")

if __name__ == "__main__":
    main()
