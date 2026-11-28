#!/usr/bin/env python3
"""Conditional parameter assignment is not a below-coverage activation gate."""
import argparse
import importlib.util
from pathlib import Path
import re
import shutil
import subprocess

HERE=Path(__file__).resolve().parent
WORK=HERE/".spp-grasp-coverage-profile-selftest"


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--cbc",default=str(HERE.parent/"src/cbc"))
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location("features",HERE/"spp-grasp-adaptive-features-selftest.py")
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    if WORK.exists():raise RuntimeError("refusing to overwrite "+str(WORK))
    WORK.mkdir()
    try:
        full=WORK/"full.mps"
        partial=WORK/"partial.mps"
        below=WORK/"below.mps"
        forced=WORK/"forced.mps"
        packed=WORK/"packed.mps"
        fractional=WORK/"fractional-native.mps"
        helper.write_mps(full,["E","E"],[(1,[(0,1)]),(2,[(1,1)]),(0,[(0,1),(1,1)])])
        helper.write_mps(partial,["E","E"],[(1,[(0,1)]),(0,[(0,1),(1,1)]),(2,[(0,1),(1,1)])])
        helper.write_mps(below,["E","E","E"],[(1,[(0,1)]),(0,[(0,1),(1,1),(2,1)]),(2,[(0,1),(1,1),(2,1)])])
        helper.write_mps(forced,["E","E"],[(0,[(0,1)]),(1,[(0,1),(1,1)]),(2,[(1,1)])],fixed=(0,))
        helper.write_mps(packed,["E","E","L"],[(1,[(0,1),(2,1)]),(1,[(1,1),(2,1)]),(2,[(0,1),(1,1)])])
        helper.write_mps(fractional,["E"]*3,[
            (0,[(0,1),(1,1)]),(0,[(1,1),(2,1)]),(0,[(0,1),(2,1)]),
            (1,[(0,1)]),(1,[(1,1)]),(1,[(2,1)]),
        ])
        profile=WORK/"tuned.params"
        profile.write_text("# experimental qualifying calls only\ndrows 1\ndgrow 0.2\nartificial 2\nmaxwork 100000\nnoimp 3\n")
        def solve(model,profilefile=None,extra=(),budget=True):
            cmd=[args.cbc,str(model),"-threads","1","-log","2","-randomSeed","17",
                 "-randomCbcSeed","17","-sec","5","-preprocess","off","-preRootLPStrenghtening","off",
                 "-feas","on"]
            if budget:cmd+=["-sppGraspWork",".001","-sppGraspNoImprove","1"]
            if profilefile:cmd+=["-sppGraspTuningProfile",str(profilefile)]
            result=subprocess.run([*cmd,*extra,"-solve","-quit"],capture_output=True,text=True,
                                  timeout=15,check=True)
            text=result.stdout+result.stderr
            assert "No match for" not in text,text
            return text
        for model in (full,partial):
            text=solve(model,profile)
            assert "SppGrasp coverage tuning applied:" in text,text
            assert "stage original; native coverage" in text,text
            assert "; drows 1; dgrow 0.2;" in text,text
            # Existing explicit budget wins over the value in the file.
            assert "; maxWork 1000; noImprove 1" in text,text
        ordinary=solve(full)
        assert "SppGrasp run result: stage original; work " in ordinary,ordinary
        phased=WORK/"phased.params"
        phased.write_text("lls 3\nartificial 2\n")
        for model in (full,partial):
            text=solve(model,phased)
            assert "SppGrasp coverage tuning applied:" in text,text
            assert "; lnsLocalSearch 3;" in text,text
        text=solve(full,profile,("-sppGraspDestroyRows","7","-sppGraspDestroyGrow","0"))
        assert "; drows 7; dgrow 0;" in text,text
        root=solve(fractional,profile,("-sppGraspRootSequence","spp-fpump",
            "-sppGraspRootKeepOriginal","off","-cuts","off"))
        assert "SppGrasp coverage tuning applied: stage root;" in root,root
        assert "SppGrasp coverage tuning result: stage root;" in root,root
        assert "stage original;" not in root,root
        excessive=WORK/"excessive.params"
        excessive.write_text("totalwork 1e12\nnoimp 0\n")
        short=solve(fractional,excessive,("-sppGraspRootSequence","spp-fpump-short",
            "-sppGraspStageFirstPercent","1"),budget=False)
        counts=re.findall(r"coverage tuning result: stage original; work ([^;]+)",short)
        assert counts and max(map(float,counts))<=15010000,short
        def trajectory(text):
            return re.findall(r"SppGrasp: (?:first .*?|no feasible solution .*?)after \d+ LNS iterations"
                              r"|Objective value:\s+\S+|Enumerated nodes:\s+\d+",text)
        for model in (below,forced,packed):
            baseline=solve(model)
            tuned=solve(model,profile)
            assert "SppGrasp coverage tuning legacy fallback:" in tuned,tuned
            assert trajectory(tuned)==trajectory(baseline),(baseline,tuned)
            assert "SppGrasp: first" in tuned or "SppGrasp: no feasible solution" in tuned,tuned
        text=solve(full,profile,("-sppGrasp","off"))
        assert "SppGrasp coverage tuning" not in text and "SppGrasp: first" not in text,text
        text=solve(partial,profile,("-sppGraspArtificial","0"))
        assert "SppGrasp coverage tuning applied" not in text and "SppGrasp: first" not in text,text
        for invalid in ("drows -1\n","unknown 4\n","lp nan\n","drows 2 extra\n",
                        "lls 6\n","lls -2\n"):
            bad=WORK/"bad.params";bad.write_text(invalid)
            text=solve(full,bad)
            assert "invalid SppGrasp coverage tuning profile" in text,text
            assert trajectory(text)==trajectory(solve(full)),text
        missing=solve(full,WORK/"missing.params")
        assert "cannot read SppGrasp coverage tuning profile" in missing,missing
        directory=solve(full,WORK)
        assert "cannot read SppGrasp coverage tuning profile" in directory,directory
        incompatible=WORK/"incompatible.params"
        incompatible.write_text("artificial 0\ndrows 1\n")
        text=solve(partial,incompatible)
        assert "incompatible artificial setting, retaining legacy model" in text,text
        assert trajectory(text)==trajectory(solve(partial)),text
        print("PASS: original/root full/partial profiles, below50/forced/packing legacy, explicit overrides/OFF, scheduler budgets and invalid-file fallback")
    finally:
        shutil.rmtree(WORK)


if __name__=="__main__":main()
