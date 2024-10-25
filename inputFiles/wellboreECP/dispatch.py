#!/opt/cray/pe/python/3.11.5/bin/python

import argparse
import os
import sys

def main(args):
  parser = argparse.ArgumentParser()
  parser.add_argument('-a','--account',required=True)
  parser.add_argument('-b','--binary',required=True)
  parser.add_argument('-c','--caliper',default="")
  parser.add_argument('-p','--problem',required=True)
  parser.add_argument('-s','--script',default="launch_frontier")
  parser.add_argument('-d','--dispatch',default="sbatch")
  parser.add_argument('-v','--verbosity',type=int,default=0)
  parser.add_argument('-l','--levels',nargs='+')
  args = parser.parse_args(args)

  if args.caliper in ("default", "def", "standard", "std"):
    profiling="\"-t runtime-report,max_column_width=200,calc.inclusive,output=stdout,mpi-report\""
  else:
    profiling=args.caliper

  avail_problems = [ d for d in os.listdir(os.getcwd()) if os.path.isdir(d) ]

  if args.problem not in avail_problems:
    print(f"Unkown problem {args.problem}")
    raise SystemExit(1)

  os.chdir(os.path.join(os.getcwd(),args.problem))
  prob_dir = os.getcwd()

  levels = list( sorted( ( d for d in os.listdir(prob_dir) if os.path.isdir(d) and d.startswith("level") ), key=lambda v: int(v[5:]) ) )
  if args.levels is not None:
    args.levels = [ int(l) for l in args.levels ]
    levels = [ l for l in levels if int(l[5:]) in args.levels ]

  if args.verbosity:
    print(f"{levels = }")

  launch = ' '.join([args.dispatch, "-A", args.account, args.script, args.binary, profiling])
  if args.verbosity:
    print(f"{launch = }")

  for level in levels:
    os.chdir(level)
    print(os.getcwd())
    os.system(launch)
    os.chdir(prob_dir)

if __name__ == "__main__":
  main(sys.argv[1:])
