
import numpy as np
import argparse


def write_table_files(tables, output_path, output_format='%1.5e'):
  for k in tables.keys():
    tmp = np.reshape(tables[k], (-1), order='F')
    np.savetxt('%s/%s.csv' % (output_path, k), tmp, fmt=output_format, delimiter=',')


def buildInputTables(output_path):
  # Pumping schedule
  q = 3.25 # kg/s
  t_shutDown = 120.0
  epsilon = 0.0001

  # Build the pumping schedule
  tables = {}
  tables['flowRate'] = [q, q, 0.0, 0.0 ]
  tables['flowRate_time'] = [0.0, t_shutDown - epsilon,  t_shutDown + epsilon, 1e9 ]

  # Write the tables
  write_table_files(tables, output_path)

def getParser():
  parser = argparse.ArgumentParser(description='Geneerate input tables for dfit test case')

  parser.add_argument('-d',
                      '--directory',
                      type=str,
                      default='./',
                      help='directory in which input tables are created.')

  return parser

def main(output_path):
  buildInputTables(output_path)
    
if __name__ == "__main__":
  parser = getParser()
  args, unknown_args = parser.parse_known_args()
  if unknown_args:
      print("unknown arguments %s" % unknown_args)

  main(args.directory)
