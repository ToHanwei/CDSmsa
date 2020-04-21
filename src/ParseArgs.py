#!coding:utf-8

import argparse
from multiprocessing import cpu_count

from src.config import TOOLNAME
from src.config import VERSION


class ParseCommand(object):
	"""Parse command line"""

	def parsecommand(self):
		"""FindOR.py command line parse"""

		parser = argparse.ArgumentParser(prog=TOOLNAME,
										 description="Olfactory receptor annotation",
										 epilog="http://zhaolab.shanghaitech.edu.cn/"
										 )
		# positional arguments
		parser.add_argument("input",
							help="String, nhmmer output file path.",
							)
		parser.add_argument('genome',
							help="String, Genomic data file path.",
							)

		# optional arguments
		parser.add_argument('-o', '--outputdir',
							help="String, Result save directory.(default:../output)",
							default="../output",
							metavar='',
							)
		parser.add_argument('-p', '--prefix',
							help="String, output file prefix.(default:ORannotation)",
							default="ORannotation",
							metavar='',
							)
		parser.add_argument('-e', '--EvalueLimit',
							type=float,
							help="Float, Sequence similarity threshold.\
								 (default:1e-60)",
							# An empirical value
							default=1e-60,
							metavar='',
							)
		parser.add_argument('-l', '--SeqLengthLimit',
							type=int,
							help="Int, An artificially set OR's sequence \
								 length threshold.(default:868)",
							# The OR sequence length is approximately 310
							default=869,
							metavar='',
							)
		parser.add_argument('-v', '--verbose',
							action='count',
							help="Print verbose information.",
							)
		parser.add_argument('-V', '--version',
							action='version',
							help='Show version message and exit.',
							version=VERSION,
							)
		args = parser.parse_args()
		return args


