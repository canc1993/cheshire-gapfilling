import cobra
import optlang
import os
import pandas as pd
from copy import deepcopy
from joblib import Parallel, delayed
import read_paras
import fba
import sys
import warnings

warnings.filterwarnings("ignore")


def validate():
    # we currently only support cplex
    solvers = [match.split("_interface")[0] for match in dir(optlang) if "_interface" in match]
    if "cplex" not in solvers:
        raise RuntimeError("cplex not found.")

    # read in arguments
    if len(sys.argv) > 2:
        raise RuntimeError('at most one parameter is supplied.')
    if len(sys.argv) == 2:
        input_file = sys.argv[1]
    else:
        if os.path.exists('input_parameters.txt'):
            input_file = 'input_parameters.txt'
        else:
            raise RuntimeError(
                'input file not specified and the default input_parameters.txt cannot be found as well.'
            )

    # read input parameters
    paras = read_paras.read(input_file)

    # load reaction pools
    universe = cobra.io.read_sbml_model(paras['REACTION_POOL'])
    universe.solver = 'cplex'

    # ***********************************************************************
    # add gapfilled reactions or random reactions selected from reaction pools
    # ***********************************************************************

    # compute fermentation flux
    print('-------------------------------------------------------')
    df_output = None
    retLst = Parallel(n_jobs=int(paras['NUM_CPUS']))(
        delayed(fba.predict_fermentation)(gem_file, universe, paras) for gem_file in paras['GEMs'].split(';')
    )
    if df_output is None:
        df_output = deepcopy(pd.concat(retLst))
    else:
        df_output = pd.concat([df_output, pd.concat(retLst)])

    # write intermediate results to file
    output_file: str = "%s/%s" % (paras['OUTPUT_DIRECTORY'], paras['OUTPUT_FILENAME'])
    df_output.to_csv(output_file, index=False)
    print('done!')
