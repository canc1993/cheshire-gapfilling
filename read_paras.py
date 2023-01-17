import os
import pandas as pd


def read(input_file):
    df = pd.read_csv(input_file, header=None)
    df.columns = ["Field", "Value"]

    # The following fields are mandatory
    mandate_fields: list[str] = [
        'CULTURE_MEDIUM',  # csv file with NAMESPACE and flux
        'REACTION_POOL',   # a universal genome scale model
        'GEM_DIRECTORY',  # directory of genome-scale metabolic models
        'GAPFILLED_RXNS_DIRECTORY',  # reactions predicted to be missing
        'NUM_GAPFILLED_RXNS_TO_ADD',  # number of reactions to add
        'ADD_RANDOM_RXNS', # 1 if adding random reactions
        'SUBSTRATE_EX_RXNS',  # csv file with NAMESPACE
    ]
    for f in mandate_fields:
        if f not in list(df.Field):
            raise RuntimeError("Field %s is mandatory." % f)

    # The following fields are optional
    optional_fields = {
        'MIN_PREDICTED_SCORES': 0.9995, # candidate reactions w/ predicted scores below this cutoff are excluded
        'NUM_CPUS': 1,  # number of cpus to use. use -1 if using all cpus
        'EX_SUFFIX': "_e",  # exchange reaction suffix
        'RESOLVE_EGC': True,  # whether resolve energy-generating cycle
        'OUTPUT_DIRECTORY': "./",  # output file directory
        'OUTPUT_FILENAME': 'suggested_gaps.csv',  # output file name
        'FLUX_CUTOFF': 1e-5,  # cutoff flux for
        'ANAEROBIC': True,  # anaerobic fermentation?
        'BATCH_SIZE': 10,  # number of reactions added in a batch
        'NAMESPACE': 'bigg'  # currently we only support bigg and modelseed
    }

    # Keep only fields that are recognizable
    df = df[df.Field.isin(mandate_fields + list(optional_fields.keys()))]

    # Convert dataframe to dictionary
    paras = {f: v for f, v in zip(df.Field, df.Value)}

    # Assign default values for optional fields
    for f in optional_fields:
        if f not in paras:
            paras[f] = optional_fields[f]

    # make sure culture medium file exists
    if not os.path.exists(paras['CULTURE_MEDIUM']):
        raise RuntimeError('cannot find culture medium file.')

    # NAMESPACE only support 'bigg' and 'modelseed'
    if paras['NAMESPACE'] not in ['bigg', 'modelseed']:
        raise RuntimeError('unrecognized namespace %s.' % (paras['NAMESPACE']))

    # Find overlaps between GEM_DIRECTORY and GAPFILLED_RXNS_DIRECTORY
    filenames_in_GEM_DIRECTORY = [f.rstrip('.xml') for f in os.listdir(paras['GEM_DIRECTORY']) if f.endswith('.xml')]
    filenames_in_GAPFILLED_RXNS_DIRECTORY = [
        f.rstrip('.csv') for f in os.listdir(paras["GAPFILLED_RXNS_DIRECTORY"]) if f.endswith('.csv')
    ]
    overlapped_gems = list(set(filenames_in_GEM_DIRECTORY).intersection(set(filenames_in_GAPFILLED_RXNS_DIRECTORY)))
    if len(overlapped_gems) == 0:
        raise RuntimeError("cannot find gapfilled reactions for any genome-scale model.")
    else:
        # print("found %d matched GEMs and gapfilled reactions." % len(overlapped_gems))
        paras['GEMs'] = ';'.join(overlapped_gems)

    # Read target exchange reactions
    df_ex = pd.read_csv(paras['SUBSTRATE_EX_RXNS'], index_col=0)
    paras['TARGET_EX_RXNS'] = ['EX_' + cpd + paras['EX_SUFFIX'] for cpd in df_ex[paras['NAMESPACE']] if
                               str(cpd) != 'nan']

    return paras
