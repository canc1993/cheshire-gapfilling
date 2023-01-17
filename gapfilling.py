import cobra
from copy import deepcopy
import numpy as np
import pandas as pd
from egc import resolve_egc
import random


# This function adds exchange reactions for targeted fermentation products
def add_ex_reactions(model, universe, target_ex_rxns):
    for rid in target_ex_rxns:
        if rid not in model.reactions:
            #print('add_ex_rxns: adding %s...' % rid)
            reaction = cobra.Reaction(rid)
            reaction.name = 'R_%s' % rid
            reaction.lower_bound = -1000.0
            reaction.upper_bound = 1000.0
            met_id = rid.lstrip('EX_')
            if met_id in model.metabolites:
                met = model.metabolites.get_by_id(met_id)
            else:
                met = universe.metabolites.get_by_id(met_id)
            reaction.add_metabolites({met: -1.0})
            model.add_reactions([reaction])
    return None


# constrain nutrient uptake
def constrain_media(model, df_cm, namespace, rxn_suffix):
    for ex in model.exchanges:
        ex_cpd = ex.id.lstrip('EX_')
        ex_cpd = ex_cpd.rstrip(rxn_suffix)
        if ex_cpd in list(df_cm[namespace]):
            ex.lower_bound = -np.abs(df_cm.loc[df_cm[namespace] == ex_cpd, 'flux'].values[0])
            ex.upper_bound = 1000.0
        else:
            ex.lower_bound = 0.0
            ex.upper_bound = 1000.0
    return None


# test if adding reactions will increase growth rate
def test_growth_inflation(model, rxns_to_add):
    max_growth = model.slim_optimize()
    model2 = deepcopy(model)
    model2.add_reactions(rxns_to_add)
    max_growth2 = model2.slim_optimize()
    if np.abs(max_growth2 - max_growth) <= 1e-6:
        return 0
    else:
        return 1


# add reactions predicted from deep learning model or randomly selected from reaction pools
def add_gapfilled_reactions(gem_file, universe, paras):
    # read parameters
    namespace = paras['NAMESPACE']
    batch_size = int(paras['BATCH_SIZE'])
    num_rxns_to_fill = int(paras['NUM_GAPFILLED_RXNS_TO_ADD'])
    add_random_rxns = int(paras['ADD_RANDOM_RXNS'])
    target_ex_rxns = paras['TARGET_EX_RXNS']
    ex_suffix = paras['EX_SUFFIX']

    # read GEM into cobrapy
    model = cobra.io.read_sbml_model("%s/%s.xml" % (paras['GEM_DIRECTORY'], gem_file))
    model.solver = 'cplex'
    print('add_gapfilled_reaction: model %s loaded' % gem_file)

    # read culture medium
    df_cm = pd.read_csv(paras['CULTURE_MEDIUM'])
    if namespace not in list(df_cm.columns):
        raise RuntimeError('cannot find %s as a column of media file.' % namespace)

    # add exchange reactions
    add_ex_reactions(model, universe, target_ex_rxns)

    # constrain boundary reactions by culture media
    constrain_media(model, df_cm, namespace, ex_suffix)

    # skip the model if it does not grow
    max_growth = model.slim_optimize()
    if max_growth == 0.0:
        raise RuntimeError("model %s cannot grow in current medium." % (gem_file.rstrip('.xml')))
    else:
        print('add_gapfilled_reaction: max growth rate = %2.2f' % max_growth)

    #  model before reaction added
    model_no_gapfill = deepcopy(model)
    model_w_gapfill = deepcopy(model)

    if add_random_rxns:
        # randomize reactions in universal pools
        candidate_reactions = list(set([r.id for r in universe.reactions if r.id not in model_w_gapfill.reactions]))
        random.shuffle(candidate_reactions)
    else:
        # read deep learning model predicted reactions with scores
        score_file = "%s/%s.csv" % (paras['GAPFILLED_RXNS_DIRECTORY'], gem_file)
        df_gapfill = pd.read_csv(score_file, index_col=0)
        # keep gapfilled reactions that are not contained in the model but included in the reaction pools
        df_gapfill = df_gapfill.loc[
            [rid for rid in df_gapfill.index if rid not in model_w_gapfill.reactions and rid in universe.reactions]
        ]
        df_gapfill = df_gapfill[df_gapfill['predicted_scores']>=float(paras['MIN_PREDICTED_SCORES'])].sort_values('similarity_scores')
        candidate_reactions = list(df_gapfill.index)
    print('add_gapfilled_reaction: find %d candidate reactions' % len(candidate_reactions))

    # if Resolve_EGC is True, add reaction one by one; otherwise add all together
    rids_added = []
    num_rxns_before = len(model_w_gapfill.reactions)
    counter = 0
    counter2 = 0
    max_counter = np.min([num_rxns_to_fill, len(candidate_reactions)])

    if int(paras['RESOLVE_EGC']):
        while counter < max_counter:

            # determine which reactions to add for this batch
            rxns_to_add = []
            for rid in candidate_reactions[counter2:]:
                if len(rxns_to_add) >= batch_size or counter + len(rxns_to_add) >= max_counter:
                    break

                assert not rid.startswith('EX_')
                rxn = universe.reactions.get_by_id(rid)

                # anaerobic condition; do not add reactions involving oxygen
                if int(paras['ANAEROBIC']):
                    if namespace == "bigg":
                        if 'o2_c' in [met.id for met in rxn.reactants] or 'o2_c' in [met.id for met in rxn.products]:
                            counter2 += 1
                            continue
                    elif namespace == "modelseed":
                        if 'cpd00007_c0' in [met.id for met in rxn.reactants] or 'cpd00007_c0' in [met.id for met in
                                                                                                   rxn.products]:
                            counter2 += 1
                            continue
                rxns_to_add.append(rxn)
                counter2 += 1

            # add reactions in batch and test growth inflation
            is_inflated = test_growth_inflation(model_w_gapfill, rxns_to_add)
            if is_inflated:
                # growth inflation found, add reaction one by one, and resolve potential egc
                for rxn in rxns_to_add:
                    assert rxn.id not in model_w_gapfill.reactions
                    is_inflated = test_growth_inflation(model_w_gapfill, [rxn])
                    if is_inflated:
                        model_w_gapfill2 = deepcopy(model_w_gapfill)
                        model_w_gapfill2.add_reactions([rxn])
                        egc_resolved = resolve_egc(model_w_gapfill2, rxn.id, namespace)
                        if egc_resolved and model_w_gapfill2.slim_optimize() < 2.81:
                            model_w_gapfill = deepcopy(model_w_gapfill2)
                            rids_added.append(rxn.id)
                            counter += 1
                            print(
                                'add_gapfilled_reaction: added %s after resolving egc, current counter = %d'
                                % (rxn.id, counter))
                        else:
                            print(
                                'add_gapfilled_reaction: failed to add %s and egc not resolved, current counter = %d'
                                % (rxn.id, counter))
                    else:
                        model_w_gapfill.add_reactions([rxn])
                        rids_added.append(rxn.id)
                        counter += 1
                        print(
                            'add_gapfilled_reaction: added %s, no growth inflation, current counter = %d'
                            % (rxn.id, counter))
            else:
                for rxn in rxns_to_add:
                    assert rxn.id not in model_w_gapfill.reactions
                model_w_gapfill.add_reactions(rxns_to_add)
                rids_added.extend([rxn.id for rxn in rxns_to_add])
                counter += len(rxns_to_add)
                print('add_gapfilled_reaction: added %d reactions, no growth inflation, current counter = %d'
                      % (len(rxns_to_add), counter))
    else:
        max_counter = np.min(paras['NUM_GAPFILLED_RXNS_TO_ADD'], len(candidate_reactions))
        rxns_to_add = candidate_reactions[0:max_counter]
        model_w_gapfill.add_reactions([rxns_to_add])

    # calculate number of reactions added
    num_rxns_after = len(model_w_gapfill.reactions)
    num_rxns_added = num_rxns_after - num_rxns_before
    assert num_rxns_added <= max_counter
    print("model %s: successfully added %d reactions." % (gem_file.replace('.xml', ''), num_rxns_added))

    return model_no_gapfill, model_w_gapfill, num_rxns_added, rids_added
