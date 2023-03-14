import cobra
from cobra.flux_analysis import flux_variability_analysis
from optlang.interface import OPTIMAL
import pandas as pd
import numpy as np
from gapfilling import add_gapfilled_reactions
from optlang.symbolics import add
from copy import deepcopy


def flux_balance_analysis(model, paras):
    # run flux balance analysis
    # try different linear programming method the default algorithm fails
    fba_solution = model.optimize()
    if model.solver.status != OPTIMAL:
        is_optimal = False
        for lp_method in ["primal", "dual", "network", "barrier", "sifting", "concurrent"]:
            model.solver.configuration.lp_method = lp_method
            fba_solution = model.optimize()
            if model.solver.status == OPTIMAL:
                is_optimal = True
                break
        assert is_optimal
    assert fba_solution.objective_value > 0.0

    # run parsimonious flux balance analysis
    # try different linear programming method the default algorithm fails
    pfba_solution = cobra.flux_analysis.pfba(model)

    # modify flux bounds to minimize input fluxes that do not contribute to growth
    # For flux > 0, set its lower bound to 0
    # For flux <= 0, set its lower bound to the flux value
    for ex in model.exchanges:
        assert ex.id.startswith('EX_')
        if pfba_solution.fluxes[ex.id] >= 0.0:
            ex.lower_bound = 0.0
        else:
            ex.lower_bound = pfba_solution.fluxes[ex.id]

    # run flux variability analysis
    fva = flux_variability_analysis(
        model,
        paras['TARGET_EX_RXNS'],
        fraction_of_optimum=0.999999,
        loopless=True
    )
    fva.index.name = 'reaction'
    fva = fva.reset_index()
    fva['biomass'] = fba_solution.objective_value
    fva['normalized_maximum'] = fva['maximum'] / fva['biomass']
    fva['phenotype'] = (fva['normalized_maximum'] >= float(paras['FLUX_CUTOFF'])).astype(int)
    return fva


def predict_fermentation(gem_file, universe, paras):
    print('predicting fermentation: %s...' % gem_file)
    # read model and add missing reactions
    model_no_gapfill, model_w_gapfill, num_rxns_added, rids_added = add_gapfilled_reactions(gem_file, universe, paras)

    # run flux balance analysis for model with and without gap filling
    fva_no_gapfill = flux_balance_analysis(model_no_gapfill, paras)
    fva_no_gapfill.columns = [c + '__no_gapfill' if c != 'reaction' else 'reaction' for c in fva_no_gapfill.columns]
    fva_w_gapfill = flux_balance_analysis(model_w_gapfill, paras)
    fva_w_gapfill.columns = [c + '__w_gapfill' if c != 'reaction' else 'reaction' for c in fva_w_gapfill.columns]
    fva = pd.merge(fva_no_gapfill, fva_w_gapfill, left_on='reaction', right_on='reaction', how='inner')

    # expand fva
    fva['gem_file'] = gem_file.rstrip('.xml')
    fva['random_rxns'] = paras['ADD_RANDOM_RXNS']
    fva['num_rxns_to_add'] = paras['NUM_GAPFILLED_RXNS_TO_ADD']
    fva['num_rxns_added'] = num_rxns_added
    fva['rxn_ids_added'] = ';'.join(rids_added)

    # find reactions that lead to phenotypic changes from 0 to 1
    key_rxns = []
    for ex, phe1, phe2 in zip(fva['reaction'], fva['phenotype__no_gapfill'], fva['phenotype__w_gapfill']):
        if phe1 == 0 and phe2 == 1:
            model_w_gapfill2 = deepcopy(model_w_gapfill)
            model_w_gapfill2.reactions.get_by_id(ex).lower_bound = 0.1  # some nontrivial small number
            indicator_vars = []
            for rid in rids_added:
                var = model_w_gapfill2.problem.Variable('indicator_var_' + rid, lb=0, ub=1, type='binary')
                indicator_vars.append(var)
                con1 = model_w_gapfill2.problem.Constraint(
                    (model_w_gapfill2.reactions.get_by_id(rid).flux_expression + 1000.0 * var).expand(),
                    name='constr1' + rid,
                    lb=0
                )
                con2 = model_w_gapfill2.problem.Constraint(
                    (model_w_gapfill2.reactions.get_by_id(rid).flux_expression - 1000.0 * var).expand(),
                    name='constr2' + rid,
                    ub=0
                )
                model_w_gapfill2.add_cons_vars([var, con1, con2])
                model_w_gapfill2.solver.update()

            model_w_gapfill2.objective = add(*indicator_vars)
            model_w_gapfill2.objective.direction = "min"
            model_w_gapfill2.solver.update()

            skip_this_prediction = False
            for tol in [1e-9, 1e-8, 1e-7, 1e-6]:
                model_w_gapfill2.solver.problem.parameters.mip.tolerances.integrality.set(tol)
                try:
                    sol = model_w_gapfill2.optimize()
                    break
                except:
                    if tol == 1e-6:
                        skip_this_prediction = True
            if skip_this_prediction:
                key_rxns.append(np.NaN)
            else:
                if sol.objective_value >= 1:
                    key_rxns2 = []
                    for var in model_w_gapfill2.variables:
                        rid = var.name.lstrip('indicator_var_')
                        if var.primal == 1:
                            print(ex, var.name)
                        if var.primal == 1 and rid in rids_added:
                            key_rxns2.append(rid)
                    key_rxns.append(';'.join(key_rxns2))
                    print('predict_fermentation: %s can be gapfilled by %s' % (ex, key_rxns[-1]))
                else:
                    key_rxns.append(np.NaN)
        else:
            key_rxns.append(np.NaN)
    fva['essential_rxns'] = key_rxns

    return fva
