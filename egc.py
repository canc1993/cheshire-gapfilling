from copy import deepcopy
from cobra import Reaction
import warnings

warnings.filterwarnings("ignore")

# energy couplers
energy_couples_modelseed = {
    'cpd00002_c0': 'cpd00008_c0',
    'cpd00052_c0': 'cpd00096_c0',
    'cpd00038_c0': 'cpd00031_c0',
    'cpd00062_c0': 'cpd00014_c0',
    'cpd00068_c0': 'cpd00090_c0',
    'cpd00005_c0': 'cpd00006_c0',
    'cpd00004_c0': 'cpd00003_c0',
    'cpd00982_c0': 'cpd00015_c0',
    'cpd01270_c0': 'cpd00050_c0',
    'cpd15561_c0': 'cpd15560_c0',
    'cpd15499_c0': 'cpd15500_c0',
    'cpd15994_c0': 'cpd15995_c0',
    'cpd23255_c0': 'cpd11606_c0',
    'cpd15353_c0': 'cpd15352_c0',
    'cpd00022_c0': 'cpd00010_c0',
    'cpd00023_c0': 'cpd00024_c0',
    'cpd00067_p0': 'cpd00067_c0'
}
energy_couples_bigg = {
    'atp_c': 'adp_c',
    'ctp_c': 'cdp_c',
    'gtp_c': 'gdp_c',
    'utp_c': 'udp_c',
    'itp_c': 'idp_c',
    'nadph_c': 'nadp_c',
    'nadh_c': 'nad_c',
    'fadh2_c': 'fad_c',
    'fmnh2_c': 'fmn_c',
    'q8h2_c': 'q8_c',
    'mql8_c': 'mqn8_c',
    'mql6_c': 'mqn6_c',
    'mql7_c': 'mqn7_c',
    '2dmmql8_c': '2dmmq8_c',
    'accoa_c': 'coa_c',
    'glu__L_c': 'akg_c',
    'h_p': 'h_c'
}


# This function identifies energy-generating cycle by maximizing flux towards a dissipation reaction
# If rid is not none, it returns EGC with non-zero flux through rid
def detect_egc(model, metabolite_id, rid=None, namespace='bigg'):
    # close all boundary reactions (including exchange reactions, demand reactions and sink reactions)
    for boundary in model.boundary:
        boundary.bounds = (0, 0)

    # get dissipation product
    energy_couples = energy_couples_bigg
    if namespace == 'modelseed':
        energy_couples = energy_couples_modelseed
    met = model.metabolites.get_by_id(metabolite_id)
    if energy_couples[met.id] in model.metabolites:
        dissipation_product = model.metabolites.get_by_id(energy_couples[met.id])
    else:
        return None

    # add a demand reaction
    dissipation_rxn = Reaction('Dissipation')
    model.add_reactions([dissipation_rxn])
    if namespace == 'bigg':
        if met.id in ['atp_c', 'ctp_c', 'gtp_c', 'utp_c', 'itp_c']:
            # build nucleotide-type dissipation reaction
            dissipation_rxn.reaction = "h2o_c --> h_c + pi_c"
        elif met.id in ['nadph_c', 'nadh_c']:
            # build nicotinamide-type dissipation reaction
            dissipation_rxn.reaction = "--> h_c"
        elif met.id in ['fadh2_c', 'fmnh2_c', 'q8h2_c', 'mql8_c',
                        'mql6_c', 'mql7_c', 'dmmql8_c']:
            # build redox-partner-type dissipation reaction
            dissipation_rxn.reaction = "--> 2 h_c"
        elif met.id == 'accoa_c':
            dissipation_rxn.reaction = "h2o_c --> h_c + ac_c"
        elif met.id == "glu__L_c":
            dissipation_rxn.reaction = "h2o_c --> 2 h_c + nh3_c"
        elif met.id == 'h_p':
            pass
    elif namespace == 'modelseed':
        if met.id in ['cpd00002_c0', 'cpd00052_c0', 'cpd00038_c0', 'cpd00062_c0', 'cpd00068_c0']:
            # build nucleotide-type dissipation reaction
            dissipation_rxn.reaction = "cpd00001_c0 --> cpd00067_c0 + cpd00009_c0"
        elif met.id in ['cpd00005_c0', 'cpd00004_c0']:
            # build nicotinamide-type dissipation reaction
            dissipation_rxn.reaction = "--> cpd00067_c0"
        elif met.id in ['cpd00982_c0', 'cpd01270_c0', 'cpd15561_c0', 'cpd15499_c0',
                        'cpd15994_c0', 'cpd23255_c0', 'cpd15353_c0']:
            # build redox-partner-type dissipation reaction
            dissipation_rxn.reaction = "--> 2 cpd00067_c0"
        elif met.id == 'cpd00022_c0':
            dissipation_rxn.reaction = "cpd00001_c0 --> cpd00067_c0 +  cpd00029_c0"
        elif met.id == 'cpd00023_c0':
            dissipation_rxn.reaction = "cpd00001_c0 --> 2 cpd00067_c0 +  cpd00013_c0"
        elif met.id == 'cpd00067_p0':
            pass
    dissipation_rxn.add_metabolites(
        {met.id: -1, dissipation_product: 1})
    model.objective = dissipation_rxn

    # it is possible that an EGC solution has zero flux through the target reaction
    # to check if a target reaction mediates formation of EGC, force fluxes through this reaction
    if rid is None:
        solution = model.optimize()
    else:
        # find EGC when forcing flux through the target reaction
        original_bounds = model.reactions.get_by_id(rid).bounds
        solution_w_enforced_flux_found = False
        for enforced_flux in ['negative', 'positive']:
            if enforced_flux == 'negative':
                if original_bounds[0] <= -0.01:
                    model.reactions.get_by_id(rid).bounds = [-1000.0, -0.01]
            if enforced_flux == 'positive:':
                if original_bounds[1] >= 0.01:
                    model.reactions.get_by_id(rid).bounds = [0.01, 1000.0]
            solution = model.optimize()
            if solution.status == "optimal":
                solution_w_enforced_flux_found = True
            else:
                model.reactions.get_by_id(rid).bounds = original_bounds
        if not solution_w_enforced_flux_found:
            solution = model.optimize()

    # return results
    if solution.status == 'infeasible':
        raise RuntimeError(
            "The model cannot be solved as the solver status is"
            "infeasible. This may be a bug."
        )
    elif solution.objective_value > 0.0:
        return solution.fluxes[solution.fluxes.abs() > 0.0].to_frame().drop(["Dissipation"])
    else:
        return None


# This function detects and resolves EGC formed when adding reaction rid
def resolve_egc(model, rid, namespace):
    energy_couples = energy_couples_bigg
    if namespace == "modelseed":
        energy_couples = energy_couples_modelseed
    for key_met in energy_couples.keys():
        if key_met in model.metabolites:  # energy couple is part of the model
            df_sol = detect_egc(deepcopy(model), key_met, rid, namespace)
            if df_sol is not None and rid in list(df_sol.index):  # EGC found and involves reaction rid
                rxn = model.reactions.get_by_id(rid)

                # The principle for resolving EGC:
                # 1. If the reaction is irreversible, then we cannot add this reaction. Skip this reaction.
                # 2. If the reaction is reversible, turn it to irreversible and detect EGC again.
                # 3. If EGC remains, then we have to skip this reaction
                if rxn.lower_bound == 0.0 or rxn.upper_bound == 0.0:  # irreversible
                    return False
                else:
                    rid_flux = df_sol.loc[rid, 'fluxes']  # flux through the reaction rid
                    if rid_flux > 0:
                        rxn.upper_bound = float(0.0)
                    else:
                        rxn.lower_bound = float(0.0)
                    df_sol = detect_egc(deepcopy(model), key_met, rid, namespace)
                    if df_sol is not None:  # EGC remains
                        return False

    return True
