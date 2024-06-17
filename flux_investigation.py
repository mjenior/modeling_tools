#!/usr/bin/env python

import os
import math
import pandas
import subprocess

import cobra, cobra.io
from cobra.flux_analysis import production_envelope, pfba
import escher

import matplotlib.pyplot as matplot
from graphviz import Digraph

#-------------------------------------------#
# draw_pathway

class Reaction(object):
    """Reaction class"""

    def __init__(self, rid=None, flux=1, dGr=0, metabolites={}, equation='', reversible=True):
        self.rid = rid
        self.flux = flux
        self.dGr = dGr
        self.metabolites = metabolites
        self.equation = equation
        self.reversible = reversible

    @property
    def reactants(self):
        if self.flux > 0:
            return [k for k,v in self.metabolites.items() if v < 0]
        else:
            return [k for k,v in self.metabolites.items() if v > 0]
    @property
    def products(self):
        if self.flux > 0:
            return [k for k,v in self.metabolites.items() if v > 0]
        else:
            return [k for k,v in self.metabolites.items() if v < 0]


    @classmethod
    def create_Reaction_list_from_dict(cls, dataDict, excludeExchangeRxn=True):

        RxnObjList = []
        for i in range(len(dataDict['reaction_id'])):
            if excludeExchangeRxn:
                if 'EX_' in dataDict['reaction_id'][i] or 'Ex_' in dataDict['reaction_id'][i]:
                    continue
            if dataDict['dGr']:
                tempRxn = cls(dataDict['reaction_id'][i], dataDict['flux'][i], dataDict['dGr'][i])
            else:
                tempRxn = cls(dataDict['reaction_id'][i], dataDict['flux'][i])

            #Get the metabolites dictionary {'C00001': -1, ...} for each reaction
            tempRxn.metabolites = draw_pathway.rxn_dict[tempRxn.rid]['sji']
            RxnObjList.append(tempRxn)
        return RxnObjList

class Pathway(object):
    """OptStoic Pathway class"""

    def __init__(self, id=None, name=None, reaction_ids=[], fluxes=[], dGr_lst=[], reactions=None,
                sourceSubstrateID=None, endSubstrateID=None, note=''):
        """
        id: pathway id
        name: pathway name
        reaction_ids: list of reaction IDs (kegg_id) in the pathway
        fluxes: list of reaction fluxes corresponding to the reaction_ids
        dGr: list of dGr for a reaction ideally from mdf calculation
        reactions : list of reaction object that form the pathway
        sourceSubstrateID: Kegg compound ID of the source metabolite of the pathway
        endSubstrateID: Kegg compound ID of the end metabolite of the pathway
        """
        self.id = id
        self.name = name
        self.note = note

        #iniatilize pathway object using list of reaction_ids, fluxes and dGr
        if reactions is None:
            #check if reaction_ids, fluxes and dGrs are list and contain the same number of item
            assert (isinstance(reaction_ids, list) == 1)
            assert (isinstance(fluxes, list) == 1)
            assert (isinstance(dGr_lst, list) == 1)
            if fluxes:
                assert len(reaction_ids) == len(fluxes), "number of reactions must equal number of fluxes!"
                # convert all flux to float
                for i in range(len(fluxes)):
                    try:
                        fluxes[i] = float(fluxes[i])
                    except:
                        raise Exception(reaction_ids[i]+' flux is not float')
            if dGr_lst:
                assert len(reaction_ids) == len(dGr_lst), "number of reactions must equal number of dGr!"
                # convert all dGr to float
                for i in range(len(dGr_lst)):
                    try:
                        dGr_lst[i] = float(dGr_lst[i])
                    except:
                        raise Exception(reaction_ids[i]+' dGr is not float')


            #change EX_h+ to EX_hplus as optstoic pulp fail to read "+" symbols
            self.reaction_ids = ["EX_hplus" if x == "EX_h+" else x for x in reaction_ids]

            self.fluxes = fluxes
            self.dGr_lst = dGr_lst
            #create list of reaction objects
            self.reactions = Reaction.create_Reaction_list_from_dict({'reaction_id': self.reaction_ids, 'flux': self.fluxes, 'dGr': self.dGr_lst})
        #iniatilize pathway object using list of reaction objects
        else:
            self.reactions = reactions
            self.fluxes = [r.flux for r in reactions]
            self.dGr_lst = [r.dGr for r in reactions]
            self.reaction_ids = [r.rid for r in reactions]

        self.sourceSubstrateID = sourceSubstrateID
        self.endSubstrateID = endSubstrateID

color_configs = {}
# for light background
color_configs['light'] = dict(COFACTOR_SHAPE= "ellipse", # "diamond"
                            OTHER_COFACTOR_COLOR = '#B6B6B6',
                            NONCOFACTOR_SHAPE = "plaintext", # "box"
                            NONCOFACTOR_COLOR = 'transparent', #"#D2EBEB"
                            REACTION_COLOR = "#512DA8",
                            RXN_NODE_COLOR = "#323232",
                            EDGE_COLOR = "#323232", #"#505050"
                            BACKGROUND_COLOR = "transparent",
                            ALL_FONT_COLOR = "black")

#for dark background
color_configs['dark'] = dict(COFACTOR_SHAPE = "ellipse", # "diamond"
                            OTHER_COFACTOR_COLOR = "#7F7F7F",
                            NONCOFACTOR_SHAPE = "plaintext", # "box"
                            NONCOFACTOR_COLOR = "transparent", # "#CCFF33"
                            REACTION_COLOR = "#FFFF00",
                            EDGE_COLOR = "#E5E5E5", #"#505050"
                            RXN_NODE_COLOR = "#E5E5E5",
                            BACKGROUND_COLOR = "transparent",
                            ALL_FONT_COLOR = "white")

def _load_global_styles(colorConfig):

	#for light background
    color_configs['light']['colorMapping'] = {draw_pathway.cfDict['ATP']:'#FF5252', draw_pathway.cfDict['ADP']:'#FF5252',\
                        draw_pathway.cfDict['NAD+']:'#FFEB3B', draw_pathway.cfDict['NADH']:'#FFEB3B',\
                       draw_pathway.cfDict['NADPH']:'#F8C674', draw_pathway.cfDict['NADP+']:'#F8C674', draw_pathway.cfDict['NH4(+)']:'#BBDEFB'\
                        }

    #for dark background
    color_configs['dark']['colorMapping'] = {draw_pathway.cfDict['ATP']:'#F05456', draw_pathway.cfDict['ADP']:'#F05456',\
                    draw_pathway.cfDict['NAD+']:'#149B76', draw_pathway.cfDict['NADH']:'#149B76',\
                    draw_pathway.cfDict['NADPH']:'#2393CB', draw_pathway.cfDict['NADP+']:'#2393CB',draw_pathway.cfDict['NH4(+)']:'#F8C674'\
                    }


    colorMapping = colorConfig['colorMapping']
    for c in draw_pathway.cfDict.values():
        if c not in colorMapping:
            colorMapping[c] = colorConfig['OTHER_COFACTOR_COLOR']

    global_styles = {
        'graph': {
            'fontsize': '30',
            'fontname': 'Helvetica',
            'bgcolor': colorConfig['BACKGROUND_COLOR'],
        },
        'nodes': {
            'fontname': 'Helvetica',
            'fontsize': '30',
            'fontcolor': colorConfig['ALL_FONT_COLOR'],
        },
        'edges': {
            'fontname': 'Helvetica',
            'fontsize': '14',
        }
    }
    return global_styles, colorMapping

def _apply_styles(graph, styles):
    graph.graph_attr.update(
        ('graph' in styles and styles['graph']) or {}
    )
    graph.node_attr.update(
        ('nodes' in styles and styles['nodes']) or {}
    )
    graph.edge_attr.update(
        ('edges' in styles and styles['edges']) or {}
    )
    return graph


def draw_pathway(Pathway, imageFileName=None, imageFormat='pdf',
                 scaleLineWidth=False, scalingFactor=200.0,
                 cleanup=True, engine='dot', darkBackgroundMode=True):

    """
    Pathway: pathway object
    cleanup: delete the ".dot" file after drawing
    dark_background_mode: change all color settings to make graph for dark background
    """
    if darkBackgroundMode:
        colorConfig = color_configs['dark']
    else:
        colorConfig = color_configs['light']

    global_styles, colorMapping = _load_global_styles(colorConfig)

    g = Digraph('G', format=imageFormat, engine=engine)
    #Layout engines: circo dot fdp neato nop1 nop2 osage patchwork sfdp twopi

    g.graph_attr['size'] ="12,8"
    if imageFormat == 'png':
        g.graph_attr['dpi'] ='300'
    elif imageFormat == 'svg':
        g.graph_attr['dpi'] ='72'

    g = _apply_styles(g, global_styles)
    r_counter = 1

    #automatically use scaling factor if flux > 10 (scaling factor is set to the nearest 100)
    all_f = [abs(f) for f in Pathway.fluxes]

    if max(all_f) > 10:
        scaleLineWidth = True
        scalingFactor = 10**(math.ceil(math.log(max(all_f),10)))

    for rxn in Pathway.reactions:

        if Pathway.dGr_lst:
            g.node(rxn.rid, shape='point', color=colorConfig['RXN_NODE_COLOR'], xlabel=rxn.rid+'; '+str(abs(rxn.flux))+'; '+str(rxn.dGr), \
                fontsize='50', fontcolor=colorConfig['REACTION_COLOR'])
        else:
            g.node(rxn.rid, shape='point', color=colorConfig['RXN_NODE_COLOR'], xlabel=rxn.rid+'; '+str(abs(rxn.flux)), \
                fontsize='50', fontcolor=colorConfig['REACTION_COLOR'])

        if scaleLineWidth:
            lineW = '%i'% (10 * abs(rxn.flux)/scalingFactor + 1)
        elif rxn.flux >= 2 or rxn.flux <= -2:
            lineW = '%s'%(abs(rxn.flux)*2)
        else:
            lineW = '3'

        for met in rxn.reactants:
            if met in draw_pathway.cfDict.values():
                cf = abs(rxn.metabolites[met])
                if cf > 1:
                    clabel = '%i %s'%(cf, draw_pathway.met_dict[met]['name']) #show stoichiometric coefficient only when it is >= 2
                else:
                    clabel = draw_pathway.met_dict[met]['name']
                g.node(met+'_'+str(r_counter), shape=colorConfig['COFACTOR_SHAPE'], color=colorMapping[met], style="filled", \
                    fontsize='40', label=clabel)
                g.edge(met+'_'+str(r_counter), rxn.rid, penwidth=lineW, weight='1', arrowhead="none", color=colorConfig['EDGE_COLOR'])
            else:
                g.node(met, shape=colorConfig['NONCOFACTOR_SHAPE'], label=draw_pathway.met_dict[met]['name'], fontname='helvetica bold', \
                    fontsize='60', style="filled", color=colorConfig['NONCOFACTOR_COLOR'])
                g.edge(met, rxn.rid, penwidth=lineW, weight='2', arrowhead="none", color=colorConfig['EDGE_COLOR'])

        for met in rxn.products:
            if met in draw_pathway.cfDict.values():
                cf = abs(rxn.metabolites[met])
                if cf > 1:
                    clabel = '%i %s'%(cf, draw_pathway.met_dict[met]['name'])
                else:
                    clabel = draw_pathway.met_dict[met]['name']
                g.node(met+'_'+str(r_counter), shape=colorConfig['COFACTOR_SHAPE'], color=colorMapping[met], style="filled", \
                    fontsize='40', label=clabel)
                g.edge(rxn.rid, met+'_'+str(r_counter), weight='1', penwidth=lineW, color=colorConfig['EDGE_COLOR'])
            else:
                g.node(met, shape=colorConfig['NONCOFACTOR_SHAPE'], label=draw_pathway.met_dict[met]['name'], fontname='helvetica bold', \
                    fontsize='60', style="filled", color=colorConfig['NONCOFACTOR_COLOR'])
                g.edge(rxn.rid, met, penwidth=lineW, weight='2', color=colorConfig['EDGE_COLOR'])
        r_counter += 1

    if imageFileName is None:
        imageFileName = Pathway.name
    g.render(imageFileName, cleanup=cleanup)

#-------------------------------------------#




def read_metadata(metadata_file):
    metadata = {}
    with open(metadata_file, 'r') as inFile:
        for line in inFile:
            line = line.split()
            metadata[line[0]] = line[1]
    
    return metadata



def _single_pfba(model, obj_id):

    with model:
        try:
            model.objective = obj_id
        except:
            raise ValueError(f'{obj_id} is not present in model')
        
        pfba_sol = pfba(model)

    return pfba_sol.fluxes


# Calculate theoretical yields
def calc_flux_distributions(model, project_metadata):

    solution_df = pandas.DataFrame(columns=['equation'])
    for rxn in model.reactions: 
        solution_df.loc[rxn.id] = ({'equation' : rxn.reaction})

    test_rxns = ['ATPM', 'EX_co2_e', 'Biomass', 
                 'EX_'+project_metadata['substrate'], 'EX_'+project_metadata['product']]
    
    for rxn in test_rxns:
        if rxn != 'none':
            try:
                solution_df['Max_'+rxn] = _single_pfba(model, rxn)
            except:
                continue
    
    out_file = f"{project_metadata['project_id']}_solution_fluxes.tsv"
    solution_df.to_csv(out_file, sep='\t', index=False)
    print('Flux distribution results saved to:', out_file)
    model.objective = project_metadata['biomass_rxn']

    return solution_df


# Find exchange reaction ID for a metabolite
def _find_exchange(model, metID):

    for rxn in model.metabolites.get_by_id(metID).reactions:
        if len(rxn.metabolites) == 1:
            return rxn.id


# Yield of product on substrate
def calc_product_mty(model, project_metadata):

    product_exch = _find_exchange(model, project_metadata['product'])
    substrate_exch = _find_exchange(model, project_metadata['substrate'])

    model.objective = product_exch
    product_fluxes = pfba(model).fluxes

    mw_product = model.metabolites.get_by_id(project_metadata['product']).formula_weight
    mw_substrate = model.metabolites.get_by_id(project_metadata['substrate']).formula_weight

    yield_mol = -product_fluxes[product_exch] / product_fluxes[substrate_exch]
    yield_g = (yield_mol * mw_product) / mw_substrate # Mass yields

    yeild_str1 = '''Theoretical yield of {product_name} ({product}) on {substrate_name} ({substrate}) is {mols} mol/mol, or {gs} g/g'''.format(
        product_name=model.metabolites.get_by_id(project_metadata['product']).name, product=project_metadata['product'], 
        substrate_name=model.metabolites.get_by_id(project_metadata['substrate']).name, substrate=project_metadata['substrate'], 
        mols=round(yield_mol,3), gs=round(yield_g,3))
    print(yeild_str1)

    model.objective = project_metadata['biomass_rxn']
    product_fluxes = pfba(model).fluxes
    yield_mol = -product_fluxes[project_metadata['biomass_rxn']] / product_fluxes[substrate_exch]
    yield_g = (yield_mol * 1000.) / 180. # Mass yields

    yeild_str2 = '''Theoretical yield of Biomass on {substrate_name} ({substrate}) is {mols} mol/mol, or {gs} g/g
    '''.format(
        substrate_name=model.metabolites.get_by_id(project_metadata['substrate']).name, substrate=project_metadata['substrate'], 
        mols=round(yield_mol,3), gs=round(yield_g,3))
    print(yeild_str2)


# Plot production envelope
def prod_env_fig(model, project_metadata, file_type='png'):

    product_exch = _find_exchange(model, project_metadata['product'])
    substrate_exch = _find_exchange(model, project_metadata['substrate'])
    model.objective = product_exch

    substrate_uptake_rate = abs(model.reactions.get_by_id(substrate_exch).lower_bound)
    prod_env = production_envelope(model,[project_metadata['biomass_rxn']], carbon_sources=[substrate_exch])
    mw_substrate = model.metabolites.get_by_id(project_metadata['substrate']).formula_weight
    prod_env['biomass_yield_gg'] = (prod_env[project_metadata['biomass_rxn']] * 1000.) / (substrate_uptake_rate * mw_substrate)

    matplot.style.use('classic')
    new_fig, ax = matplot.subplots()
    product_label = model.metabolites.get_by_id(project_metadata['product']).name

    ax.plot(prod_env['biomass_yield_gg'], prod_env['mass_yield_maximum'], color='red')
    ax.set_xlabel('Biomass Yield (g/g)', fontsize=15)
    ax.set_ylabel(f'Product Yield (g/g)', fontsize=15)
    ax.set_title(f'{product_label} Production Envelope', fontsize=18)
    ax.fill_between(prod_env['biomass_yield_gg'], prod_env['mass_yield_maximum'], 0, facecolor='red', alpha=0.1)

    plot_file = f"{project_metadata['product']}_prodenv.{file_type}"
    new_fig.savefig(plot_file, dpi=300, format=file_type)
    print('Production envelope figure saved to', plot_file)
    model.objective = project_metadata['biomass_rxn']

    return prod_env


# Create Escher map
def plot_escher_map(model, project_metadata):

    product_exch = _find_exchange(model, project_metadata['product'])
    model.objective = product_exch
    product_fluxes = pfba(model).fluxes

    builder = escher.Builder(map_json=project_metadata['escher_map'], model=model, height=700, identifiers_on_map='bigg_id', show_gene_reaction_rules=False)
    builder.reaction_data = product_fluxes
    builder.reaction_scale = [
        { 'type': 'value','value': 0.0, 'color': 'gray',  'size': 0.5 },
        { 'type': 'value','value': 0.1, 'color': 'black', 'size': 8 },
        { 'type': 'value','value': 0.5, 'color': 'black', 'size': 15 },
        { 'type': 'value','value': 2.0, 'color': 'black', 'size': 20 },
    ]

    escher_base = project_metadata['escher_map'].split('/')[-1].rstrip('.json').replace('_','')
    escher_html = f"{project_metadata['strain_id']}_{escher_base}_escher.html"
    builder.save_html(escher_html)
    print('Escher pathway html saved to', escher_html, '\n')
    model.objective = project_metadata['biomass_rxn']

    return builder



# Plot precursor acquisition pathway
def product_precursor_path(model, project_metadata, namespace='bigg', file_type='pdf'):

    product_exch = _find_exchange(model, project_metadata['product'])
    model.objective = product_exch
    product_solution = pfba(model)

    draw_pathway.met_dict = {met.id:{'name':met.name} for met in model.metabolites}
    draw_pathway.rxn_dict = {rxn.id:{'name':rxn.name, 'sji':{x.id:y for x, y in rxn.metabolites.items()}} for rxn in model.reactions}
    if namespace == 'bigg':
        draw_pathway.cfDict = {'H2O': 'h2o_c', 'ATP': 'atp_c', 'NAD+': 'nad_c', 'NADH': 'nadh_c', 'NADPH': 'nadph_c', 'NADP+': 'nadp_c', 
                'O2': 'o2_c', 'ADP': 'adp_c', 'Orthophosphate': 'pi_c', 'CoA': 'coa_c', 'CO2': 'co2_c', 'Diphosphate': 'ppi_c', 
                'UDP': 'udp_c', 'FAD': 'fad_c', 'S-Adenosyl-L-methionine': 'amet_c', 'AMP': 'amp_c', 'GDP': 'gdp_c', 'GTP': 'gtp_c', 
                'CTP': 'ctp_c', 'UTP': 'utp_c', 'H+': 'h_c', 'ITP': 'itp_c', 'IDP': 'idp_c', 'UMP': 'ump_c', 'CDP': 'cdp_c', 
                'GMP': 'gmp_c', 'CMP': 'cmp_c', 'dATP': 'datp_c', 'dADP': 'datdp_c', 'Ubiquinone (Mitochondria)': 'q6_m', 
                'Ubiquinol (Mitochondria)': 'q6h2_m', 'IMP': 'imp_c', 'dUTP': 'dutp_c', 'NH4(+)': 'nh4_c', 'H2O (Mitochondria)': 'h2o_m', 
                'NAD+ (Mitochondria)': 'nad_m', 'NADH (Mitochondria)': 'nadh_m', 'NADPH (Mitochondria)': 'nadph_m', 
                'NADP+ (Mitochondria)': 'nadp_m', 'O2 (Mitochondria)': 'o2_m', 'ATP (Mitochondria)': 'atp_m', 'ADP (Mitochondria)': 'adp_m', 
                'Orthophosphate (Mitochondria)': 'pi_m', 'CoA (Mitochondria)': 'coa_m', 'CO2 (Mitochondria)': 'co2_m', 
                'Diphosphate (Mitochondria)': 'ppi_m', 'H+ (Mitochondria)': 'h_m', 'NH4(+) (Mitochondria)': 'nh4_m', 'FADH2 (Mitochondria)': 'fadh2_m', 
                'FAD (Mitochondria)': 'fad_m', 'H2O (Peroxisome)': 'h2o_x', 'O2 (Peroxisome)': 'o2_x', 'NAD+ (Peroxisome)': 'nad_x', 
                'NADH (Peroxisome)': 'nadh_x', 'NADP+ (Peroxisome)': 'nadp_x', 'NADPH+ (Peroxisome)': 'nadph_x', 'CoA (Peroxisome)': 'coa_x', 
                'H+ (Peroxisome)': 'h_x'}
    elif namespace == 'wei':
        draw_pathway.cfDict = {'H2O': 'M_m32', 'ATP': 'M_m141', 'NAD+': 'M_m122', 'NADH': 'M_m123', 'NADPH': 'M_m40', 'NADP+': 'M_m41', 
              'O2': 'M_m109', 'ADP': 'M_m143', 'Orthophosphate': 'M_m35', 'CoA': 'M_m69', 'CO2': 'M_m82', 'Diphosphate': 'M_m203', 
              'UDP': 'M_m11', 'FAD': 'M_m528', 'S-Adenosyl-L-methionine': 'M_m256', 'AMP': 'M_m86', 'GDP': 'M_m268', 'GTP': 'M_m266', 
              'CTP': 'M_m406', 'UTP': 'M_m439', 'H+': 'M_m10', 'ITP': 'M_m843', 'IDP': 'M_m839', 'UMP': 'M_m149', 'CDP': 'M_m500', 
              'GMP': 'M_m93', 'CMP': 'M_m95', 'dATP': 'M_m606', 'dADP': 'M_m607', 'Ubiquinone (Mitochondria)': 'M_m468', 
              'Ubiquinol (Mitochondria)': 'M_m471', 'IMP': 'M_m91', 'dUTP': 'M_m462', 'NH4(+)': 'M_m38', 'H2O (Mitochondria)': 'M_m26', 
              'NAD+ (Mitochondria)': 'M_m27', 'NADH (Mitochondria)': 'M_m30', 'NADPH (Mitochondria)': 'M_m176', 
              'NADP+ (Mitochondria)': 'M_m178', 'O2 (Mitochondria)': 'M_m64', 'ATP (Mitochondria)': 'M_m46', 'ADP (Mitochondria)': 'M_m197', 
              'Orthophosphate (Mitochondria)': 'M_m58', 'CoA (Mitochondria)': 'M_m73', 'CO2 (Mitochondria)': 'M_m44', 
              'Diphosphate (Mitochondria)': 'M_m204', 'H+ (Mitochondria)': 'M_m28', 'NH4(+) (Mitochondria)': 'M_m581', 'FADH2 (Mitochondria)': 'M_m566', 
              'FAD (Mitochondria)': 'M_m563', 'H2O (Peroxisome)': 'M_m375', 'O2 (Peroxisome)': 'M_m221', 'NAD+ (Peroxisome)': 'M_m116', 
              'NADH (Peroxisome)': 'M_m119', 'NADP+ (Peroxisome)': 'M_m740', 'NADPH+ (Peroxisome)': 'M_m742', 'CoA (Peroxisome)': 'M_m182', 
              'H+ (Peroxisome)': 'M_m118'}
    
    flux_df = product_solution[(product_solution.fluxes > 1e-3) | (product_solution.fluxes < -1e-3)]
    max_product_pathway = Pathway(reaction_ids=list(flux_df.index.values), fluxes=list(flux_df.round(3)))

    draw_pathway(max_product_pathway, imageFileName=f"{project_metadata['product']}_pathway", darkBackgroundMode=False, imageFormat=file_type)
    print('Substrate precursor pathway plot saved to:', f"{project_metadata['strain_id']}_{project_metadata['product']}_pathway.{file_type}")
    model.objective = project_metadata['biomass_rxn']
